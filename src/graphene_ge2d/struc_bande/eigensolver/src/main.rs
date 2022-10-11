extern crate lapack_src;
extern crate lapack_sys;

use std::fs;
use std::io::Write;
use std::path::Path;
use scilib::math::bessel;
use rgsl::exponential;
use num_complex::Complex;
use lapack::zhpevd;

const PI: f64 = std::f64::consts::PI;

fn coquilles(n_coq: u32) -> (Vec<[f64; 2]>, Vec<[i64; 2]>) {
    let g_max_sq: f64 = (4.0 / 3.0) * (n_coq as f64);

    let mut list_perm: Vec<[f64; 2]> = Vec::new();
    let mut list_ind: Vec<[i64; 2]> = Vec::new();
    let INV_SQRT3 = 1.0 / f64::sqrt(3.0);
    let INV_2SQRT3 = 2.0 / f64::sqrt(3.0);
    for i in (- (n_coq as i64))..((n_coq + 1) as i64) {
        for j in (- (n_coq as i64))..((n_coq + 1) as i64) {
            let x = (j as f64 * INV_SQRT3) + (i as f64 * INV_2SQRT3);
            if x*x + ((j*j) as f64) < g_max_sq {
                list_perm.push([j as f64, x]);
                list_ind.push([i, j]);
            }
        }
    }
    (list_perm, list_ind)
}

fn potentiel(g: [f64; 2],
             EPS: f64,
             R_0_BAR: f64,
             ETA: f64,
             Z_BAR: f64,
             V_0: f64) -> Complex::<f64> {
    let EPS = f64::sqrt(3.0) / 2.0;
    let norm_g = f64::sqrt(g[0]*g[0] + g[1]*g[1]);
    if norm_g < 1e-50 {
        return Complex::from(V_0 * R_0_BAR / (EPS * ETA));
    }
    let bessi = bessel::j(2.0 * PI * norm_g * R_0_BAR / ETA, 1) / norm_g;
    let expo = exponential::exp(- 2.0 * PI * norm_g * Z_BAR / ETA);
    return V_0 * R_0_BAR * bessi * expo / (EPS * ETA);
}

fn main() {
    let path_parametres = Path::new("raw_data/parametres.txt");
    let chemin = fs::read_to_string(path_parametres)
        .expect("Data for parameters does not exist, check that I look at the right place.");
    let params = chemin.split(" ").collect::<Vec<&str>>();
    let V_0: f64 = params[0].parse().expect("V_0 must be a real float.");
    let R_0_BAR: f64 = params[1].parse().expect("R_0_BAR must be a real float.");
    let ETA: f64 = params[2].parse().expect("ETA must be a real float.");
    let N_COQ: u32 = params[3].parse().expect("N_COQ must be a real INT.");
    let Z_BAR: f64 = params[4].parse().expect("Z_BAR must be a real float.");
    let E_0: f64 = params[5].parse().expect("E_0 must be a real float.");
    let _EPS: f64 = params[6].parse().expect("EPS must be a real float.");
    let potentiel_param_set = |g: [f64; 2]| potentiel(g, _EPS, R_0_BAR, ETA, Z_BAR, V_0);

    let (vecteurs_g, indices_g) = coquilles(N_COQ);
    let n_vect = indices_g.len();
    let deff_val = Complex::from(0.0);

    // Let's get the potential for all those gs
    let mut pot: Vec<Complex::<f64>> = Vec::new();
    for i in 0..n_vect {
        pot.push(potentiel_param_set(vecteurs_g[i]));
        //pot.push(Complex::from(0.0));
    }

    // Let's gen the potential for all non 0 Gs.
    let mut matrice: Vec<Complex::<f64>> = vec![Complex::from(0.0); n_vect * (n_vect + 1) / 2];
    for i in 0..n_vect {
        let g1 = indices_g[i];
        matrice[i * (i + 3)/2] = Complex::from(vecteurs_g[i][0]*vecteurs_g[i][0] + vecteurs_g[i][1]*vecteurs_g[i][1]);
        for j in 0..i {
            let g2 = indices_g[j];
            for t in 0..n_vect {
                if (g2[0] - g1[0] == indices_g[t][0]) && (g2[1] - g1[1] == indices_g[t][1]) {
                   matrice[j + (i*(i as i32 + 1)as usize)/2] = - pot[t] / Complex::from(E_0);
                  break;
                }
                //if t == n_vect - 1 {
                //    matrice[j + (i*(i as i32 + 1) as usize)/2] = - potentiel(
                //        [vecteurs_g[j][0] - vecteurs_g[i][0], vecteurs_g[j][1] - vecteurs_g[i][1]]) / Complex::from(E_0);
                    //matrice.push(deff_val);
                //}
            }
        }
    }
    //println!("{:?}", matrice);
    // First let's allocate the working space for LAPACK
    let ldz: i32 = n_vect as i32;
    let lwork: i32 = 2 * n_vect as i32;
    let lrwork: i32 = 1 + 5 * n_vect as i32 + 2 * (n_vect*n_vect) as i32;
    let liwork: i32 = 3 + 5 * n_vect as i32;
    let mut work: Vec<Complex::<f64>> = vec![Complex::from(0.0); lwork.try_into().unwrap()];
    let mut rwork: Vec<f64> = vec![0.0; lrwork.try_into().unwrap()];
    let mut iwork: Vec<i32> = vec![0; liwork.try_into().unwrap()];
    let mut info: i32 = 0;

    let mut eigen_val: Vec<f64> = vec![0.0; n_vect];
    let mut zvec: Vec<Complex::<f64>> = vec![Complex::from(0.0); n_vect];

    // Now let's diagonalise it
    // We need to take points k on a path to modify the diagonal.
    let path_chemin = Path::new("raw_data/chemin.txt");
    let chemin = fs::read_to_string(path_chemin)
        .expect("Data for path in the Brillouin zone does not exist! Check if it is at the right place.");
    let mut output_eigval = fs::File::create("raw_data/list_eigvals.txt").expect("Couldn't create output_eigval");
    for k in chemin.split("\n") {
        let next_k: Vec<&str> = k.split(" ").collect::<Vec<&str>>();
        if next_k[0].is_empty() { break; }  // Assume End Of File
        let next_k: [f64; 2] = [next_k[0].parse().expect("Invalid format in chemin.txt"),
                                next_k[1].parse().expect("Invalid format in chemin.txt")];
        diagonalise(&matrice, n_vect as i32, &mut eigen_val, &mut zvec, &ldz,
                    &mut work, &lwork, &mut rwork, &lrwork, &mut iwork, &liwork, &mut info);
        let mut strings: Vec<String> = Vec::new();
        for i in 0..n_vect {
            let curr_g = vecteurs_g[i];
            // The diagonal elements are always i away from the last one
            matrice[i*(i + 3) / 2] = Complex::from((curr_g[0] + next_k[0]).powf(2.0) + (curr_g[1] + next_k[1]).powf(2.0));
            strings.push((eigen_val[i] * E_0).to_string());
        }
        writeln!(output_eigval, "{}", strings.join(" "));
    }
}

fn diagonalise(
    matrice: &[Complex::<f64>],  // Reference to the matrix you want to diagonalise
    n_vect: i32,                 // Dimension of the matrix (n_vect*n_vect)
    eigen_val: &mut [f64],       // A mutable reference to the output vector of eigenvalues
    zvec: &mut [Complex::<f64>], // A mutable reference to the output matrix of eigenvectors
    ldz: &i32,                   // The dimension of zvec
    work: &mut [Complex::<f64>], // zhpevd working memory for all below
    lwork: &i32,
    rwork: &mut [f64],
    lrwork: &i32,
    iwork: &mut [i32],
    liwork: &i32,
    info: &mut i32               // An int representing if it converged or not
) {  // Function changes by-reference the vector eigen_val
    let mut vec_mat = (*matrice).to_vec();
    unsafe {
        zhpevd(b"N"[0], b"U"[0], n_vect as i32, &mut vec_mat, eigen_val, zvec, *ldz,
                work, *lwork, rwork, *lrwork, iwork, *liwork, info);
    }
}
