use scilib::math::bessel;
use rgsl::exponential;
use num_complex::Complex;
use lapack::zhpevd;

const V_0: f64 = - 100.0;
const R_0_BAR: f64 = 10.0 / 100.0;
const ETA: f64 = 1.0;
const N_COQ: u32 = 15;
const PI: f64 = std::f64::consts::PI;
const Z_BAR: f64 = 40.0 / 100.0;

fn coquilles(n_coq: u32) -> (Vec<[f64; 2]>, Vec<[i64; 2]>) {
    let EPS = f64::sqrt(3.0) / 2.0;
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

fn potentiel(g: [f64; 2]) -> Complex::<f64> {
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
    let EPS = f64::sqrt(3.0) / 2.0;
    let (vecteurs_g, indices_g) = coquilles(N_COQ);
    let n_vect = indices_g.len();
    let deff_val = Complex::from(V_0 * R_0_BAR / (EPS * ETA));

    // Let's get the potential for all those gs
    let mut pot: Vec<Complex::<f64>> = Vec::new();
    for i in 0..n_vect {
        pot.push(potentiel(vecteurs_g[i]));
    }

    // Let's gen the potential for all non 0 Gs.
    let mut matrice: Vec<Complex::<f64>> = Vec::new();
    for i in 0..n_vect {
        let g1 = indices_g[i];
        matrice.push(Complex::from(vecteurs_g[i][0]*vecteurs_g[i][0] + vecteurs_g[i][1]*vecteurs_g[i][1]));
        for j in (i + 1)..n_vect {
            let g2 = indices_g[j];
            for t in 0..n_vect {
                if (g1[0] - g2[0] == indices_g[t][0]) && (g1[1] - g2[1] == indices_g[t][1]) {
                    matrice.push(pot[t]);
                    break;
                }
                if t == n_vect - 1 {
                    matrice.push(deff_val);
                }
            }
        }
    }
    println!("{:?}", matrice);
    // First let's allocate the working space for LAPACK
    let mut work: Vec<Complex::<f64>> = Vec::with_capacity(n_vect);
    let mut rwork: Vec<f64> = Vec::with_capacity(n_vect);
    let mut iwork: Vec<i32> = Vec::with_capacity(n_vect);
    let ldz: i32 = n_vect as i32;
    let lwork: i32 = n_vect as i32;
    let lrwork: i32 = n_vect as i32;
    let liwork: i32 = n_vect as i32;
    let mut info: i32 = 0;

    let mut eigen_val: Vec<f64> = Vec::with_capacity(n_vect);
    let mut zvec: Vec<Complex::<f64>> = Vec::with_capacity(n_vect);

    // Now let's diagonalise it
    unsafe {
        zhpevd(b"N"[0], b"U"[0], n_vect as i32, &mut matrice, &mut eigen_val, &mut zvec, ldz,
                &mut work, lwork, &mut rwork, lrwork, &mut iwork, liwork, &mut info);
    }
}
