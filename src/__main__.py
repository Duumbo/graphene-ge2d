#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Diag de la matrice pour la structure de bande
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.linalg as linalg
from scipy import sparse
import scipy.special as sp
from numba import jit, types
import numba_scipy

# Constantes
SIZE_MATRICES = 10000
PATH = [100, 50, 100]  # Nombres de points dans le parcours de prise de données. En ordre selon la partie
a = 130  # nm
nc = 2
a1 = np.array([a / 2, - a * np.sqrt(3) / 2])
a2 = np.array([a, 0])
c1 = np.array([0, 0])
c2 = np.array([0, -1 / np.sqrt(3)])
eta = np.sqrt(nc)
eps = np.sqrt(3)/2
v_0 = 1  # eV
r_0 = 30  # nm
r_0_bar = r_0 / a
distance_plaque = 45  # nm
b1 = np.array([0, - 4 * np.pi / (3)])
b2 = np.array([(2 * np.pi) / (np.sqrt(3)), (2 * np.pi) / (3)])
e_0_inv = a * a / 22_450
# Points importants de la zone de Brillouin
# Le point M
POINT_M = np.array([np.pi, np.pi / (np.sqrt(3))])
# Le point K
POINT_K =  np.array([4 * np.pi / 3, 0])
# Le point Gamma
GAMMA = np.array([0, 0])


@jit(nopython=True)
def coeffs(G, z):
    """
    Calcule le terme qui change avec le signe de G
    G: Les vecteurs G qu'on veut le coefficient
    z: La distance du 2DEG au potentiel appliqué
    """
    # La somme sur i
    return 1 + np.exp(-1j * 2 * np.pi * np.dot(G, c2) / eta)


@jit(nopython=True)
def coeffs_sim(G, z):
    """
    Calcule la partie du coefficient de fourrier qui est symétrique dans la matrice
    G: Les vecteurs G qu'on veut le coefficient
    z: La distance au 2DEG qu'on applique le potentiel
    """
    z_bar = z / a
    norm_g = np.linalg.norm(G)
    if norm_g <= 1e-5:
        return v_0 * r_0_bar / (eps * eta)
    bessel_boi = sp.j1(2 * np .pi * norm_g * r_0_bar / eta) / norm_g
    petit_exp = np.exp(- 2 * np.pi * norm_g * z_bar / eta)
    return v_0 * (r_0_bar / (eps * eta)) * bessel_boi * petit_exp


@jit(nopython=True)
def gen_gs(number):
    """
    Génère des vecteurs G de l'espace réciproque
    number : 4*nombre*nombre de vecteur voulu
    """
    # Génère 4*number*number vecteurs
    output = []
    for i in range(-number, number):
        current_g1 = i * b1
        for j in range(-number, number):
            output.append(current_g1 + (j * b2))
    return output


@jit(nopython=True)
def get_matrice(n, k, z):
    """
    Get la matrice à diagonaliser
    n: matrice n x n (n doit être n = 4m^2, m entier, sinon arrondi)
    k: un vecteur de la zone de Brillouin
    z: hauteur du potentiel appliqué
    """
    m = int(np.sqrt(n / 4))
    vecteurs = gen_gs(m)
    matrix = np.empty((m*m*4, 3), dtype = np.complex128)
    for i, curr_g in enumerate(vecteurs):
        matrix[i, 0] = np.dot(k + curr_g, k + curr_g) / nc
        if i < len(vecteurs) - 1:
            other_g = vecteurs[i + 1]
            simetric_coef = coeffs_sim(curr_g - other_g, z)
            changing_coef = coeffs(curr_g - other_g, z)
            matrix[i, 1] = ( take_the_conjugate := - (simetric_coef * changing_coef))
            matrix[i, 2] = take_the_conjugate.real - take_the_conjugate.imag
    return matrix


def joli_graphique(fig, ax, *args):
    """
    Mise en forme du graphique
    fig: l'objet figure
    ax: l'axe actif
    """
    chemin_norm = args[0]
    position_m = chemin_norm[PATH[0] - 1]
    position_k = chemin_norm[PATH[0] + PATH[1] - 1]
    ax.set_xticks([0., position_m, position_k, chemin_norm[-1]], [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"])
    return fig, ax


def __main__():
    fig, ax = plt.subplots()

    # Chemin 1, c'est de Gamma à M
    partie_chemin_1 = np.linspace(GAMMA, POINT_M, PATH[0], endpoint = False)
    pas_1 = np.linalg.norm(partie_chemin_1[0, :] - partie_chemin_1[1, :])
    chemin_norm_1 = np.arange(PATH[0] * pas_1, step = pas_1)

    # Chemin 2, c'est de M à K
    partie_chemin_2 = np.linspace(POINT_M, POINT_K, PATH[1], endpoint = False)
    pas_2 = np.linalg.norm(partie_chemin_2[0, :] - partie_chemin_2[1, :])
    chemin_norm_2 = np.arange(PATH[0] * pas_1, (PATH[0] * pas_1) + PATH[1] * pas_2, step = pas_2)

    # Chemin 3, c'est de K à Gamma
    partie_chemin_3 = np.linspace(POINT_K, GAMMA, PATH[2], endpoint = False)
    pas_3 = np.linalg.norm(partie_chemin_3[0, :] - partie_chemin_3[1, :])
    chemin_norm_3 = np.arange((PATH[0] * pas_1) + (PATH[1] * pas_2), (PATH[0] * pas_1) + (PATH[1] * pas_2) + (PATH[2] * pas_3), step = pas_3)

    # On pondère l'axe des x du graphique par la variation de la norme
    chemin_norm = np.append(chemin_norm_1, chemin_norm_2)
    chemin_norm = np.append(chemin_norm, chemin_norm_3)

    # On combine les chemins entre-eux
    chemin = np.append(partie_chemin_1, partie_chemin_2, axis = 0)
    chemin = np.append(chemin, partie_chemin_3, axis = 0)
    # Tableau où l'on place les valeurs propres obtenues
    list_eigvals = np.empty((chemin.shape[0], 10))
    try:
        for i, k in tqdm(enumerate(chemin)):
            diagonals = get_matrice(SIZE_MATRICES, k, distance_plaque)
            matrix_to_diag = sparse.diags([diagonals[:, 0], diagonals[:-1, 1], diagonals[:-1, 2]], [0, 1, -1])
            if False:
                """
                Debug, pour afficher la heat map des valeurs absolues des matrices.
                """
                np.savetxt(f"{i}.txt", matrix_to_diag)
                sns.heatmap(np.abs(matrix_to_diag), ax = ax)
                plt.show()
                fig.savefig("heatmap_matrice_F.png")

            eigval = sparse.linalg.eigsh(matrix_to_diag, k = 10, return_eigenvectors = False)
            if np.any(eigval > 1e90):
                print(i)
                break
            list_eigvals[i, :] = eigval
        for i in range(10):
            #print(list_eigvals)
            ax.plot(chemin_norm, list_eigvals[:, i], "o", markersize=2)

            # Beauté du graphique
            fig, ax = joli_graphique(fig, ax, *[chemin_norm])

        plt.show()
    except Exception as e:
        print("Error, safety save")
        np.savetxt("safety_save.txt.gz", list_eigvals)
        print(e)


if __name__ == "__main__":
    __main__()

