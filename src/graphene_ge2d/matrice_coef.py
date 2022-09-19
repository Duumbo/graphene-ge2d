#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Génération des matrices à diagonaliser
"""
import numpy as np
import scipy.special as sp
from numba import jit, types
import numba_scipy
from ._constantes import *


@jit(nopython=True)
def coeffs(G):
    """
    Calcule le terme qui change avec le signe de G
    G: Les vecteurs G qu'on veut le coefficient
    """
    # La somme sur i
    return 1 + np.exp(-1j * 2 * np.pi * np.dot(G, c2) / eta)


@jit(nopython=True)
def coeffs_sim(G):
    """
    Calcule la partie du coefficient de fourrier qui est symétrique dans la matrice
    G: Les vecteurs G qu'on veut le coefficient
    """
    norm_g = np.linalg.norm(G)
    if norm_g <= 1e-100:
        return v_0 * r_0_bar / (eps * eta)
    bessel_boi = sp.j1(2 * np.pi * norm_g * r_0_bar / eta) / norm_g
    petit_exp = np.exp(- 2 * np.pi * norm_g * z_bar / eta)
    return v_0 * (r_0_bar / (eps * eta)) * bessel_boi * petit_exp


@jit(nopython=True, parallel=True)
def gen_gs(ordre):
    """
    Génère des vecteurs G de l'espace réciproque
    order : ordre du réseau (ordre-voisins) de vecteur voulu
    """
    # Génère 4*number*number vecteurs
    output = []
    n_possibles_1 = range(-ordre, ordre + 1)
    n_possibles_2 = range(-ordre, ordre + 1)
    n_possibles = [[i, j] for i in n_possibles_1 for j in n_possibles_2]
    n_poss = []
    [n_poss.append(x) for x in n_possibles if x not in n_poss]
    for n in n_poss:
        output.append(n[0] * b1 + n[1] * b2)
    return output, len(n_poss)


@jit(nopython=True, parallel=True)
def get_matrice(taille_matrice, k, ndiag, vecteurs):
    """
    Get la matrice à diagonaliser
    n: matrice générée à partir de l'ordre n (n doit être n = 4m^2, m entier, sinon arrondi)
    k: un vecteur de la zone de Brillouin
    z: hauteur du potentiel appliqué
    ndiag: nombre de diagonales voulues (Le plus près possible de la diagonale centrale)
    """
    matrix = np.empty((taille_matrice, ndiag), dtype = np.complex128)
    for i, curr_g in enumerate(vecteurs):
        matrix[i, 0] = np.dot(k + curr_g, k + curr_g) / nc
        for j in range((ndiag - 1) // 2):
            if i < len(vecteurs) - (1 + j):
                other_g = vecteurs[i + j + 1]
                simetric_coef = coeffs_sim(curr_g - other_g)
                changing_coef = coeffs(curr_g - other_g)
                matrix[i, 2 * j + 1] = ( take_the_conjugate := - (simetric_coef * changing_coef))
                matrix[i, 2 * j + 2] = np.conjugate(take_the_conjugate)
    return matrix * e_0_inv

