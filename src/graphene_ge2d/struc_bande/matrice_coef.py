#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Génération des matrices à diagonaliser
"""
import numpy as np
import scipy.special as sp
from numba import jit
import numba_scipy  # noqa: F401
from . import _constantes as cons  # c2, eta, v_0, r_0_bar, eps, z_bar, b1, b2, nc


#@jit(nopython=True)
def coeffs(G):
    """
    Calcule le terme qui change avec le signe de G
    G: Les vecteurs G qu'on veut le coefficient
    """
    # La somme sur i
    return 1 #+ np.exp(-1j * 2 * np.pi * np.dot(G, cons.c2) / cons.eta)


#@jit(nopython=True)
def coeffs_sim(G):
    """
    Calcule la partie du coefficient de fourrier qui est symétrique
    dans la matrice
    G: Les vecteurs G qu'on veut le coefficient
    """
    norm_g = np.linalg.norm(G)
    if norm_g <= 1e-100:
        return cons.v_0 * cons.r_0_bar / (cons.eps * cons.eta)
    bessel_boi = sp.j1(2 * np.pi * norm_g * cons.r_0_bar / cons.eta) / norm_g
    petit_exp = np.exp(- 2 * np.pi * norm_g * cons.z_bar / cons.eta)
    return cons.v_0 * (cons.r_0_bar / (cons.eps * cons.eta)) * bessel_boi * petit_exp


#@jit(nopython=True, parallel=True)
def get_matrice(k, ndiag, vecteurs):
    """
    Get la matrice à diagonaliser
    n: matrice générée à partir de l'ordre n (n doit être n = 4m^2, m entier,
                                              sinon arrondi)
    k: un vecteur de la zone de Brillouin
    z: hauteur du potentiel appliqué
    ndiag: nombre de diagonales voulues (Le plus près possible de la diagonale
                                         centrale)
    """
    taille_matrice = cons.n_mat
    matrix = np.zeros((taille_matrice, taille_matrice), dtype=np.complex128)
    for i, curr_g in enumerate(vecteurs):
        matrix[i, i] = np.dot(k + curr_g, k + curr_g) / cons.nc
        for j in range(taille_matrice):
            if j > i < len(vecteurs):
                other_g = vecteurs[j]
                simetric_coef = coeffs_sim(curr_g - other_g)
                changing_coef = coeffs(curr_g - other_g)
                matrix[i, j] = - (simetric_coef * changing_coef)
    matrix = matrix + np.conjugate(matrix.transpose())
    return matrix
