#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Test de la construction des matrices.
Tests unitaire comportant:
    Une matrice calculée explicitement à la main à comparer
"""
import numpy as np
from scipy.special import j1
import pytest
import graphene_ge2d.struc_bande.matrice_coef as matrix
from graphene_ge2d.struc_bande._constantes import r_0_bar, eta, eps, v_0, z_bar, e_0_inv, nc


def test_matrice_explicite():
    """
    Test de construction d'une matrice 7x7 en la comparant à celles codée
    avec k = (0, 0)
    """
    # Vecteurs de la base
    b_1 = np.array([0, - 2 / np.sqrt(3)])
    b_2 = np.array([1, 1 / np.sqrt(3)])
    c2 = np.array([0, - 1 / np.sqrt(3)])

    # Les G à utiliser dans l'ordre
    g_1 = np.array([0, 0])
    g_2 = b_2
    g_3 = - b_1
    g_4 = - b_2 - b_1
    g_5 = - b_2
    g_6 = b_1
    g_7 = b_2 + b_1
    g_8 = - b_1 + b_1
    g_9 = b_2 - b_1
    g_array = np.array([g_1, g_2, g_3, g_4, g_5, g_6, g_7, g_8, g_9])

    # Fonction des coefficients de fourier

    def coefficients_u(g):
        norm_g = np.linalg.norm(g)
        if norm_g <= 1e-100:
            return v_0 * r_o_bar / (eps * eta)
        bessel_boi = j1(2 * np.pi * norm_g * r_0_bar / eta) / norm_g
        exponentielle = np.exp(- 2 * np.pi * norm_g * z_bar / eta)
        somme = 1 + np.exp(- 1j * 2 * np.pi * np.dot(g, c2) / eta)
        return v_0 * r_0_bar * bessel_boi * exponentielle * somme / (eps * eta)

    diag_principale = e_0_inv * np.array([
            g_1[0]*g_1[0] + g_1[1]*g_1[1],
            g_2[0]*g_2[0] + g_2[1]*g_2[1],
            g_3[0]*g_3[0] + g_3[1]*g_3[1],
            g_4[0]*g_4[0] + g_4[1]*g_4[1],
            g_5[0]*g_5[0] + g_5[1]*g_5[1],
            g_6[0]*g_6[0] + g_6[1]*g_6[1],
            g_7[0]*g_7[0] + g_7[1]*g_7[1],
            g_8[0]*g_8[0] + g_8[1]*g_8[1],
            g_9[0]*g_9[0] + g_9[1]*g_9[1],
    ]) / nc
    test_diag = e_0_inv * np.array([
            np.dot(g_1, g_1),
            np.dot(g_2, g_2),
            np.dot(g_3, g_3),
            np.dot(g_4, g_4),
            np.dot(g_5, g_5),
            np.dot(g_6, g_6),
            np.dot(g_7, g_7),
            np.dot(g_8, g_8),
            np.dot(g_9, g_9),
    ]) / nc
    ##print(f"{diag_principale=}")
    #print(f"{test_diag=}")
    assert np.isclose(test_diag, diag_principale).all()
    diag_sup = e_0_inv * np.array([
            - coefficients_u(g_1 - g_2),
            - coefficients_u(g_2 - g_3),
            - coefficients_u(g_3 - g_4),
            - coefficients_u(g_4 - g_5),
            - coefficients_u(g_5 - g_6),
            - coefficients_u(g_6 - g_7),
            - coefficients_u(g_7 - g_8),
            - coefficients_u(g_8 - g_9),
    ])
    diag_inf = e_0_inv * np.array([
            - coefficients_u(- g_1 + g_2),
            - coefficients_u(- g_2 + g_3),
            - coefficients_u(- g_3 + g_4),
            - coefficients_u(- g_4 + g_5),
            - coefficients_u(- g_5 + g_6),
            - coefficients_u(- g_6 + g_7),
            - coefficients_u(- g_7 + g_8),
            - coefficients_u(- g_8 + g_9),
    ])
    #print(f"{g_array=}")
    diags = matrix.get_matrice(9, np.array([0, 0]), 3, g_array)
    #print(f"{diag_principale=}")
    print(f"{diag_sup=}")
    print(f"{diag_inf=}")
    #print(f"{diags[:, 0]=}")
    print(f"{diags[:, 1]=}")
    print(f"{diags[:, 2]=}")
    assert np.isclose(diag_principale, diags[:, 0]).all()
    assert np.isclose(diag_sup, diags[:-1, 1]).all()
    assert np.isclose(diag_inf, diags[:-1, 2]).all()


def test_norm_linalg():
    g = np.array([1+1j, 1+1j])
    #assert np.isclose(np.linalg.norm(g), )
    assert np.isclose(np.dot(g, g), 4j)

