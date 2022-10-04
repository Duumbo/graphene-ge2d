#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Définition des constantes
"""
import sys
import numpy as np
from math import factorial

cons = sys.modules[__name__]

# Constantes
cons.SIZE_MATRICES = 3
cons.n_mat = 37
#cons.n_mat = 225
# Nombres de points dans le parcours de prise de données.
# En ordre selon la partie
cons.DENS_POINT = 400
cons.longueur_tot = (1 / np.sqrt(3)) + (1 / 3) + (2 / 3)
cons.PATH = [int(DENS_POINT * (1 / np.sqrt(3)) / longueur_tot), int(DENS_POINT * (1 / 3) / longueur_tot), int(DENS_POINT * (2 / 3) / longueur_tot)]
# Nombre de diagonales au centre de la matrice creuse (doit être impair)
cons.NDIAG = cons.n_mat
cons.a = 130  # nm
cons.nc = 1
cons.a1 = np.array([1 / 2, - np.sqrt(3) / 2])
cons.a2 = np.array([1, 0])
cons.c1 = np.array([0, 0])
#cons.c2 = np.array([0, -1 / np.sqrt(3)])
cons.c2 = np.array([0, 0])
cons.eta = np.sqrt(nc)
cons.eps = np.sqrt(3)/2
cons.v_0 = -1  # eV
cons.r_0 = 30  # nm
cons.r_0_bar = r_0 / a
cons.distance_plaque = 45 # nm
cons.b1 = np.array([0, - 2 / np.sqrt(3)])
cons.b2 = np.array([1, 1 / np.sqrt(3)])
cons.e_0 = 22_450 / (a * a)
cons.e_0_inv = a * a / 22_450
cons.z_bar = distance_plaque / a
# Points importants de la zone de Brillouin
# Le point M
cons.POINT_M = np.array([1 / 2, 1 / (2 * np.sqrt(3))])
# Le point K
cons.POINT_K = np.array([2 / 3, 0])
# Le point Gamma
cons.GAMMA = np.array([0, 0])

with open("raw_data/constantes.log", "w") as fp:
    lines = [
            f"Dimension matrices: {n_mat}",
            f"Nombre de points dans le chemin: {sum(PATH)}",
            f"Nombre de diagonales: {NDIAG}",
            f"Pas du réseau: {a}",
            f"Nombre de points dans une cellule élémentaire: {nc}",
            f"Epsilon: {eps}",
            f"Potentiel appliqué: {v_0}",
            f"Distance r_0: {r_0}",
            f"Distance des plaques: {distance_plaque}",
            f"Vecteurs de base:",
            f"{a1=}, {a2=}",
            f"{b1=}, {b2=}",
            #f"{c1=}, {c2=}",
            f"Points importants de la zone de Brillouin:",
            f"{POINT_M=}",
            f"{POINT_K=}",
            f"{GAMMA=}"
    ]
    fp.writelines(s + "\n" for s in lines)
