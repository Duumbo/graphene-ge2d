#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Définition des constantes
"""
import numpy as np

# Constantes
SIZE_MATRICES = 8
PATH = [200, 200, 200]  # Nombres de points dans le parcours de prise de données. En ordre selon la partie
NDIAG = 3  # Nombre de diagonales au centre de la matrice creuse (doit être impair)
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
z_bar = distance_plaque / a
# Points importants de la zone de Brillouin
# Le point M
POINT_M = np.array([np.pi, np.pi / (np.sqrt(3))])
# Le point K
POINT_K =  np.array([4 * np.pi / 3, 0])
# Le point Gamma
GAMMA = np.array([0, 0])

