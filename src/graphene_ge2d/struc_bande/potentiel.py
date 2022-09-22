#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Potentiel.py
Programme pour faire une carte du potentiel selon le z
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from .matrice_coef import coeffs, coeffs_sim
from . import _constantes as cons  # v_0, distance_plaque


def gen_heat_map():
    nx, ny = 2000, 2000
    xinf, xsup = -10, 10
    yinf, ysup = -10, 10

    positions_x = np.linspace(xinf, xsup, nx)
    positions_y = np.linspace(yinf, ysup, ny)
    xv, yv = np.meshgrid(positions_x, positions_y)
    out = np.zeros(xv.shape)

    vecteurs = np.load("raw_data/vecteurs.npy")

    coefficients = coeffs(vecteurs) * coeffs_sim(vecteurs)

    for n, g in tqdm(enumerate(vecteurs)):
        out += np.real(coefficients[n] * np.exp(1j * (xv * g[0] + yv * g[1])))

    return out


def __main__():
    fig, ax = plt.subplots()
    data = gen_heat_map()
    np.save("raw_data/potentiel.npy", data)
    sns.heatmap(data)
    ax.set_title(rf"Potentiel à $V_0={cons.v_0}$ et $z={cons.distance_plaque}$")
    ax.set_aspect(1)

    plt.show()
    fig.savefig("Images/map_potentiel.png")
