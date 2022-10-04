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
    nx, ny = 200, 200
    xinf, xsup = -1, 1
    yinf, ysup = -1, 1

    positions_x = np.linspace(xinf, xsup, nx)
    positions_y = np.linspace(yinf, ysup, ny)
    xv, yv = np.meshgrid(positions_x, positions_y)
    out = np.zeros(xv.shape, dtype=np.complex128)

    vecteurs = np.load("raw_data/vecteurs.npy")

    coefficients = coeffs(vecteurs) * coeffs_sim(vecteurs)

    for n, g in tqdm(enumerate(vecteurs)):
        out += coefficients * np.exp(1j * 2 * np.pi * (xv * g[0] + yv * g[1]) / cons.eta)

    return np.real(out)


def __main__():
    fig, ax = plt.subplots()
    data = gen_heat_map()
    np.save("raw_data/potentiel.npy", data)
    sns.heatmap(data)
    ax.set_title(rf"Potentiel à $V_0={cons.v_0}$ et $z={cons.distance_plaque}$")
    ax.set_aspect(1)
    ax.set_xticks([x * 20 for x in range(11)])
    ax.set_xticklabels([f"{np.round(x, 2)}" for x in np.linspace(-1, 1, 11)])
    ax.set_yticks([y * 20 for y in range(11)])
    ax.set_yticklabels([f"{np.round(y, 2)}" for y in np.linspace(-1, 1, 11)])

    plt.show()
    fig.savefig("Images/map_potentiel.png")
