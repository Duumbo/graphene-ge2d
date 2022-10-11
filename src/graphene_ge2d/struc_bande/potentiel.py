#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Potentiel.py
Programme pour faire une carte du potentiel selon le z
"""
import sys
from tqdm import tqdm
import numpy as np
from matplotlib import cbook, cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import seaborn as sns
from .matrice_coef import coeffs, coeffs_sim
from . import _constantes as cons  # v_0, distance_plaque


this = sys.modules[__name__]


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

    return np.real(out), xv, yv


def __main__():
    this.fig, this.ax = plt.subplots(subplot_kw=dict(projection="3d"))
    ls = LightSource(270, 45)
    data, x_data, y_data = gen_heat_map()
    rgb = ls.shade(data, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
    this.surf = this.ax.plot_surface(x_data, y_data, data,
                           rstride=1, cstride=1,
                           facecolors=rgb, linewidth=0,
                           antialiased=False, shade=False)
    ax.set_title("Carte du potentiel")
    np.save("raw_data/potentiel.npy", data)

    def update_pot_map():
        data, x_data, y_data = gen_heat_map()
        this.ax.clear()
        this.surf = this.ax.plot_surface(x_data, y_data, data,
                           rstride=1, cstride=1,
                           facecolors=rgb, linewidth=0,
                           antialiased=False, shade=False)
        this.fig.canvas.draw_idle()
    this.update_pot_map = update_pot_map
    fig.savefig("Images/map_potentiel.png")
