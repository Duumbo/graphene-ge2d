#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Test du chemin pris dans la zone de Brillouin.
"""
import numpy as np
import matplotlib.pyplot as plt
from graphene_ge2d.struc_bande.write_matrix import tracer_chemin


def test_chemin():
    path, pond = tracer_chemin()

    fig, ax = plt.subplots()
    ax.plot(path[:, 0], path[:, 1], "o", markersize=2)
    ratio = path[1:, 1] / path[1:, 0]
    section1 = np.append(np.array([True]), np.isclose(ratio, ratio[0]))
    section3 = np.append(np.array([False]), ratio == 0.)
    section2 = np.logical_not(np.logical_or(section1, section3))
    point_m = path[section1, 0][-1], path[section1, 1][-1]
    point_k = path[section2, 0][-1], path[section2, 1][-1]
    point_gamma = path[:, 0][-1], path[:, 1][-1]
    angle_gamma = np.arctan(ratio[0])

    plt.annotate("Point M", point_m)
    plt.annotate("Point K", point_k)
    plt.annotate(r"Point $\Gamma$", point_gamma, (point_gamma[0], point_gamma[1] + 0.03))
    plt.annotate(r"$\alpha=" + str(np.round(angle_gamma * 360 / (2 * np.pi), 4)) + "$", point_gamma)
    fig.savefig("Images/chemin.png")

