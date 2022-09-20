#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Test du chemin pris dans la zone de Brillouin.
"""
import matplotlib.pyplot as plt
from graphene_ge2d.struc_bande.write_matrix import tracer_chemin


def test_chemin():
    path, pond = tracer_chemin()

    fig, ax = plt.subplots()
    ax.plot(path[:, 0], path[:, 1], "o", markersize=2)
    fig.savefig("Images/chemin.png")

