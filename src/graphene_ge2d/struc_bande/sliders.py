#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Module sliders
    S'occupe de gérer les sliders pour tous les graphiques.
"""
import sys
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from . import _constantes as cons
from . import _runner as run


widgets = sys.modules[__name__]


def _activate():
    fig, ax = plt.subplots()
    fig.clf()

    _ax_potentiel = plt.axes([0.25, 0.2, 0.65, 0.03])
    widgets.potentiel = Slider(_ax_potentiel, r"Potentiel $(V_0)$", -300, 300, cons.v_0)

    def update_potentiel(_):
        v_0 = widgets.potentiel.val
        cons.v_0 = v_0
        cons.update_constantes()
        run.run_all()

    _ax_r_0 = plt.axes([0.25, 0.15, 0.65, 0.03])
    widgets.r_0 = Slider(_ax_r_0, r"Dimension trou $(r_0)$", 0, 50, cons.r_0)

    def update_r_0(_):
        cons.r_0 = widgets.r_0.val
        cons.update_constantes()
        run.run_all()

    _ax_z = plt.axes([0.25, 0.25, 0.65, 0.03])
    widgets.z = Slider(_ax_z, r"Distance plaque $(z)$", 0, 50, cons.distance_plaque)

    def update_z(_):
        cons.distance_plaque = widgets.z.val
        cons.update_constantes()
        run.run_all()

    _ax_a = plt.axes([0.25, 0.30, 0.65, 0.03])
    widgets.a = Slider(_ax_a, r"Pas du réseau $(a)$", 0, 50, cons.a)

    def update_a(_):
        cons.a = widgets.a.val
        cons.update_constantes()
        run.run_all()

    widgets.update_potentiel = update_potentiel
    widgets.update_r_0 = update_r_0
    widgets.update_z = update_z
    widgets.update_a = update_a

    _ax_update_all = plt.axes([0.8, 0.025, 0.1, 0.04])
    widgets.update_button = Button(_ax_update_all, "Update")

    _ax_update_pot = plt.axes([0.6, 0.025, 0.1, 0.04])
    widgets.update_pot_button = Button(_ax_update_pot, "Update Potentiel")

    def update_pot(event):
        run.run_pot_map()

    _ax_update_struct = plt.axes([0.4, 0.025, 0.1, 0.04])
    widgets.update_struct_button = Button(_ax_update_struct, "Update Structure de Bande")

    def update_struct(event):
        run.run_struct_bande()

    #widgets.update_all = update_all
    widgets.update_struct = update_struct
    widgets.update_pot = update_pot

    #widgets.update_button.on_clicked(update_all)
    widgets.potentiel.on_changed(widgets.update_potentiel)
    widgets.r_0.on_changed(widgets.update_r_0)
    widgets.z.on_changed(widgets.update_z)
    widgets.a.on_changed(widgets.update_a)
