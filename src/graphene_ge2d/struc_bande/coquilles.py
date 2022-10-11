#!/usr/bin/env python
# -*- coding-UFT-8 -*-
"""
Génère les G par coquiles
"""
import numpy as np
import matplotlib.pyplot as plt
from . import _constantes as cons


def __main__():
    n_coq = 15
    g_max = (4 / 3) * n_coq

    # Gen toutes les permutations de composantes possibles
    list_perm = []
    for i in range(- n_coq, n_coq + 1):
        for j in range(- n_coq, n_coq + 1):
            x = j
            y = (j / np.sqrt(3)) + (2 * i /np.sqrt(3))
            if x*x + y*y < g_max:
                list_perm.append(np.array([j, (j / np.sqrt(3)) + (2 * i /np.sqrt(3))]))
                #print(np.array([j, (j / np.sqrt(3)) + (2 * i /np.sqrt(3))]))

    np.save("raw_data/vecteurs.npy", np.array(list_perm))
    list_perm = np.array(list_perm)
    print(list_perm)

    fig, ax = plt.subplots()
    ax.plot(list_perm[:, 0], list_perm[:, 1], "o")
    plt.show()


if __name__ == "__main__":
    __main__()
