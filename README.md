# Graphène Artificiel sur GE2D
Simulation de la structure du graphène en utilisant un gaz d'électrons 2D en forme de réseau honeycomb.

## Calcul de la structure de bande
Le calcul de la structure de bande est fait en diagonalisant la matrice représentant l'équation de Schrödinger
pour obtenir l'énergie selon certaine valeurs de **k**. Pour ce faire, on crée une grille dans l'espace réciproque
des différents vecteurs **G**=n**b~1~**+m**b~2~** et on écrit la matrice pour chaque point **k** du chemin à évaluer.
### Spécifications
Les matrices à diagonaliser sont de la même taille que le nombre de vecteurs **G**. Les termes de cette matrices diminuent
à zéro assez rapidement pour justifier l'utilisation de matrices multi-diagonales; le choix du nombre de diagonales utilisées
devrait être testé. Le chemin de **k** utilisé se limite à des droites entre les points importants de la zone de Brillouin.

## Structure du projet
Dans le fichier `src` se trouve le code source, écrit en python. Le script de base se trouve dans `\_\_main\_\_.py`, dans le
but de pouvoir appeler le script dans le futur avec des arguments. Le dossier `tests` regroupe les tests unitaires, gérés par
_pytest_. Le dossier `Images` est composé des différents graphiques et images générés par le code. Le dossier `raw_data` s'occupe
de garder les sauvegardes de données computationnellement dispendieuses.

## Installation
Ce projet utilise _poetry_ pour gérer les dépendances et installation des _packages_ dans les environnements virtuels. Ceci
est important, comme ce projet dépend de `numba-scipy`, qui dépend à ce jour d'une version spécifique de `scipy` et `numba`
comme ce projet est encore en développement. Pour installer les dépendances, on recommande de créer un environnement virtuel
au préalable afin de ne pas mélanger ces packages à votre installation globale. J'utilise `virtualenv` à la racine du projet:
```
python -m virtualenv .venv
```
et ensuite on peut utiliser _poetry_ pour installer les dépendances:
```
poetry install
```
_Poetry_ va automatiquement trouver l'environnement virtuel `.venv` créé précédemment. Pour utiliser ce environnement pour rouler
les programmes de ce projet, on peut activer l'environnement avec:
```
poetry shell
python src/__main__.py
```
