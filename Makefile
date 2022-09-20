.PHONY: clean

Images/structure_de_bande.png: raw_data/chemin_norm.npy raw_data/list_eigvals.npy src/graphene_ge2d/struc_bande/generate_graph.py
	python src/graphene_ge2d/ -g

raw_data/diagonals.npy: src/graphene_ge2d/struc_bande/_constantes.py src/graphene_ge2d/struc_bande/write_matrix.py src/graphene_ge2d/struc_bande/matrice_coef.py
	python src/graphene_ge2d/ -m

raw_data/list_eigvals.npy: raw_data/diagonals.npy raw_data/chemin.npy src/graphene_ge2d/struc_bande/compute_eigenvalues.py
	python src/graphene_ge2d/ -e

clean:
	rm raw_data/chemin.npy raw_data/chemin_norm.npy raw_data/diagonals.npy raw_data/list_eigvals.npy
