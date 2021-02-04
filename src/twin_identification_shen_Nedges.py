import numpy as np
from scipy.stats import pearsonr
import functions.twin_IDs as tid
import random
import functions.general_functions as fun

# Monozygotic twins
with open('./../data/Rest1/list_of_ID_MZ_R1R2') as fp:
    lines_mz = fp.read().split("\n")
number_of_MZ_twins = len(lines_mz[0:len(lines_mz) - 1])

# Dizygotic twins
with open('./../data/Rest1/list_of_ID_DZ_R1R2') as fp:
    lines_dz = fp.read().split("\n")
number_of_DZ_twins = len(lines_dz[0:len(lines_dz) - 1])

# All sample
number_of_individuals = number_of_DZ_twins + number_of_MZ_twins

# Subjects_ID
for i in range(0, 380):
    with open('./../data/Rest2/list_of_ID_all_MZ_DZ_selected') as fp:
        lines = fp.read().split("\n")

# Individuals ID
ID_individuals = lines[0:len(lines) - 1]

# Twin ID
twin_ID = tid.twin_pairs_ID

# Connectomes
finn_matrices_r1 = []
finn_matrices_r2 = []


# shen_parcelation_MZ_twin
for i in range(0, 246):
    # Rest1
    rel_path_1 = './../data/Rest1/All_MZ_R1/' + lines[
        i] + '/mat_conn_finn_r1.csv'
    connectome_1 = np.reshape(np.loadtxt(rel_path_1, delimiter=','),(268,268))
    connectome_1 = connectome_1[np.triu_indices(268,1)]
    finn_matrices_r1.append(connectome_1)

    # Rest2
    rel_path_2 = './../data/Rest2/All_MZ_R2/' + lines[
        i] + '/mat_conn_finn_r2.csv'
    connectome_2 = np.reshape(np.loadtxt(rel_path_2, delimiter=','),(268,268))
    connectome_2 = connectome_2[np.triu_indices(268,1)]
    finn_matrices_r2.append(connectome_2)


# shen_parcelation_DZ_twin
for i in range(246, 380):
    # Rest1
    rel_path_1 = './../data/Rest1/All_DZ_R1/' + lines[
        i] + '/mat_conn_finn_r1.csv'
    connectome_1 = np.reshape(np.loadtxt(rel_path_1, delimiter=','),(268,268))
    connectome_1 = connectome_1[np.triu_indices(268,1)]
    finn_matrices_r1.append(connectome_1)

    # Rest2
    rel_path_2 = './../data/Rest2/All_DZ_R2/' + lines[
        i] + '/mat_conn_finn_r2.csv'
    connectome_2 = np.reshape(np.loadtxt(rel_path_2, delimiter=','),(268,268))
    connectome_2 = connectome_2[np.triu_indices(268,1)]
    finn_matrices_r2.append(connectome_2)




iterations=100
edges = [10,50,100,500,1000,5000,10000]

accuracies_MZ=np.zeros((iterations, len(edges)))
accuracies_DZ=np.zeros((iterations, len(edges)))
accuracies_ID=np.zeros((iterations, len(edges)))

for k in range(iterations):
    print(k)
    for j in range(len(edges)):

        # Computation for all individuals identification
        corr_r1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
        corr_r2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
        corr_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
        corr_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
        ID_corr_r1xr2 = np.zeros(
            shape=(number_of_individuals, number_of_individuals))
        ID_corr_r2xr1 = np.zeros(
            shape=(number_of_individuals, number_of_individuals))

        # Defition of the list of connectivity matrices

        random_edges_index = random.sample(range(0, 35778), edges[j])

        matrices_r1 = np.array(finn_matrices_r1).T[random_edges_index].T
        matrices_r2 = np.array(finn_matrices_r2).T[random_edges_index].T


        # Twin identification by ROIsxROIs matrices
        # comparisons
        # Correlation matrix among the subjects
        # Rest1xRest1
        fun.corr(matrices_r1, corr_r1)

        # Rest2xRest2
        fun.corr(matrices_r2, corr_r2)

        # Rest1xRest2
        for i in range(0, len(corr_r1)):
            corr_row = []
            for m in range(0, len(corr_r1)):
                corr_row.append(
                    pearsonr(matrices_r1[i].reshape(-1),
                             matrices_r2[m].reshape(-1))[0])
            corr_r1xr2[i] = corr_row
            ID_corr_r1xr2[i] = corr_row
        np.fill_diagonal(corr_r1xr2, -1)

        # Rest2xRest1
        for i in range(0, len(corr_r1)):
            corr_row = []
            for m in range(0, len(corr_r1)):
                corr_row.append(
                    pearsonr(matrices_r2[i].reshape(-1),
                             matrices_r1[m].reshape(-1))[0])
            corr_r2xr1[i] = corr_row
            ID_corr_r2xr1[i] = corr_row
        np.fill_diagonal(corr_r2xr1, -1)

        corr_matrices = [corr_r1, corr_r2, corr_r1xr2, corr_r2xr1]

        # Individual ID
        corr_matrices_ID = {
            'n0_SI': [ID_corr_r1xr2, ID_corr_r2xr1]}
        identification_accuracy = fun.accuracy(corr_matrices_ID['n0_SI'],
                                               ID_individuals,
                                               ID_individuals)

        accuracies_ID[k, j] = np.mean(identification_accuracy)


        # Monozygotic twin identifications
        MZ_id_accuracies = []
        # for i in range(0, len(corr_matrices_whole_MZ)):
        for i in range(0, len(corr_matrices)):
            MZ_id_accuracies.append(fun.accuracy_t_MZ(
                corr_matrices[i][:246], ID_individuals,
                tid.twin_pairs_ID[:246]))

        # Dizygotic twins identifications
        DZ_id_accuracies = []
        for i in range(0, len(corr_matrices)):
            DZ_id_accuracies.append(
                fun.accuracy_t_DZ(corr_matrices[i][246:], ID_individuals,
                                  tid.twin_pairs_ID[246:]))

        accuracies_MZ[k, j] = np.mean(MZ_id_accuracies)
        accuracies_DZ[k, j] = np.mean(DZ_id_accuracies)

    print(accuracies_ID[k])
    print(accuracies_MZ[k])
    print(accuracies_DZ[k])

# np.savez('./../outputs/accuracies_SI_slices.npz',list = accuracies_ID)
# np.savez('./../outputs/accuracies_MZ_slices.npz',list = accuracies_MZ)
# np.savez('./../outputs/accuracies_DZ_slices.npz',list = accuracies_DZ)
