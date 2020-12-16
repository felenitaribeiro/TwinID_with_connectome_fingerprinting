import numpy as np
from scipy.stats import pearsonr
import functions.twin_IDs as tid
import functions.networks_gordon as ng
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

# Computation for all individuals identification
corr_r1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_r2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))

# Defition of the list of connectivity matrices
gordon_matrices_r1 = []
gordon_matrices_r2 = []

# Subjects_ID
for i in range(0, 380):
    with open('./../data/Rest2/list_of_ID_all_MZ_DZ_selected') as fp:
        lines = fp.read().split("\n")

# gordon_parcelation_MZ_twin
for i in range(0, 246):
    # Rest1
    rel_path_1 = './../data/Rest1/All_MZ_R1/' + lines[
        i] + '/mat_conn_gordon_r1.csv'
    gordon_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    # Rest2
    rel_path_2 = './../data/Rest2/All_MZ_R2/' + lines[
        i] + '/mat_conn_gordon_r2.csv'
    gordon_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

# gordon_parcelation_DZ_twin
for i in range(246, 380):
    # Rest1
    rel_path_1 = './../data/Rest1/All_DZ_R1/' + lines[
        i] + '/mat_conn_gordon_r1.csv'
    gordon_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    # Rest2
    rel_path_2 = './../data/Rest2/All_DZ_R2/' + lines[
        i] + '/mat_conn_gordon_r2.csv'
    gordon_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

# Individuals ID
ID_individuals = lines[0:len(lines) - 1]

# Twin ID
twin_ID = tid.twin_pairs_ID

# Twin identification by ROIsxROIs matrices
# comparisons
# Correlation matrix among the subjects
# Rest1xRest1
fun.corr(gordon_matrices_r1, corr_r1)

# Rest2xRest2
fun.corr(gordon_matrices_r2, corr_r2)

# Rest1xRest2
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(gordon_matrices_r1[i].reshape(-1),
                     gordon_matrices_r2[m].reshape(-1))[0])
    corr_r1xr2[i] = corr_row
np.fill_diagonal(corr_r1xr2, -1)

# Rest2xRest1
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(gordon_matrices_r2[i].reshape(-1),
                     gordon_matrices_r1[m].reshape(-1))[0])
    corr_r2xr1[i] = corr_row
np.fill_diagonal(corr_r2xr1, -1)

corr_matrices = [corr_r1, corr_r2, corr_r1xr2, corr_r2xr1]

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

# Permutation - randon twin_pair_ID
accuracy_random_MZ = []
accuracy_random_DZ = []
for j in range(0, 1000):
    accuracies_MZ = []
    accuracies_DZ = []
    for i in range(0, len(corr_matrices)):
        accuracies_MZ.append(fun.accuracy_t_MZ(
            corr_matrices[i][:246], ID_individuals,
            np.ndarray.tolist(
                np.random.permutation(tid.twin_pairs_ID[:246]))))
    accuracy_random_MZ.append(np.mean(accuracies_MZ))

    for i in range(0, len(corr_matrices)):
        accuracies_DZ.append(fun.accuracy_t_DZ(
            corr_matrices[i][246:], ID_individuals,
            np.ndarray.tolist(
                np.random.permutation(tid.twin_pairs_ID[246:]))))
    accuracy_random_DZ.append(np.mean(accuracies_DZ))

# Twins identification for individual nets
matrices_network_r1 = {'network1': [], 'network2': [], 'network3': [],
                       'network4': [], 'network5': [], 'network6': [],
                       'network7': [], 'network8': [], 'network9': [],
                       'network10': [], 'network11': [], 'network12': []
                       }
matrices_network_r2 = {'network1': [], 'network2': [], 'network3': [],
                       'network4': [], 'network5': [], 'network6': [],
                       'network7': [], 'network8': [], 'network9': [],
                       'network10': [], 'network11': [], 'network12': []
                       }

corr_network = {}

# Combine rest1xrest1 and rest2xrest1
corr_matrices_n1 = []
corr_matrices_n2 = []
corr_matrices_n3 = []
corr_matrices_n4 = []
corr_matrices_n5 = []
corr_matrices_n6 = []
corr_matrices_n7 = []
corr_matrices_n8 = []
corr_matrices_n9 = []
corr_matrices_n10 = []
corr_matrices_n11 = []
corr_matrices_n12 = []

# accuracy
accuracies = {'n1_MZ': [], 'n1_DZ': [], 'n2_MZ': [], 'n2_DZ': [],
              'n3_MZ': [], 'n3_DZ': [], 'n4_MZ': [], 'n4_DZ': [],
              'n5_MZ': [], 'n5_DZ': [], 'n6_MZ': [], 'n6_DZ': [],
              'n7_MZ': [], 'n7_DZ': [], 'n8_MZ': [], 'n8_DZ': [],
              'n9_MZ': [], 'n9_DZ': [], 'n10_MZ': [], 'n10_DZ': [],
              'n11_MZ': [], 'n11_DZ': [], 'n12_MZ': [], 'n12_DZ': []
              }

accuracies_random = {'n1_MZ': [], 'n1_DZ': [], 'n2_MZ': [], 'n2_DZ': [],
                     'n3_MZ': [], 'n3_DZ': [], 'n4_MZ': [], 'n4_DZ': [],
                     'n5_MZ': [], 'n5_DZ': [], 'n6_MZ': [], 'n6_DZ': [],
                     'n7_MZ': [], 'n7_DZ': [], 'n8_MZ': [], 'n8_DZ': [],
                     'n9_MZ': [], 'n9_DZ': [], 'n10_MZ': [], 'n10_DZ': [],
                     'n11_MZ': [], 'n11_DZ': [], 'n12_MZ': [], 'n12_DZ': []
                     }

number_of_networks_gordon = 12
for j in range(1, number_of_networks_gordon + 1):
    fun.network(eval('ng.network' + str(j)), gordon_matrices_r1,
                matrices_network_r1['network' + str(j)])
    fun.network(eval('ng.network' + str(j)), gordon_matrices_r2,
                matrices_network_r2['network' + str(j)])

    # Correlation matrix among the subjects
    corr_network['n' + str(j) + '_r1xr1'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))
    corr_network['n' + str(j) + '_r2xr2'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))
    corr_network['n' + str(j) + '_r1xr2'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))
    corr_network['n' + str(j) + '_r2xr1'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))
    # Rest1xRest1
    fun.corr(matrices_network_r1['network' + str(j)],
             corr_network['n' + str(j) + '_r1xr1'])
    # Rest2xRest2
    fun.corr(matrices_network_r2['network' + str(j)],
             corr_network['n' + str(j) + '_r2xr2'])

    # Rest1xRest2
    for i in range(0, len(corr_network['n' + str(j) + '_r1xr1'])):
        corr_row = []
        for m in range(0, len(corr_network['n' + str(j) + '_r1xr1'])):
            corr_row.append(pearsonr(
                matrices_network_r1['network' + str(j)][i].reshape(-1),
                matrices_network_r2['network' + str(j)][m].reshape(-1))[0])
        corr_network['n' + str(
            j) + '_r1xr2'][i] = corr_row
    np.fill_diagonal(corr_network['n' + str(j) + '_r1xr2'], -1)

    # Rest2xRest1
    for i in range(0, len(corr_network['n' + str(j) + '_r1xr1'])):
        corr_row = []
        for m in range(0, len(corr_network['n' + str(j) + '_r1xr1'])):
            corr_row.append(pearsonr(
                matrices_network_r2['network' + str(j)][i].reshape(-1),
                matrices_network_r1['network' + str(j)][m].reshape(-1))[0])
        corr_network['n' + str(j) + '_r2xr1'][i] = corr_row
    np.fill_diagonal(corr_network['n' + str(j) + '_r2xr1'], -1)

    eval('corr_matrices_n' + str(j)).append(
        corr_network['n' + str(j) + '_r1xr1'])
    eval('corr_matrices_n' + str(j)).append(
        corr_network['n' + str(j) + '_r2xr2'])
    eval('corr_matrices_n' + str(j)).append(
        corr_network['n' + str(j) + '_r1xr2'])
    eval('corr_matrices_n' + str(j)).append(
        corr_network['n' + str(j) + '_r2xr1'])

    # Identification analyses
    # Monozygotic
    print('MZ identification accuracies for network' + str(j))
    for i in range(0, len(eval('corr_matrices_n' + str(j)))):
        accuracies['n' + str(j) + '_MZ'].append(
            fun.accuracy_t_MZ(eval('corr_matrices_n' + str(j))[i][0:246, :],
                              ID_individuals, tid.twin_pairs_ID))
    print(accuracies['n' + str(j) + '_MZ'])
    # print np.mean(accuracies['n' + str(j) + '_MZ'])
    # print np.std(accuracies['n' + str(j) + '_MZ'])

    # Dizygotic
    print('DZ identification accuracies for network' + str(j))
    for i in range(0, len(eval('corr_matrices_n' + str(j)))):
        accuracies['n' + str(j) + '_DZ'].append(
            fun.accuracy_t_DZ(eval('corr_matrices_n' + str(j))[i][246:, :],
                              ID_individuals, tid.twin_pairs_ID[246:]))
    print(accuracies['n' + str(j) + '_DZ'])
    # print np.mean(accuracies['n' + str(j) + '_DZ'])
    # print np.std(accuracies['n' + str(j) + '_DZ'])

    # Permutation - random twin_pair_ID
    for k in range(0, 1000):
        accuracies_MZ = []
        accuracies_DZ = []
        for i in range(0, len(eval('corr_matrices_n' + str(j)))):
            accuracies_MZ.append(
                fun.accuracy_t_MZ(
                    eval('corr_matrices_n' + str(j))[i][0:246, :],
                    ID_individuals, np.ndarray.tolist(
                        np.random.permutation(tid.twin_pairs_ID[:246]))))
        accuracies_random['n' + str(j) + '_MZ'].append(np.mean(accuracies_MZ))

        for i in range(0, len(eval('corr_matrices_n' + str(j)))):
            accuracies_DZ.append(
                fun.accuracy_t_DZ(eval('corr_matrices_n' + str(j))[i][246:, :],
                                  ID_individuals, np.ndarray.tolist(
                        np.random.permutation(tid.twin_pairs_ID[246:]))))
        accuracies_random['n' + str(j) + '_DZ'].append(np.mean(accuracies_DZ))

accuracies['n0_MZ'] = MZ_id_accuracies
accuracies['n0_DZ'] = DZ_id_accuracies
accuracies_random['n0_MZ'] = accuracy_random_MZ
accuracies_random['n0_DZ'] = accuracy_random_DZ

for j in range(13):
    print('Network' + str(j) + '=' + str(np.sum(
        accuracies_random['n' + str(j) + '_MZ'] > np.mean(
            accuracies['n' + str(j) + '_MZ']))))
    print('Network' + str(j) + '=' + str(np.sum(
        accuracies_random['n' + str(j) + '_DZ'] > np.mean(
            accuracies['n' + str(j) + '_DZ']))))

# Saving final results
np.savez('./../outputs/accuracies_twin_id_gordon', dict=accuracies)
np.savez('./../outputs/accuracies_perm_twin_id_gordon', dict=accuracies_random)
