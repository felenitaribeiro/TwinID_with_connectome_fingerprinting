import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import networks_gordon as ns
import xlrd
import functions as fun
import pathlib

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
corr_network = {
    'n0_r1xr2': np.zeros(shape=(number_of_individuals, number_of_individuals)),
    'n0_r2xr1': np.zeros(shape=(number_of_individuals, number_of_individuals))}

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

# Fischer z-transformed correlation matrix plot
# ax = sns.heatmap(gordon_matrices_r1[0],
#                  cbar_kws={'label': 'Fischer z-transformed correlation'})
#
# plt.yticks(np.arange(1, 268, step=50),np.arange(1, 268, step=50))
# plt.xticks(np.arange(1, 268, step=50),np.arange(1, 268, step=50),rotation=45)
# plt.show()

# Individuals ID
ID_individuals = lines[0:len(lines) - 1]

# Individuals identification by ROIsxROIs matrices
# comparisons

# Rest1xRest2
for i in range(0, len(corr_network['n0_r1xr2'])):
    corr_row = []
    for m in range(0, len(corr_network['n0_r1xr2'])):
        corr_row.append(
            pearsonr(gordon_matrices_r1[i].reshape(-1),
                     gordon_matrices_r2[m].reshape(-1))[0])
    corr_network['n0_r1xr2'][i] = corr_row

# Rest2xRest1
for i in range(0, len(corr_network['n0_r1xr2'])):
    corr_row = []
    for m in range(0, len(corr_network['n0_r1xr2'])):
        corr_row.append(
            pearsonr(gordon_matrices_r2[i].reshape(-1),
                     gordon_matrices_r1[m].reshape(-1))[0])
    corr_network['n0_r2xr1'][i] = corr_row

# # To plot the heatmap of the matrix subjectsXsubjects just uncomment the
# # following lines
# ax=sns.heatmap(corr_r1xr2,cmap='Blues',cbar_kws={'label':
# 'Pearson\'s correlation'},vmax=1)
# plt.yticks(np.arange(1, 380, step=50),np.arange(1, 380, step=50))
# plt.xticks(np.arange(1, 380, step=50),np.arange(1, 380, step=50),rotation=45)
# plt.show()


# Identifications
corr_matrices = {'n0_SI': [corr_network['n0_r1xr2'], corr_network['n0_r2xr1']]}
identification_accuracy = {
    'n0_SI': fun.accuracy(corr_matrices['n0_SI'], ID_individuals,
                          ID_individuals)}

# Non-Parametric test
accuracies_list = []
for i in range(0, 10000):
    permuted_ID = np.random.permutation(ID_individuals)
    permuted_ID = permuted_ID.tolist()
    accuracies_list.append(
        np.mean(
            fun.accuracy(corr_matrices['n0_SI'], permuted_ID, ID_individuals)))

accuracies_list = np.array(accuracies_list)
accuracies_random = {'n0_SI': accuracies_list}

# Lists of correlation values for violin plots (to save as lists, use the np
# function savez, such as
# np.savez('name_to_be_saved', list=list_of_corr_values)
# List of twin pairs
workbook = xlrd.open_workbook(
    './../data/List_twin_pairs.xlsx')
sheet = workbook.sheet_by_index(0)
twin_pairs_ID = []
for rowx in range(sheet.nrows):
    cols = sheet.row_values(rowx)
    twin_pairs_ID.append(cols)

# Correlation values of rest1xrest2
corr_distributions = {'list_of_corr_DZ_gordon': [],
                      'list_of_corr_MZ_gordon': [],
                      'list_of_corr_SI_gordon': [],
                      'list_of_corr_nonrelated_gordon': []}
for i in range(0, 246):
    for j in range(0, len(corr_network['n0_r1xr2'])):
        if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[
            i]:
            corr_distributions['list_of_corr_MZ_gordon'].append(
                corr_network['n0_r1xr2'][i][j])
        elif i == j:
            corr_distributions['list_of_corr_SI_gordon'].append(
                corr_network['n0_r1xr2'][i][j])
        else:
            corr_distributions['list_of_corr_nonrelated_gordon'].append(
                corr_network['n0_r1xr2'][i][j])

for i in range(246, 380):
    for j in range(0, len(corr_network['n0_r1xr2'])):
        if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[
            i]:
            corr_distributions['list_of_corr_DZ_gordon'].append(
                corr_network['n0_r1xr2'][i][j])
        elif i == j:
            corr_distributions['list_of_corr_SI_gordon'].append(
                corr_network['n0_r1xr2'][i][j])
        else:
            corr_distributions['list_of_corr_nonrelated_gordon'].append(
                corr_network['n0_r1xr2'][i][j])

np.savez('./../outputs/corr_distributions/list_of_corr_SI_gordon_n0',
         list=corr_distributions['list_of_corr_SI_gordon'])
np.savez('./../outputs/corr_distributions/list_of_corr_MZ_gordon_n0',
         list=corr_distributions['list_of_corr_MZ_gordon'])
np.savez('./../outputs/corr_distributions/list_of_corr_DZ_gordon_n0',
         list=corr_distributions['list_of_corr_DZ_gordon'])
np.savez('./../outputs/corr_distributions/list_of_corr_nonrelated_gordon_n0',
         list=corr_distributions['list_of_corr_nonrelated_gordon'])

# Individuals identification for individual nets
matrices_network_r1 = {'network1': [], 'network2': [], 'network3': [],
                       'network4': [], 'network5': [], 'network6': [],
                       'network7': [], 'network8': [], 'network9': [],
                       'network10': [], 'network11': [], 'network12': [],
                       'network13': []
                       }
matrices_network_r2 = {'network1': [], 'network2': [], 'network3': [],
                       'network4': [], 'network5': [], 'network6': [],
                       'network7': [], 'network8': [], 'network9': [],
                       'network10': [], 'network11': [], 'network12': [],
                       'network13': []
                       }

number_of_networks_gordon = 12
# Identification analysis for each functional network


pathlib.Path('./../outputs/corr_distributions').mkdir(parents=True, exist_ok=True)
for j in range(1, number_of_networks_gordon + 1):
    fun.network(eval('ns.network' + str(j)), gordon_matrices_r1,
                matrices_network_r1['network' + str(j)])
    fun.network(eval('ns.network' + str(j)), gordon_matrices_r2,
                matrices_network_r2['network' + str(j)])

    # print len(eval('matrices_network'+str(j)+'_r1'))
    # print len(eval('matrices_network' + str(j) + '_r2'))

    # Correlation matrix among the subjects
    corr_network['n' + str(j) + '_r1xr2'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))
    corr_network['n' + str(j) + '_r2xr1'] = np.zeros(
        shape=(number_of_individuals,
               number_of_individuals))

    # Rest1xRest2
    for i in range(0, len(corr_network['n' + str(
            j) + '_r1xr2'])):
        corr_row = []
        for m in range(0, len(corr_network['n' + str(
                j) + '_r1xr2'])):
            corr_row.append(pearsonr(
                matrices_network_r1['network' + str(j)][i].reshape(-1),
                matrices_network_r2['network' + str(j)][m].reshape(-1))[0])
        corr_network['n' + str(
            j) + '_r1xr2'][i] = corr_row

    # Rest2xRest1
    for i in range(0, len(corr_network['n' + str(
            j) + '_r1xr2'])):
        corr_row = []
        for m in range(0, len(corr_network['n' + str(
                j) + '_r1xr2'])):
            corr_row.append(pearsonr(
                matrices_network_r2['network' + str(j)][i].reshape(-1),
                matrices_network_r1['network' + str(j)][m].reshape(-1))[0])
        corr_network['n' + str(j) + '_r2xr1'][i] = corr_row

    corr_matrices['n' + str(j) + '_SI'] = [
        corr_network['n' + str(j) + '_r1xr2'],
        corr_network['n' + str(j) + '_r2xr1']
    ]

    # Identifications
    identification_accuracy['n' + str(j) + '_SI'] = fun.accuracy(
        corr_matrices['n' + str(j) + '_SI'], ID_individuals,
        ID_individuals)

    # Permutation - non-parametric test
    accuracies_list = []
    for i in range(0, 1000):
        permuted_ID = np.random.permutation(ID_individuals)
        permuted_ID = permuted_ID.tolist()
        accuracies_list.append(
            np.mean(
                fun.accuracy(corr_matrices['n' + str(j) + '_SI'], permuted_ID,
                             ID_individuals)))

    accuracies_list = np.array(accuracies_list)
    accuracies_random['n' + str(j) + '_SI'] = accuracies_list

    # Lists of correlation values for violin plots
    corr_distributions['list_of_corr_MZ_gordon_n' + str(j)] = []
    corr_distributions['list_of_corr_DZ_gordon_n' + str(j)] = []
    corr_distributions['list_of_corr_nonrelated_gordon_n' + str(j)] = []
    corr_distributions['list_of_corr_SI_gordon_n' + str(j)] = []

    for i in range(0, 246):
        for l in range(0, len(corr_matrices['n' + str(j) + '_SI'][0])):
            if [int(ID_individuals[i]), int(ID_individuals[l])] == \
                    twin_pairs_ID[i]:
                corr_distributions[
                    'list_of_corr_MZ_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])
            elif i == l:
                corr_distributions[
                    'list_of_corr_SI_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])
            else:
                corr_distributions[
                    'list_of_corr_nonrelated_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])

    for i in range(246, 380):
        for l in range(0, len(corr_matrices['n' + str(j) + '_SI'][0])):
            if [int(ID_individuals[i]), int(ID_individuals[l])] == \
                    twin_pairs_ID[i]:
                corr_distributions[
                    'list_of_corr_DZ_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])
            elif i == l:
                corr_distributions[
                    'list_of_corr_SI_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])
            else:
                corr_distributions[
                    'list_of_corr_nonrelated_gordon_n' + str(j)].append(
                    corr_matrices['n' + str(j) + '_SI'][0][i][l])

# To save these list, just uncomment the line below
for k in range(1, number_of_networks_gordon + 1):
    np.savez('./../outputs/corr_distributions/list_of_corr_SI_gordon_n' + str(k),
             list=corr_distributions['list_of_corr_SI_gordon_n' + str(k)])
    np.savez('./../outputs/corr_distributions/list_of_corr_MZ_gordon_n' + str(k),
             list=corr_distributions['list_of_corr_MZ_gordon_n' + str(k)])
    np.savez('./../outputs/corr_distributions/list_of_corr_DZ_gordon_n' + str(k),
             list=corr_distributions['list_of_corr_DZ_gordon_n' + str(k)])
    np.savez(
        './../outputs/corr_distributions/list_of_corr_nonrelated_gordon_n' + str(k),
        list=corr_distributions['list_of_corr_nonrelated_gordon_n' + str(k)])

np.savez('./../outputs/accuracies_ind_id_gordon.npz',
         dict=identification_accuracy)
np.savez('./../outputs/accuracies_perm_ind_id_gordon.npz', dict=accuracies_random)

print(identification_accuracy)

for j in range((len(identification_accuracy))):
    print('Network' + str(j) + '=' + str(np.sum(
        accuracies_random['n' + str(j) + '_SI'] > np.mean(
            identification_accuracy['n' + str(j) + '_SI']))))
