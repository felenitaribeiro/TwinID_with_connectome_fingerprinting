#!/usr/bin/env python3
"""
Code for individual identification based on functional connectivity
"""
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import xlrd

import twin_IDs as tid
import networks_finn as nf
import functions_fernanda as fun

# Monozygotic twins
with open('./Dados/Rest1/list_of_ID_MZ_R1R2') as fp:
    lines_mz = fp.read().split("\n")
number_of_MZ_twins = len(lines_mz[0:len(lines_mz) - 1])

# Dizygotic twins
with open('./Dados/Rest1/list_of_ID_DZ_R1R2') as fp:
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
finn_matrices_r1 = []
finn_matrices_r2 = []

# Subjects_ID
for i in range(0, 380):
    with open('./Dados/Rest2/list_of_ID_all_MZ_DZ_selected') as fp:
        lines = fp.read().split("\n")

# finn_parcelation_MZ_twin
for i in range(0, 246):
    # Rest1
    rel_path_1 = './Dados/Rest1/All_MZ_R1/' + lines[i] + '/mat_conn_finn_r1.csv'
    finn_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    # Rest2
    rel_path_2 = './Dados/Rest2/All_MZ_R2/' + lines[i] + '/mat_conn_finn_r2.csv'
    finn_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

# finn_parcelation_DZ_twin
for i in range(246, 380):
    # Rest1
    rel_path_1 = './Dados/Rest1/All_DZ_R1/' + lines[i] + '/mat_conn_finn_r1.csv'
    finn_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    # Rest2
    rel_path_2 = './Dados/Rest2/All_DZ_R2/' + lines[i] + '/mat_conn_finn_r2.csv'
    finn_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

# Individuals ID
ID_individuals = lines[0:len(lines) - 1]

########################Individuals identification by ROIsxROIs matrices comparisons########################


# Correlation matrix among the subjects
# Rest1xRest1
fun.corr(finn_matrices_r1, corr_r1)
# Rest2xRest2
fun.corr(finn_matrices_r2, corr_r2)
# Rest1xRest2
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(finn_matrices_r1[i].reshape(-1), finn_matrices_r2[m].reshape(-1))[0])
    corr_r1xr2[i] = corr_row

# Rest2xRest1
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(finn_matrices_r2[i].reshape(-1), finn_matrices_r1[m].reshape(-1))[0])
    corr_r2xr1[i] = corr_row

# TODO: Make pyplot save figure in the outputs directory
# To plot the heatmap of the matrix subjectsXsubjects just uncomment the following lines
# plt.imshow(corr_r1xr2, cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()


# Identifications
corr_matrices = [corr_r1, corr_r2, corr_r1xr2, corr_r2xr1]
fun.accuracy(corr_matrices, ID_individuals)

'''
#Non-Parametric test
accuracies_list=[]
i=0
for i in range(0,10000):
    permuted_ID=np.random.permutation(lines)
    permuted_ID=permuted_ID.tolist()
    accuracies_list.append(fun.accuracy(corr_matrices,permuted_ID,lines))
    i=i+1
    #print i

accuracies_list=np.array(accuracies_list)
print accuracies_list[accuracies_list>Identification_accuracy]
'''

# Lists of correlation values for violin plots and estimation of genetic factors contribution (to save as lists, use the np function savez, such as np.savez('name_to_be_saved', list=list_of_corr_values)
# List of twin pairs
workbook = xlrd.open_workbook('/home/hardleste/Desktop/HCP_DATA_12_2017/python_scripts_ESTUDO1/List_twin_pairs.xlsx')
sheet = workbook.sheet_by_index(0)
twin_pairs_ID = []
for rowx in range(sheet.nrows):
    cols = sheet.row_values(rowx)
    twin_pairs_ID.append(cols)
############################################################
list_of_corr_DZ_finn = []
list_of_corr_MZ_finn = []  # Correlation values of twin(rest1xrest2)
list_of_corr_nonrelated_finn = []  # Correlation values of nonrelated individuals (rest1xrest2)
list_of_corr_SI_finn = []  # Correlation values of the same individual. SI (rest1xrest2)
for i in range(0, 246):
    for j in range(0, len(corr_r1xr2)):
        if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[i]:
            list_of_corr_MZ_finn.append(corr_r1xr2[i][j])
        elif i == j:
            list_of_corr_SI_finn.append(corr_r1xr2[i][j])
        else:
            list_of_corr_nonrelated_finn.append(corr_r1xr2[i][j])

for i in range(246, 380):
    for j in range(0, len(corr_r1xr2)):
        if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[i]:
            list_of_corr_DZ_finn.append(corr_r1xr2[i][j])
        elif i == j:
            list_of_corr_SI_finn.append(corr_r1xr2[i][j])
        else:
            list_of_corr_nonrelated_finn.append(corr_r1xr2[i][j])

#######################Individuals identification by networks matrices comparisons########################
matrices_network1_r1 = []
matrices_network2_r1 = []
matrices_network3_r1 = []
matrices_network4_r1 = []
matrices_network5_r1 = []
matrices_network6_r1 = []
matrices_network7_r1 = []
matrices_network8_r1 = []

matrices_network1_r2 = []
matrices_network2_r2 = []
matrices_network3_r2 = []
matrices_network4_r2 = []
matrices_network5_r2 = []
matrices_network6_r2 = []
matrices_network7_r2 = []
matrices_network8_r2 = []

corr_network1_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network2_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network3_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network4_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network5_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network6_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network7_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network8_r1xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))

corr_network1_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network2_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network3_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network4_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network5_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network6_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network7_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network8_r2xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))

corr_network1_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network2_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network3_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network4_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network5_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network6_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network7_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network8_r1xr2 = np.zeros(shape=(number_of_individuals, number_of_individuals))

corr_network1_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network2_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network3_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network4_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network5_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network6_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network7_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))
corr_network8_r2xr1 = np.zeros(shape=(number_of_individuals, number_of_individuals))

# Juntar de matrizes rest1xrest1 com rest2xrest1
corr_matrices_n1 = []
corr_matrices_n2 = []
corr_matrices_n3 = []
corr_matrices_n4 = []
corr_matrices_n5 = []
corr_matrices_n6 = []
corr_matrices_n7 = []
corr_matrices_n8 = []

# Lists of corr for each network
list_of_corr_MZ_finn_n1 = []
list_of_corr_DZ_finn_n1 = []
list_of_corr_nonrelated_finn_n1 = []
list_of_corr_SI_finn_n1 = []

list_of_corr_MZ_finn_n2 = []
list_of_corr_DZ_finn_n2 = []
list_of_corr_nonrelated_finn_n2 = []
list_of_corr_SI_finn_n2 = []

list_of_corr_MZ_finn_n3 = []
list_of_corr_DZ_finn_n3 = []
list_of_corr_nonrelated_finn_n3 = []
list_of_corr_SI_finn_n3 = []

list_of_corr_MZ_finn_n4 = []
list_of_corr_DZ_finn_n4 = []
list_of_corr_nonrelated_finn_n4 = []
list_of_corr_SI_finn_n4 = []

list_of_corr_MZ_finn_n5 = []
list_of_corr_DZ_finn_n5 = []
list_of_corr_nonrelated_finn_n5 = []
list_of_corr_SI_finn_n5 = []

list_of_corr_MZ_finn_n6 = []
list_of_corr_DZ_finn_n6 = []
list_of_corr_nonrelated_finn_n6 = []
list_of_corr_SI_finn_n6 = []

list_of_corr_MZ_finn_n7 = []
list_of_corr_DZ_finn_n7 = []
list_of_corr_nonrelated_finn_n7 = []
list_of_corr_SI_finn_n7 = []

list_of_corr_MZ_finn_n8 = []
list_of_corr_DZ_finn_n8 = []
list_of_corr_nonrelated_finn_n8 = []
list_of_corr_SI_finn_n8 = []

number_of_networks_finn = 8
# Identification analysis for each functional network a priori defined
# plt.figure()
for j in range(1, number_of_networks_finn + 1):
    fun.network(eval('nf.network' + str(j)), finn_matrices_r1, eval('matrices_network' + str(j) + '_r1'))
    fun.network(eval('nf.network' + str(j)), finn_matrices_r2, eval('matrices_network' + str(j) + '_r2'))

    # print len(eval('matrices_network'+str(j)+'_r1'))
    # print len(eval('matrices_network' + str(j) + '_r2'))

    # Correlation matrix among the subjects
    # Rest1xRest1
    fun.corr(eval('matrices_network' + str(j) + '_r1'), eval('corr_network' + str(j) + '_r1xr1'))
    # Rest2xRest2
    fun.corr(eval('matrices_network' + str(j) + '_r2'), eval('corr_network' + str(j) + '_r2xr2'))
    # Rest1xRest2
    for i in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
        corr_row = []
        for m in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
            corr_row.append(pearsonr(eval('matrices_network' + str(j) + '_r1')[i].reshape(-1),
                                     eval('matrices_network' + str(j) + '_r2')[m].reshape(-1))[0])
        eval('corr_network' + str(j) + '_r1xr2')[i] = corr_row
    # Rest2xRest1
    for i in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
        corr_row = []
        for m in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
            corr_row.append(pearsonr(eval('matrices_network' + str(j) + '_r2')[i].reshape(-1),
                                     eval('matrices_network' + str(j) + '_r1')[m].reshape(-1))[0])
        eval('corr_network' + str(j) + '_r2xr1')[i] = corr_row

    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r1xr1'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r2xr2'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r1xr2'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r2xr1'))

    # Heatmap
    # fig = plt.subplot(2, 4, j)
    # im = fig.imshow(eval('corr_network' + str(j) + '_r1xr2'), cmap='hot', interpolation='nearest')
    # plt.title(j)

    # Identifications
    print(fun.accuracy(eval('corr_matrices_n' + str(j)), ID_individuals))

# plt.show()


# Lists of correlation values for violin plots and estimation of genetic factors contribution
for k in range(1, number_of_networks_finn + 1):
    for i in range(0, 246):
        for j in range(0, len(corr_r1xr2)):
            if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[i]:
                eval('list_of_corr_MZ_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])
            elif i == j:
                eval('list_of_corr_SI_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])
            else:
                eval('list_of_corr_nonrelated_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])

    for i in range(246, 380):
        for j in range(0, len(corr_r1xr2)):
            if [int(ID_individuals[i]), int(ID_individuals[j])] == twin_pairs_ID[i]:
                eval('list_of_corr_DZ_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])
            elif i == j:
                eval('list_of_corr_SI_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])
            else:
                eval('list_of_corr_nonrelated_finn_n' + str(k)).append(eval('corr_network' + str(k) + '_r1xr2')[i][j])

# To save these list, just uncomment the line below
'''for k in range(1,number_of_networks_finn+1):
    np.savez('list_of_corr_SI_finn_n' + str(k), list=eval('list_of_corr_SI_finn_n' + str(k)))
    np.savez('list_of_corr_MZ_finn_n' + str(k), list=eval('list_of_corr_MZ_finn_n' + str(k)))
    np.savez('list_of_corr_DZ_finn_n' + str(k), list=eval('list_of_corr_DZ_finn_n' + str(k)))
    np.savez('list_of_corr_nonrelated_finn_n' + str(k), list=eval('list_of_corr_nonrelated_finn_n' + str(k)))'''
