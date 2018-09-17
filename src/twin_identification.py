import numpy as np
import h5py
from scipy.stats import pearsonr
import os
import matplotlib.pyplot as plt

import twin_IDs as tid
import networks_finn as nf
import functions_fernanda as fun

#Monozygotic twins
with open('./Dados/Rest1/list_of_ID_MZ_R1R2') as fp:
    lines_mz = fp.read().split("\n")
number_of_MZ_twins=len(lines_mz[0:len(lines_mz)-1])

#Dizygotic twins
with open('./Dados/Rest1/list_of_ID_DZ_R1R2') as fp:
    lines_dz = fp.read().split("\n")
number_of_DZ_twins=len(lines_dz[0:len(lines_dz)-1])

#All sample
number_of_individuals=number_of_DZ_twins+number_of_MZ_twins

#Computation for all individuals identification
corr_r1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_r2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))

#Defition of the list of connectivity matrices
finn_matrices_r1=[]
finn_matrices_r2=[]

#Subjects_ID
for i in range(0,380):
    with open('./Dados/Rest2/list_of_ID_all_MZ_DZ_selected') as fp:
        lines = fp.read().split("\n")

#finn_parcelation_MZ_twin
for i in range(0,246):
    #Rest1
    rel_path_1 = './Dados/Rest1/All_MZ_R1/'+lines[i]+'/mat_conn_finn_r1.csv'
    finn_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    #Rest2
    rel_path_2 = './Dados/Rest2/All_MZ_R2/'+lines[i]+'/mat_conn_finn_r2.csv'
    finn_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

#finn_parcelation_DZ_twin
for i in range(246,380):
    #Rest1
    rel_path_1 = './Dados/Rest1/All_DZ_R1/'+lines[i]+'/mat_conn_finn_r1.csv'
    finn_matrices_r1.append(np.loadtxt(rel_path_1, delimiter=','))
    #Rest2
    rel_path_2 = './Dados/Rest2/All_DZ_R2/'+lines[i]+'/mat_conn_finn_r2.csv'
    finn_matrices_r2.append(np.loadtxt(rel_path_2, delimiter=','))

#Individuals ID
ID_individuals=lines[0:len(lines)-1]

#Twin ID
twin_ID=tid.twin_pairs_ID








########################Twin identification by ROIsxROIs matrices comparisons########################
#Correlation matrix among the subjects
#Rest1xRest1
fun.corr(finn_matrices_r1,corr_r1)
        #Monozygotic
corr_r1_MZ=corr_r1[0:246,:]
        #Dizygotic
corr_r1_DZ=corr_r1[246:380,:]

#Rest2xRest2
fun.corr(finn_matrices_r2,corr_r2)
        #Monozygotic
corr_r2_MZ=corr_r2[0:246,:]
        #Dizygotic
corr_r2_DZ=corr_r2[246:380,:]

#Rest1xRest2
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(finn_matrices_r1[i].reshape(-1), finn_matrices_r2[m].reshape(-1))[0])
    corr_r1xr2[i] = corr_row
np.fill_diagonal(corr_r1xr2, -1)
        #Monozygotic
corr_r1xr2_MZ=corr_r1xr2[0:246,:]
        #Dizygotic
corr_r1xr2_DZ=corr_r1xr2[246:380,:]


#Rest2xRest1
for i in range(0, len(corr_r1)):
    corr_row = []
    for m in range(0, len(corr_r1)):
        corr_row.append(
            pearsonr(finn_matrices_r2[i].reshape(-1), finn_matrices_r1[m].reshape(-1))[0])
    corr_r2xr1[i] = corr_row
np.fill_diagonal(corr_r2xr1, -1)
        #Monozygotic
corr_r2xr1_MZ=corr_r2xr1[0:246,:]
        #Dizygotic
corr_r2xr1_DZ=corr_r2xr1[246:380,:]

#Heatmap
#plt.imshow(corr_r1xr2, cmap='hot', interpolation='nearest')
#plt.show()

corr_matrices_whole_MZ=[corr_r1_MZ,corr_r2_MZ,corr_r1xr2_MZ,corr_r2xr1_MZ]
corr_matrices_whole_DZ=[corr_r1_DZ,corr_r2_DZ,corr_r1xr2_DZ,corr_r2xr1_DZ]

#Monozygotic twin identifications
def accuracy_t_MZ(corr_matrices,list_subject_ID,list_of_twin):
    acc = 0
    identified_pair=[]
    for i in range(0, len(corr_matrices)):
        identified_pair.append([int(list_subject_ID[i]),int(list_subject_ID[np.argmax(corr_matrices[i])])])
        if identified_pair[i]==list_of_twin[i]:
           acc = acc + 1.
    accuracy_t = (acc / len(corr_matrices)) * 100.
    return accuracy_t

MZ_id_accuracies=[]
for i in range(0, len(corr_matrices_whole_MZ)):
    MZ_id_accuracies.append(accuracy_t_MZ(corr_matrices_whole_MZ[i],ID_individuals,tid.twin_pairs_ID))

#Dizygotic
def accuracy_t_DZ(corr_matrices,list_subject_ID,list_of_twin):
    acc = 0
    identified_pair=[]
    for i in range(0, len(corr_matrices)):
        identified_pair.append([int(list_subject_ID[i+246]),int(list_subject_ID[np.argmax(corr_matrices[i])])])
        if identified_pair[i]==list_of_twin[i]:
           acc = acc + 1.
    #print len(identified_pair)
    accuracy_t = (acc / len(corr_matrices)) * 100.
    return accuracy_t

DZ_id_accuracies=[]
for i in range(0, len(corr_matrices_whole_DZ)):
    DZ_id_accuracies.append(accuracy_t_DZ(corr_matrices_whole_DZ[i],ID_individuals,tid.twin_pairs_ID[246:]))



####Permutation - randon twin_pair_ID
accuracy_random_MZ=[]
accuracy_random_DZ=[]
for j in range(0,1000):
    accuracies_MZ = []
    accuracies_DZ = []
    for i in range(0, len(corr_matrices_whole_MZ)):
        accuracies_MZ.append(accuracy_t_MZ(corr_matrices_whole_MZ[i], ID_individuals, np.ndarray.tolist(
np.random.permutation(tid.twin_pairs_ID))))
    accuracy_random_MZ.append(np.mean(accuracies_MZ))

    for i in range(0, len(corr_matrices_whole_DZ)):
        accuracies_DZ.append(accuracy_t_DZ(corr_matrices_whole_DZ[i], ID_individuals, np.ndarray.tolist(
np.random.permutation(tid.twin_pairs_ID[246:]))))
    accuracy_random_DZ.append(np.mean(accuracies_DZ))



#######################Individuals identification by networks matrices comparisons########################
matrices_network1_r1=[]
matrices_network2_r1=[]
matrices_network3_r1=[]
matrices_network4_r1=[]
matrices_network5_r1=[]
matrices_network6_r1=[]
matrices_network7_r1=[]
matrices_network8_r1=[]

matrices_network1_r2=[]
matrices_network2_r2=[]
matrices_network3_r2=[]
matrices_network4_r2=[]
matrices_network5_r2=[]
matrices_network6_r2=[]
matrices_network7_r2=[]
matrices_network8_r2=[]

corr_network1_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network2_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network3_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network4_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network5_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network6_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network7_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network8_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))

corr_network1_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network2_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network3_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network4_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network5_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network6_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network7_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network8_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))

corr_network1_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network2_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network3_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network4_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network5_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network6_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network7_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network8_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))

corr_network1_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network2_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network3_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network4_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network5_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network6_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network7_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
corr_network8_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))


#Juntar de matrizes rest1xrest1 com rest2xrest1
corr_matrices_n1=[]
corr_matrices_n2=[]
corr_matrices_n3=[]
corr_matrices_n4=[]
corr_matrices_n5=[]
corr_matrices_n6=[]
corr_matrices_n7=[]
corr_matrices_n8=[]

# accuracy
accuracies_n1_MZ = []
accuracies_n2_MZ = []
accuracies_n3_MZ = []
accuracies_n4_MZ = []
accuracies_n5_MZ = []
accuracies_n6_MZ = []
accuracies_n7_MZ = []
accuracies_n8_MZ = []

accuracies_n1_DZ = []
accuracies_n2_DZ = []
accuracies_n3_DZ = []
accuracies_n4_DZ = []
accuracies_n5_DZ = []
accuracies_n6_DZ = []
accuracies_n7_DZ = []
accuracies_n8_DZ = []

number_of_networks_finn=8

#plt.figure()
for j in range(1,number_of_networks_finn+1):
    fun.network(eval('nf.network'+str(j)),finn_matrices_r1,eval('matrices_network'+str(j)+'_r1'))
    fun.network(eval('nf.network' + str(j)), finn_matrices_r2, eval('matrices_network' + str(j) +'_r2'))

    #print len(eval('matrices_network'+str(j)+'_r1'))
    #print len(eval('matrices_network' + str(j) + '_r2'))

    # Correlation matrix among the subjects
    # Rest1xRest1
    fun.corr(eval('matrices_network' + str(j) + '_r1'), eval('corr_network' + str(j)+ '_r1xr1'))
    # Rest2xRest2
    fun.corr(eval('matrices_network' + str(j) + '_r2'), eval('corr_network' + str(j)+ '_r2xr2'))
    # Rest1xRest2
    for i in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
        corr_row = []
        for m in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
            corr_row.append(pearsonr(eval('matrices_network' + str(j) + '_r1')[i].reshape(-1),
                                     eval('matrices_network' + str(j) + '_r2')[m].reshape(-1))[0])
        eval('corr_network' + str(j) + '_r1xr2')[i] = corr_row
    np.fill_diagonal(eval('corr_network' + str(j) + '_r1xr2'), -1)
    #Rest2xRest1
    for i in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
        corr_row = []
        for m in range(0, len(eval('corr_network' + str(j) + '_r1xr1'))):
            corr_row.append(pearsonr(eval('matrices_network' + str(j) + '_r2')[i].reshape(-1),
                                     eval('matrices_network' + str(j) + '_r1')[m].reshape(-1))[0])
        eval('corr_network' + str(j) + '_r2xr1')[i] = corr_row
    np.fill_diagonal(eval('corr_network' + str(j) + '_r2xr1'), -1)


    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j)+ '_r1xr1'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j)+ '_r2xr2'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r1xr2'))
    eval('corr_matrices_n' + str(j)).append(eval('corr_network' + str(j) + '_r2xr1'))

    # Heatmap
    #fig = plt.subplot(2, 4, j)
    #im = fig.imshow(eval('corr_network' + str(j) + '_r1xr2'), cmap='hot', interpolation='nearest')
    #plt.title(j)


    #Identification analyses
    #Monozygotic
    print('MZ identification accuracies for network'+str(j))
    for i in range(0,len(eval('corr_matrices_n' + str(j)))):
        eval('accuracies_n'+str(j)+'_MZ').append(accuracy_t_MZ(eval('corr_matrices_n' + str(j))[i][0:246,:],ID_individuals,tid.twin_pairs_ID))
    print eval('accuracies_n'+str(j)+'_MZ')
    #print np.mean(eval('accuracies_n'+str(j)+'_MZ'))
    #print np.std(eval('accuracies_n'+str(j)+'_MZ'))

    #Dizygotic
    print('DZ identification accuracies for network' + str(j))
    for i in range(0, len(eval('corr_matrices_n' + str(j)))):
        eval('accuracies_n' + str(j)+'_DZ').append(
            accuracy_t_DZ(eval('corr_matrices_n' + str(j))[i][246:,:], ID_individuals, tid.twin_pairs_ID[246:]))
    print eval('accuracies_n' + str(j) + '_DZ')
    #print np.mean(eval('accuracies_n' + str(j)+'_DZ'))
    #print np.std(eval('accuracies_n' + str(j)+'_DZ'))

    ####Permutation - random twin_pair_ID
    accuracy_random_MZ_n1=[]
    accuracy_random_MZ_n2 = []
    accuracy_random_MZ_n3 = []
    accuracy_random_MZ_n4 = []
    accuracy_random_MZ_n5 = []
    accuracy_random_MZ_n6 = []
    accuracy_random_MZ_n7 = []
    accuracy_random_MZ_n8 = []
    accuracy_random_DZ_n1=[]
    accuracy_random_DZ_n2 = []
    accuracy_random_DZ_n3 = []
    accuracy_random_DZ_n4 = []
    accuracy_random_DZ_n5 = []
    accuracy_random_DZ_n6 = []
    accuracy_random_DZ_n7 = []
    accuracy_random_DZ_n8 = []

    for k in range(0,1000):
        accuracies_MZ = []
        accuracies_DZ = []
        for i in range(0, len(eval('corr_matrices_n' + str(j)))):
            accuracies_MZ.append(accuracy_t_MZ(eval('corr_matrices_n' + str(j))[i][0:246,:], ID_individuals, np.ndarray.tolist(
    np.random.permutation(tid.twin_pairs_ID[:246]))))
        eval('accuracy_random_MZ_n' + str(j)).append(np.mean(accuracies_MZ))

        for i in range(0, len(corr_matrices_whole_DZ)):
            accuracies_DZ.append(accuracy_t_DZ(eval('corr_matrices_n' + str(j))[i][246:,:], ID_individuals, np.ndarray.tolist(
    np.random.permutation(tid.twin_pairs_ID[246:]))))
        eval('accuracy_random_DZ_n' + str(j)).append(np.mean(accuracies_DZ))

#plt.show()
