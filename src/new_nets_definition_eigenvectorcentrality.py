import numpy as np
import h5py
import os
import bct
from operator import itemgetter
import functions_fernanda as fun
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import twin_IDs as tid


#Defition of the list of connectivity matrices
finn_matrices_r1=[]
finn_matrices_r2=[]

#Definition of accuracy function for individual identification
def accuracy(corr_matrix,list_of_ID):
    acc = 0
    for j in range(0, len(corr_matrix)):
            #print(np.argmax(corr_matrix[j]),list_of_ID.index(list_of_ID[j]))
            if np.argmax(corr_matrix[j]) == list_of_ID.index(list_of_ID[j]):
                acc = acc + 1.
                #print acc
    # len(corr)=number of individuals
    accuracy = (acc / len(corr_matrix[0])) * 100.
    # print len(corr_matrices[0])
    # print accuracy
    return accuracy

#Definition of accuracy function for twins identification
#Monozygotic
def accuracy_t_MZ(corr_matrices,list_subject_ID,list_of_twin):
    acc = 0
    identified_pair=[]
    for i in range(0, len(corr_matrices)):
        identified_pair.append([int(list_subject_ID[i]),int(list_subject_ID[np.argmax(corr_matrices[i])])])
        if identified_pair[i]==list_of_twin[i]:
           acc = acc + 1.
    accuracy_t = (acc / len(corr_matrices)) * 100.
    return accuracy_t
#Dizygotic
def accuracy_t_DZ(corr_matrices,list_subject_ID,list_of_twin):
    acc = 0
    identified_pair=[]
    for i in range(0, len(corr_matrices)):
        identified_pair.append([int(list_subject_ID[i+246]),int(list_subject_ID[np.argmax(corr_matrices[i])])])
        if identified_pair[i]==list_of_twin[i]:
           acc = acc + 1.
    accuracy_t = (acc / len(corr_matrices)) * 100.
    return accuracy_t



#Individuals_ID
for i in range(0,380):
    with open('./Dados/Rest2/list_of_ID_all_MZ_DZ_selected') as fp:
        lines = fp.read().split("\n")

ID_individuals = lines[0:len(lines) - 1]


#Twin ID
twin_ID=tid.twin_pairs_ID


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

number_of_individuals=len(finn_matrices_r1)


#Determination of the most central nodes
#Average connectivity matrix
ACM_r1=np.mean(finn_matrices_r1, axis=0)
ACM_r2=np.mean(finn_matrices_r2, axis=0)

#Absolute values ACM
abs_ACM_r1=np.abs(ACM_r1)
abs_ACM_r2=np.abs(ACM_r2)

#Eigenvector centrality
EVC_r1=bct.eigenvector_centrality_und(abs_ACM_r1)
EVC_r2=bct.eigenvector_centrality_und(abs_ACM_r2)

#Eigenvector + node-label
EVC_r1_labels=[]
EVC_r2_labels=[]
for i in range(0,len(EVC_r1)):
    EVC_r1_labels.append([EVC_r1[i], i + 1])
    EVC_r2_labels.append([EVC_r2[i], i + 1])

#Ordering nodes
EVC_r1_labels_descending=sorted(EVC_r1_labels, key=itemgetter(0), reverse=True)
EVC_r2_labels_descending=sorted(EVC_r2_labels, key=itemgetter(0), reverse=True)
     #For increasing numkber of nodes from least to most central
     #EVC_r1_labels_descending = sorted(EVC_r1_labels, key=itemgetter(0))
     #EVC_r2_labels_descending = sorted(EVC_r2_labels, key=itemgetter(0))

#Selecting nodes with highest eigenvector centrality -> New network definition -> Identification accuracy
net_acc_ind_r1xr2=[]
net_acc_ind_r2xr1=[]
net_acc_MZ_r1xr2=[]
net_acc_DZ_r1xr2=[]

#For functional networks with fixed size, use this code
for j in range(0,len(EVC_r1_labels_descending)-60): #In this case, the network size is of 60 nodes
#For increasing number of nodes, use the following line instead
#for j in range(0, len(EVC_r1_labels_descending)):

    new_network=[]
    matrices_net_r1=[]
    matrices_net_r2=[]

    #corr_network1_r1xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))
    #corr_network1_r2xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
    corr_network_r1xr2= np.zeros(shape=(number_of_individuals,number_of_individuals))
    corr_network_r2xr1= np.zeros(shape=(number_of_individuals,number_of_individuals))


    #Network definition
    for i in range(0,60):
        new_network.append(EVC_r1_labels_descending[i+j][1])
    #For increasing number of nodes, use the following lines instead
        #for i in range(0, j):
        #   new_network.append(EVC_r1_labels_descending[i][1])

    fun.network(new_network, finn_matrices_r1, matrices_net_r1)
    fun.network(new_network, finn_matrices_r2, matrices_net_r2)


    #Correlation among individuals connectivity matrices
    for i in range(0, len(corr_network_r1xr2)):
        corr_row = []
        for m in range(0, len(corr_network_r1xr2)):
            corr_row.append(pearsonr(matrices_net_r1[i].reshape(-1),
                                     matrices_net_r2[m].reshape(-1))[0])
        corr_network_r1xr2[i] = corr_row




    #Identification analyses
    #Individual identification
    net_acc_ind_r1xr2.append(accuracy(corr_network_r1xr2,ID_individuals))
    net_acc_ind_r2xr1.append(accuracy(corr_network_r2xr1,ID_individuals))

    #Twin identification
    np.fill_diagonal(corr_network_r1xr2, -1)
    np.fill_diagonal(corr_network_r2xr1, -1)
    # Monozygotic
    corr_r1xr2_MZ = corr_network_r1xr2[0:246, :]
    corr_r2xr1_MZ = corr_network_r2xr1[0:246, :]
    # Dizygotic
    corr_r1xr2_DZ = corr_network_r1xr2[246:380, :]
    corr_r2xr1_DZ = corr_network_r2xr1[246:380, :]

    net_acc_MZ_r1xr2.append(accuracy_t_MZ(corr_r1xr2_MZ,ID_individuals,tid.twin_pairs_ID))
    net_acc_DZ_r1xr2.append(accuracy_t_DZ(corr_r1xr2_DZ,ID_individuals,tid.twin_pairs_ID[246:]))

    #Plotting accuracy values while in the loop
    #plt.scatter(j, net_acc_ind_r1xr2[j])
    #plt.pause(0.05)

#plt.show()

