#!/usr/bin/env python3
"""
Script to transform CONN outputs to CSV files according the parcellation scheme.
"""
import scipy.io
import numpy as np
import os

for i in range(0, 246):
    with open('./Dados/Rest2/list_of_ID_MZ_R1R2') as fp:
        lines = fp.read().split("\n")

    lines = lines[0:len(lines) - 1]

    # Loading CONN output files
    mat = scipy.io.loadmat('./Dados/Rest2/All_MZ_R2/resultsROI_' + lines[i] + '_R2.mat')
    mat_conn = mat["Z"]
    mat_name_lines = mat["names"]
    mat_name_colunm = mat["names2"]

    # Slicing matrix according to the parcellation scheme
    finn_slice = slice(133, 401)
    gordon_slice = slice(401, 734)
    networks_slice = slice(734, 766)

    mat_conn_finn = mat_conn[finn_slice, finn_slice]
    mat_conn_gordon = mat_conn[gordon_slice, gordon_slice]
    mat_conn_networks = mat_conn[networks_slice, networks_slice]

    # Filling diagonal of connectivity matrix with 1s
    np.fill_diagonal(mat_conn_finn, 1)
    np.fill_diagonal(mat_conn_gordon, 1)
    np.fill_diagonal(mat_conn_networks, 1)

    # Creating directory for each subject
    if not os.path.exists("./Dados/Rest2/All_MZ_R2/" + lines[i]):
        os.makedirs("./Dados/Rest2/All_MZ_R2/" + lines[i])

    # Creating CSV files with functional connectivity matrix of each subject with each parcellation scheme
    np.savetxt("./Dados/Rest2/All_MZ_R2/" + lines[i] + "/mat_conn_finn_r2.csv", mat_conn_finn, delimiter=",")
    np.savetxt("./Dados/Rest2/All_MZ_R2/" + lines[i] + "/mat_conn_gordon_r2.csv", mat_conn_gordon, delimiter=",")
    np.savetxt("./Dados/Rest2/All_MZ_R2/" + lines[i] + "/mat_conn_networks_r2.csv", mat_conn_networks, delimiter=",")
