import os

with open('/home/hardleste/Desktop/HCP_DATA_12_2017/python_scripts_ESTUDO1/Rest2/list_of_ID') as fp:
    lines = fp.read().split("\n")

print lines
for i in range(0,50):
    os.rename('/home/hardleste/Desktop/HCP_DATA_12_2017/python_scripts_ESTUDO1/Rest2/Preprocessed_subjects/rest2_401_422/results/firstlevel/ANALYSIS_01/resultsROI_Subject'+str(i+1).zfill(3)+'_Condition001.mat','/home/hardleste/Desktop/HCP_DATA_12_2017/python_scripts_ESTUDO1/Rest2/Preprocessed_subjects/rest2_401_422/results/firstlevel/ANALYSIS_01/resultsROI_'+lines[i+395]+'_R2.mat')