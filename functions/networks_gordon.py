import numpy as np
import xlrd

workbook = xlrd.open_workbook('./../data/gordon_333_parcellation_networklabels.xlsx')
sheet = workbook.sheet_by_index(0)

labels=[]
for rowx in range(sheet.nrows):
    cols = sheet.row_values(rowx)
    labels.append(cols)

#print len(labels)

#def networks
network1=[]
network2=[]
network3=[]
network4=[]
network5=[]
network6=[]
network7=[]
network8=[]
network9=[]
network10=[]
network11=[]
network12=[]
network13=[]


for i in range (0,len(labels)):
    if labels[i][1]=='Default':
        network1.append(labels[i][0])
    elif labels[i][1]=='SMhand':
        network2.append(labels[i][0])
    elif labels[i][1]=='SMmouth':
        network3.append(labels[i][0])
    elif labels[i][1]=='Visual':
        network4.append(labels[i][0])
    elif labels[i][1]=='FrontoParietal':
        network5.append(labels[i][0])
    elif labels[i][1]=='Auditory':
        network6.append(labels[i][0])
    elif labels[i][1]=='CinguloParietal':
        network7.append(labels[i][0])
    elif labels[i][1]=='RetrosplenialTemporal':
        network8.append(labels[i][0])
    elif labels[i][1]=='CinguloOperc':
        network9.append(labels[i][0])
    elif labels[i][1]=='VentralAttn':
        network10.append(labels[i][0])
    elif labels[i][1]=='Salience':
        network11.append(labels[i][0])
    elif labels[i][1]=='DorsalAttn':
        network12.append(labels[i][0])
    elif labels[i][1]=='None':
        network13.append(labels[i][0])

n_nodes=[len(network1),len(network2),len(network3),len(network4),len(network5),len(network6),len(network7),len(network8),len(network9),len(network10),len(network11),len(network12)]
