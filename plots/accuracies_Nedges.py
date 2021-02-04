import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

accuracies_ID_slices = np.load('./../outputs/accuracies_SI_slices.npz')['list']
accuracies_MZ_slices = np.load('./../outputs/accuracies_MZ_slices.npz')['list']
accuracies_DZ_slices = np.load('./../outputs/accuracies_DZ_slices.npz')['list']

iterations = 100
edges = [10, 50, 100, 500, 1000, 5000, 10000]
identication = ['SI', 'MZ', 'DZ']

accuracies_ID_slices = np.reshape(accuracies_ID_slices, (-1))
accuracies_MZ_slices = np.reshape(accuracies_MZ_slices, (-1))
accuracies_DZ_slices = np.reshape(accuracies_DZ_slices, (-1))

data = np.concatenate([[np.array(accuracies_DZ_slices).T, edges * iterations,
                        iterations * len(edges) * [identication[2]]],
                       [np.array(accuracies_ID_slices).T, edges * iterations,
                        iterations * len(edges) * [identication[0]]],
                       [np.array(accuracies_MZ_slices).T, edges * iterations,
                        iterations * len(edges) * [identication[1]]]
                       ],
                      axis=1)
data = pd.DataFrame(columns=['Accuracy', 'N_edges', 'Group'],
                    data=data.T)
data['Accuracy'] = data['Accuracy'].astype(float)
data['N_edges'] = data['N_edges'].astype(float)

sns.set_style("darkgrid")
fig = sns.lineplot(data=data, x="N_edges", y="Accuracy", hue="Group",
                   palette='colorblind', ci='sd')
fig.set(xscale="log")
plt.ylabel('Accuracy (%)')
plt.xlabel('Number of edges')
# plt.savefig('fingerprint_accuracyVsNEdges.pdf')
plt.show()
