import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

parcellation=['shen','gordon']
nets=[9,13]

# Create data
for j in range(len(parcellation)):
    accuracies_ind = np.load('./../outputs/accuracies_ind_id_' + parcellation[j] + '.npz')['dict']
    accuracies_twin = np.load('./../outputs/accuracies_twin_id_' + parcellation[j] + '.npz')['dict']
    effect_size = np.load('./../outputs/effect_sizes_' + parcellation[j] + '.npz')['dict']
    results={'acc_SI' : [], 'error_SI' : [], 'effect_size_SI' : [],
             'acc_MZ': [], 'error_MZ': [], 'effect_size_MZ': [],
             'acc_DZ': [], 'error_DZ': [], 'effect_size_DZ': []
             }

    for i in range(nets[j]):
        results['acc_SI'].append(np.mean(accuracies_ind.item()['n' + str(i) + '_SI']))
        results['error_SI'].append(np.std(accuracies_ind.item()['n' + str(i) + '_SI']))
        results['effect_size_SI'].append(np.float(effect_size.item()['n'+str(i)+'_SIvsUN']))

        results['acc_MZ'].append(np.mean(accuracies_twin.item()['n' + str(i) + '_MZ']))
        results['error_MZ'].append(np.std(accuracies_twin.item()['n' + str(i) + '_MZ']))
        results['effect_size_MZ'].append(np.float(effect_size.item()['n' + str(i) + '_MZvsUN']))

        results['acc_DZ'].append(np.mean(accuracies_twin.item()['n' + str(i) + '_DZ']))
        results['error_DZ'].append(np.std(accuracies_twin.item()['n' + str(i) + '_DZ']))
        results['effect_size_DZ'].append(np.float(effect_size.item()['n' + str(i) + '_DZvsUN']))


    data = np.concatenate([[results['effect_size_DZ'], results['acc_DZ'], len(results['effect_size_DZ']) * ['DZ']],
                           [results['effect_size_SI'], results['acc_SI'], len(results['effect_size_SI']) * ['SI']],
                           [results['effect_size_MZ'], results['acc_MZ'], len(results['effect_size_MZ']) * ['MZ']]],
                          axis=1)

    df = pd.DataFrame(columns=['Effect size', 'Accuracy', 'Group'], data=data.T)
    df['Effect size'] = df['Effect size'].astype(float)
    df['Accuracy'] = df['Accuracy'].astype(float)

    # Figure
    fig = plt.figure()
    sns.set_style("darkgrid")
    ax = sns.scatterplot(x='Effect size', y='Accuracy', hue='Group', data=df,
                         palette='colorblind')
    plt.ylim(0, 100)
    plt.xlim(0, 1)
    plt.ylabel('Prediction accuracy (%)', fontsize=15)
    plt.xlabel("Effect size (Cliff's delta)", fontsize=15)
    # plt.savefig('scatterplot_' + parcellation[j] + '.pdf')

    plt.show()
