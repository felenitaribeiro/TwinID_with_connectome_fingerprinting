import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

parcellation = ['shen', 'gordon']
nets = [9, 13]

labels_dict_shen = {'n0': ['All, 268', 'Whole-brain'],
                    'n1': ['MF, 29', 'Medial frontal'],
                    'n2': ['FP, 34', 'Frontoparietal'],
                    'n3': ['DMN, 20', 'Default mode'],
                    'n4': ['SC, 90', 'Subcortical-cerebellum'],
                    'n5': ['M, 50', 'Motor'],
                    'n6': ['VI, 18', 'Visual I'],
                    'n7': ['VII, 9', 'Visual II'],
                    'n8': ['VA, 18', 'Visual association']}
nodes_shen = [268, 29, 34, 20, 90, 50, 18, 9, 18]

labels_dict_gordon = {'n0': ['All, 333', 'Whole-brain'],
                      'n1': ['DAN, 41', 'Default mode'],
                      'n2': ['SMh, 38', 'Somato-sensory hand'],
                      'n3': ['SMm, 8', 'Somato-sensory mouth'],
                      'n4': ['V, 39', 'Visual'],
                      'n5': ['FP, 24', 'Frontoparietal'],
                      'n6': ['Au, 24', 'Auditory'],
                      'n7': ['CP, 5', 'Cingulo Parietal'],
                      'n8': ['RT, 8', 'Retrosplenial Temporal'],
                      'n9': ['CO, 40', 'Cingulo Opercular'],
                      'n10': ['VAN, 23', 'Ventral Attention'],
                      'n11': ['S, 4', 'Salience'],
                      'n12': ['DAN, 32', 'Dorsal Attention'],
                      }
nodes_gordon = [333, 41, 38, 8, 39, 24, 24, 5, 8, 40, 23, 4, 32]

falconers_h2 = np.load('./../outputs/falconers_h2.npz')['dict']
ACE_h2_shen = np.load('./../outputs/heritability_ACE_shen.npz')['list']
# cholesky_common_h2_shen = np.load('./../outputs/heritability_shen_cholesky_common.npz')['list']
# cholesky_edge_h2_shen = np.load('./../outputs/heritability_shen_cholesky_edge.npz')['list']
ACE_h2_gordon = np.load('./../outputs/heritability_ACE_gordon.npz')['list']

for j in range(len(parcellation)):
    accuracies_ind = \
        np.load('./../outputs/accuracies_ind_id_' + parcellation[j] + '.npz')[
            'dict']
    accuracies_twin = \
        np.load('./../outputs/accuracies_twin_id_' + parcellation[j] + '.npz')[
            'dict']


    results = {'SI_acc_mean': [], 'SI_acc_std': [],
               'MZ_acc_mean': [], 'MZ_acc_std': [],
               'DZ_acc_mean': [], 'DZ_acc_std': [],
               'Title': [], 'Falconers formula':[],
               'ACE':[]
               }

    for i in range(nets[j]):
        results['SI_acc_mean'].append(
            np.mean(accuracies_ind.item()['n' + str(i) + '_SI']))
        results['SI_acc_std'].append(
            np.std(accuracies_ind.item()['n' + str(i) + '_SI']))
        results['MZ_acc_mean'].append(
            np.mean(accuracies_twin.item()['n' + str(i) + '_MZ']))
        results['MZ_acc_std'].append(
            np.std(accuracies_twin.item()['n' + str(i) + '_MZ']))
        results['DZ_acc_mean'].append(
            np.mean(accuracies_twin.item()['n' + str(i) + '_DZ']))
        results['DZ_acc_std'].append(
            np.std(accuracies_twin.item()['n' + str(i) + '_DZ']))
        results['Title'].append(eval('labels_dict_' + parcellation[j])['n' + str(i)][0])
        results['Falconers formula'].append(
            falconers_h2.item()[parcellation[j]][i])
        results['ACE'].append(eval('ACE_h2_'+parcellation[j])[i])


    # Excel file
    # pd.DataFrame.from_dict(results).to_excel('./../outputs/identification_results_' + parcellation[j] + '.xlsx')
    df = pd.DataFrame(results)
    df=df.sort_values(by=['SI_acc_mean'],ascending=False)



    print(stat.pearsonr(results['MZ_acc_mean'],
                        eval('ACE_h2_'+parcellation[j])))
    print(stat.pearsonr(results['MZ_acc_mean'],
                        falconers_h2.item()[parcellation[j]]))
    print(stat.pearsonr(falconers_h2.item()[parcellation[j]],
                        eval('nodes_'+parcellation[j])))
    labels = df['Title']

    # Figure
    if nets[j] == 9:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig, ax = plt.subplots(figsize=(14, 5))
    plt.xlabel('Functional networks, n of nodes', fontsize=20)
    plt.ylabel('Identification accuracy', fontsize=20)
    bar_width = 0.25
    plt.xticks(range(nets[j]), labels, rotation=45, fontsize=15)
    plt.bar(np.arange(nets[j]) - bar_width, df['SI_acc_mean'], bar_width, align='center',
            yerr=df['SI_acc_std'], error_kw=dict(elinewidth=2, ecolor='k'), color='k',
            label='Individual identification')
    plt.bar(np.arange(nets[j]), df['MZ_acc_mean'], bar_width, align='center',
            yerr=df['MZ_acc_std'],
            error_kw=dict(elinewidth=2, ecolor='k'), color='dimgray',
            label='Monozygotic twin identification')
    plt.bar(np.arange(nets[j]) + bar_width, df['DZ_acc_mean'], bar_width, align='center',
            yerr=df['DZ_acc_std'], error_kw=dict(elinewidth=2, ecolor='k'),
            color='darkgray', label='Dizygotic twin identification')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.ylim(0, 100)
    # plt.savefig('IDaccuracies_' + parcellation[j] + '.pdf')
    plt.show()
