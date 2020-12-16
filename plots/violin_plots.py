import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

parcellation = ['shen', 'gordon']
nets = [9, 13]
title_shen = ['Whole-brain', 'Medial frontal', 'Frontoparietal',
              'Default mode', 'Subcortical-cerebellum', 'Motor', 'Visual I',
              'Visual II', 'Visual association']
title_gordon = ['Whole-brain', 'Default mode', 'Somato-sensory hand',
                'Somato-sensory mouth', 'Visual', 'Frontoparietal', 'Auditory',
                'Cingulo Parietal', 'Retrosplenial Temporal',
                'Cingulo Opercular', 'Ventral Attention', 'Salience',
                'Dorsal Attention']
nodes_shen = [268, 29, 34, 20, 90, 50, 18, 9, 18]
nodes_gordon = [333, 41, 38, 8, 39, 24, 24, 5, 8, 40, 23, 4, 32]

# Violin plots
falconer_h2_shen = {'Mean_MZ': [], 'Mean_DZ': [], 'Falconers formula':[],'n_nodes':[]}
falconer_h2_gordon = {'Mean_MZ': [], 'Mean_DZ': [], 'Falconers formula':[],'n_nodes':[]}
falconer_h2 = {'shen':[], 'gordon':[]}
labels = ['SI', 'MZ', 'DZ', 'UN']
for i in range(len(parcellation)):
    sns.set_style("darkgrid")
    if nets[i] == 9:
        fig = plt.figure(figsize=(15, 15))
    else:
        fig = plt.figure(figsize=(15, 25))

    for k in range(0, nets[i]):

        plot = [np.load(
            './../outputs/corr_distributions/list_of_corr_SI_' + parcellation[
                i] + '_n' + str(k) + '.npz')[
                    'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_MZ_' + parcellation[
                i] + '_n' + str(k) + '.npz')[
                    'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_DZ_' + parcellation[
                i] + '_n' + str(k) + '.npz')[
                    'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_nonrelated_' +
            parcellation[i] + '_n' + str(
                k) + '.npz')['list']]
        if nets[i] == 9:
            fig = plt.subplot(3, 3, k + 1)
            fig = sns.violinplot(data=plot, palette="Blues_r",
                                 inner='quart')
        else:
            fig = plt.subplot(5, 3, k + 1)
            fig = sns.violinplot(data=plot, palette="BuGn_r", inner='quart')
        plt.xticks(range(4), labels, rotation='vertical', fontsize=15)
        fig.set_title(eval('title_' + parcellation[i])[k], fontsize=20)
        # fig.grid(color='gray', linestyle='--', linewidth=0.3)
        fig.set_ylim([0.3, 1])
        plt.yticks(fontsize=12)
        plt.ylabel('Correlation', fontsize=15)

        eval('falconer_h2_' + parcellation[i])['Mean_MZ'].append(
            np.mean(plot[1]))
        eval('falconer_h2_' + parcellation[i])['Mean_DZ'].append(
            np.mean(plot[2]))
        eval('falconer_h2_' + parcellation[i])['Falconers formula'].append(
            2 * (np.mean(plot[1]) - np.mean(plot[2])))
        eval('falconer_h2_' + parcellation[i])['n_nodes'].append(
            eval('nodes_'+parcellation[i])[k])

        falconer_h2[parcellation[i]].append(
            2 * (np.mean(plot[1]) - np.mean(plot[2])))


    pd.DataFrame.from_dict(pd.DataFrame(eval('falconer_h2_' + parcellation[i]))).to_excel(
        './../outputs/Falconers_heritability_' + parcellation[i] + '.xlsx')

    # plt.savefig('corr_distributions_' + parcellation[i] + '.pdf')
    plt.show()


np.savez('falconers_h2.npz', dict = falconer_h2)