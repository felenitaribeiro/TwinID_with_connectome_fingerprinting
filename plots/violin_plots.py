import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt
import seaborn as sns

parcellation=['shen','gordon']
nets=[9,13]
title_shen = ['Whole-brain', 'Medial frontal', 'Frontoparietal',
              'Default mode', 'Subcortical-cerebellum', 'Motor', 'Visual I',
              'Visual II', 'Visual association']
title_gordon = ['Whole-brain', 'Default mode', 'Somato-sensory hand',
                'Somato-sensory mouth', 'Visual', 'Frontoparietal', 'Auditory',
                'Cingulo Parietal', 'Retrosplenial Temporal',
                'Cingulo Opercular', 'Ventral Attention', 'Salience',
                'Dorsal Attention']


# Violin plots
labels = ['SI', 'MZ', 'DZ', 'UN']
for i in range(len(parcellation)):
    sns.set_style("darkgrid")
    if nets[i] == 9:
        fig = plt.figure(figsize=(15, 15))
    else:
        fig = plt.figure(figsize=(15, 25))

    for k in range(0, nets[i]):

        plot = [np.load(
            './../outputs/corr_distributions/list_of_corr_SI_' + parcellation[i] + '_n' + str(k) + '.npz')[
                         'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_MZ_' + parcellation[i] + '_n' + str(k) + '.npz')[
                         'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_DZ_' + parcellation[i] + '_n' + str(k) + '.npz')[
                         'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_nonrelated_' + parcellation[i] + '_n' + str(
                k) + '.npz')['list']]
        if nets[i]==9:
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

    plt.show()



