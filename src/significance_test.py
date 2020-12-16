import numpy as np
from scipy.stats import pearsonr
import functions.networks_shen as ns
import functions.networks_gordon as ng

parcellation = ['shen', 'gordon']
nets = [9, 13]

labels_dict_shen = {'n0': ['All', 268, 'Whole-brain'],
                    'n1': ['MF', 29, 'Medial frontal'],
                    'n2': ['FP', 34, 'Frontoparietal'],
                    'n3': ['DMN', 20, 'Default mode'],
                    'n4': ['SC', 90, 'Subcortical-cerebellum'],
                    'n5': ['M', 50, 'Motor'],
                    'n6': ['VI', 18, 'Visual I'],
                    'n7': ['VII', 9, 'Visual II'],
                    'n8': ['VA', 18, 'Visual association']}

labels_dict_gordon = {'n0': ['All', 333, 'Whole-brain'],
                      'n1': ['DAN', 41, 'Default mode'],
                      'n2': ['SMh', 38, 'Somato-sensory hand'],
                      'n3': ['SMm', 8, 'Somato-sensory mouth'],
                      'n4': ['V', 39, 'Visual'],
                      'n5': ['FP', 24, 'Frontoparietal'],
                      'n6': ['Au', 24, 'Auditory'],
                      'n7': ['CP', 5, 'Cingulo Parietal'],
                      'n8': ['RT', 8, 'Retrosplenial Temporal'],
                      'n9': ['CO', 40, 'Cingulo Opercular'],
                      'n10': ['VAN', 23, 'Ventral Attention'],
                      'n11': ['S', 4, 'Salience'],
                      'n12': ['DAN', 32, 'Dorsal Attention'],
                      }
n_nodes_shen = [268] + ns.n_nodes
n_nodes_gordon = [333] + ng.n_nodes

significance = {}
corr = {}
for j in range(len(parcellation)):

    accuracies_ind = \
        np.load('./../outputs/accuracies_ind_id_' + parcellation[j] + '.npz')[
            'dict']
    accuracies_random_ind = \
        np.load(
            './../outputs/accuracies_perm_ind_id_' + parcellation[j] + '.npz')[
            'dict']

    accuracies_twin = \
        np.load('./../outputs/accuracies_twin_id_' + parcellation[j] + '.npz')[
            'dict']
    accuracies_random_twin = \
        np.load(
            './../outputs/accuracies_perm_twin_id_' + parcellation[j] + '.npz')[
            'dict']

    SI_acc = []
    SI_acc_random = []
    MZ_acc = []
    MZ_acc_random = []
    DZ_acc = []
    DZ_acc_random = []

    for i in range(nets[j]):
        SI_acc.append(np.mean(accuracies_ind.item()['n' + str(i) + '_SI']))
        SI_acc_random.append(
            accuracies_random_ind.item()['n' + str(i) + '_SI'])
        MZ_acc.append(np.mean(accuracies_twin.item()['n' + str(i) + '_MZ']))
        MZ_acc_random.append(
            accuracies_random_twin.item()['n' + str(i) + '_MZ'])
        DZ_acc.append(np.mean(accuracies_twin.item()['n' + str(i) + '_DZ']))
        DZ_acc_random.append(
            accuracies_random_twin.item()['n' + str(i) + '_DZ'])


        #print('Max n' + str(i) + 'random id accuracy (SI):' + str(np.max(SI_acc_random[i])))

        significance['n' + str(i) + '_SI_' + parcellation[j]] = np.sum(
            SI_acc_random[i] >
            SI_acc[i]) / 1000
        significance['n' + str(i) + '_MZ_' + parcellation[j]] = np.sum(
            MZ_acc_random[i] >
            MZ_acc[i]) / 1000
        significance['n' + str(i) + '_DZ_' + parcellation[j]] = np.sum(
            DZ_acc_random[i] >
            DZ_acc[i]) / 1000


    # Correlation between identification accuracy and number of nodes
    corr['SI_' + parcellation[j]] = pearsonr(SI_acc, eval(
        'n_nodes_' + parcellation[j]))
    corr['MZ_' + parcellation[j]] = pearsonr(MZ_acc,
                                             eval(
                                                 'n_nodes_' + parcellation[
                                                     j]))
    corr['DZ_' + parcellation[j]] = pearsonr(DZ_acc,
                                             eval(
                                                 'n_nodes_' + parcellation[
                                                     j]))

# np.savez('./../outputs/corr_idacc_nodes.npz', corr)
# np.savez('./../outputs/significance_idacc.npz', significance)
