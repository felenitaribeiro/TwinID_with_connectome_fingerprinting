import cliffsDelta.cliffsDelta as cd
import numpy as np
import pandas as pd

parcellation = ['shen', 'gordon']
nets = [9, 13]
for j in range(len(parcellation)):
    distributions = {}
    for k in range(0, nets[j]):
        distributions['n' + str(k)] = [np.load(
            './../outputs/corr_distributions/list_of_corr_SI_' + parcellation[j] + '_n' + str(k) + '.npz')[
                                           'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_MZ_' + parcellation[j] + '_n' + str(k) + '.npz')[
                                           'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_DZ_' + parcellation[j] + '_n' + str(k) + '.npz')[
                                           'list'], np.load(
            './../outputs/corr_distributions/list_of_corr_nonrelated_' + parcellation[j] + '_n' + str(
                k) + '.npz')['list']]

    cliffsdelta = {}
    results = {'cd_SIvsUN': [], 'cd_MZvsUN': [], 'cd_DZvsUN': []}

    for i in range(nets[j]):
        cliffsdelta['n' + str(i) + '_SIvsUN'] = cd.cliffsDelta(
            distributions['n' + str(i)][0], distributions['n' + str(i)][3])[0]
        cliffsdelta['n' + str(i) + '_MZvsUN'] = cd.cliffsDelta(
            distributions['n' + str(i)][1], distributions['n' + str(i)][3])[0]
        cliffsdelta['n' + str(i) + '_DZvsUN'] = cd.cliffsDelta(
            distributions['n' + str(i)][2], distributions['n' + str(i)][3])[0]

        results['cd_SIvsUN'].append(cliffsdelta['n' + str(i) + '_SIvsUN'])
        results['cd_MZvsUN'].append(cliffsdelta['n' + str(i) + '_MZvsUN'])
        results['cd_DZvsUN'].append(cliffsdelta['n' + str(i) + '_DZvsUN'])

    # np.savez('./../outputs/effect_sizes_' + parcellation[j] + '.npz', dict=cliffsdelta)

    # Excel file
    pd.DataFrame.from_dict(results).to_excel(
        './../outputs/cliffDelta_' + parcellation[j] + '.xlsx')
