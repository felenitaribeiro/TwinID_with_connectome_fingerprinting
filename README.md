# Connectome fingerprint analyses

In here, I uploaded the main scripts for individual identification based on functional connectivity profiles that we wrote on Python.
All scripts are commented and further details of the analyses are described here (REF). Each file in briefly described below.

## Getting started

Scripts are briefly described below and if you any question concerning my code (which can a bit messy :stuck_out_tongue_closed_eyes:), just contact me! I will try to describe then in order of 'creation', so that will be easier to follow my thought line when coding.

### First step codes:
#### Rename_subjects
As we pre-processed our data in batches (50 subjects per round), and in each round subjects were named from subject001 to subject050 by CONN, we had to rename our subject properly by using ID number. So, this is the little code that we used to rename the files.

#### Parcellation_files
CONN toolbox generates a .mat file containing correlation scores between fMRI signals from pairs of brain regions, which were defined by different parcellation schemas, for each subject.
Thus, this code segments the .mat file from each subject into different .csv files, each containing the functional connectivity matrix determined by a specific parcellation schema. We used 2 parcellation schemas: Finn (parcellation schema used by Finn et al. (2015), previously described by Shen et al. (2013)) and Gordon (Gordon et al. (2014)).

#### Networks_finn 
Script to classify nodes based on their attribution to specific functional networks a priori defined.

#### Functions_fernanda
Here it defined some functions that will be used in the following codes.

#### Twin_IDs
Little code to transfor a list of twin identities on excel into a list on Python. Note that, each line correspond to a pair of identities (paired twins). For our study, in this list, we sorted (first column) the monozygotic twins identities from the smallest to the largest ID number (corresponding to the identities from 1 to 246) and, in the following, we sorted in the same way the dizygotic twin identitis (corresponding to the identities from 247 to 380).

### Analyses' codes:
#### Fingerprint_ind_identification
Code for individual identification based on functional connectivity profiles previously described (Finn et al., 2015). Identification analyses are performed considering the whole-brain as well as each functional network a priori defined. For individual predictions, mean identification accuracy is reported from the comparison across databases (REST1xREST2 and REST2xREST1). Also, correlation score were grouped based the relationship between pair of individuals for further analyses.

#### Twin_identification:
Previous analysis was extend for twin identification. In here, both monozygotic and dizygotic twin identification were performed considering the whole-brain as well as each functional network.

### Pretty plots codes:

## References
Finn, E. S. et al. Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. Nat. Neurosci. 18, 1664–1671 (2015).

Gordon, E. M. et al. Generation and Evaluation of a Cortical Area Parcellation from Resting-State Correlations. Cereb. Cortex 26, 288–303 (2014).

My own paper, which will be available soon.

Shen, X., Tokoglu, F., Papademetris, X. & Constable, R. T. Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. Neuroimage 82, 403–415 (2013).
