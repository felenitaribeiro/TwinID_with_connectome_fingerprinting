# Connectome fingerprint analyses

In here, I uploaded the main scripts for individual identification based on functional connectivity profiles that we wrote on Python.
All scripts are commented and further details of the analyses are described here (REF). Each file in briefly described below.

## Getting started

Scripts are briefly described below and if you any question concerning my code (which can a bit messy), just contact me! I will try to describe then in order of 'creation', so that will be easier to follow my thought line when coding.

### Rename_subjects:
As we pre-processed our data in batches (50 subjects per round), and in each round subjects were named from subject001 to subject050 by CONN, we had to rename our subject properly by using ID number. So, this is the little code that we used to rename the files.

### Parcellation_files:
CONN toolbox generates a .mat file containing correlation scores between fMRI signals from pairs of brain regions, which were defined by different parcellation schemas, for each subject.
Thus, this code segments the .mat file from each subject into different .csv files, each containing the functional connectivity matrix determined by a specific parcellation schema. We used 2 parcellation schemas: Finn (parcellation schema used by Finn et al. (2015), previously described by Shen et al. (2013)) and Gordon (Gordon et al. (2014)).

### Networks_finn: 
Script to classify nodes based on their attribution to specific functional networks a priori defined.

### Functions_fernanda:
Here it defined some functions that will be used in the following codes.
Fingerprint_ind_identification:
Code for individual identification based on functional connectivity profiles previously described (Finn et al., 2015). Identification analyses are performed considering the whole-brain as well as each functional network a priori defined. For individual predictions, mean identification accuracy is reported from the comparison across databases (REST1xREST2 and REST2xREST1). Also, correlation score were grouped based the relationship between pair of individuals for further analyses.

### Twin_identification:
Previous analysis was extend for twin identification. In here, both monozygotic and dizygotic twin identification were performed considering the whole-brain as well as each functional network.
