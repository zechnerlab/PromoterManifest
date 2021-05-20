# PromoterManifest
This repository provides all Matlab code used for "Promoters adopt distinct dynamic manifestations depending on transcription factor context" by AS Hansen and C Zechner, Mol Syst Biol (2021).

# Folder structure
-Common: contains general helper functions and code.
-Data: contains Msn2 gene expression trajectories and artificially generated data.
-ModelCalibration: files related to the model calibration step.
-StateReconstruction: files related to the Bayesian state reconstruction analysis.
-Manifestations: scripts for plotting the manifestations reported in the paper.
-Tutorial: contains tutorial-like example code to plot arbitrary transcriptional features inferred from the data.
-ComplexManifestationExample: simulation code for a hypothetical toy example exhibiting context-dependent manifestations.

# Manifestation database
In the directory "Tutorial" we provide simple tutorial-like example code to plot the inferred transcriptional features for different promoters and conditions. The folder contains two files: 
1) PlotPromoterFeatures.m: this file plots two transcriptional features against each others for all single-pulse conditions and a particular promoter. The results show population averages and standard errors.
2) PlotTemporalReconstructions.m this file plots the temporal reconstructions of promoter state probabilities and transcription rates as a function of time for a single condition and promoter. The results show population averages and standard errors (when applicable).

# Reproducing manifestation figures
All manifestations reported in Hansen&Zechner can be reproduced using the scripts in the subfolder ‘Manifestations’. The naming of the scripts contain the figure numbers in which they are shown. Results correspond to averages over 5 independent, complete runs of the analysis pipeline. The calibrated models and state inference results can be found in the “ModelCalibration” and “StateReconstruction” folders, respectively. The subfolders containing the individual runs are named ‘results_1’, … ‘results_5’.

# Performing model calibration and state reconstruction
To perform a full analysis run, a few scripts should be performed in the following order:
1. CalibrateModels.m: Obtains estimates of the model parameters from 50% of the data using a moment-based inference scheme (Zechner et al., PNAS (2012)). The results will be stored in the folder ModelCalibration/results/
2. RunStateReconstruction.m: Takes the calibrated models and performs state reconstruction using a hybrid SMC method on the remaining 50% of cells (those not used for calibration). The results are stored in the folder StateReconstruction/results/
3. EvaluateStateReconstruction.m: Post-processes the reconstructions and extracts the promoter features.
4. (optional) AverageStateReconstructions.m: In case multiple runs are analyzed, this script can be used to average over multiple different runs. Each run must be in a separate subfolder in the folder ‘StateReconstruction’ (e.g., results_1, results_2, …). The averaged results are stored in the folder ‘StateReconstruction/resultsAvrg’. Subsequently, the results can be plotted using the scripts from the folder ‘Manifestations’ and/or ‘Tutorial’. For further information on the usage of the individual scripts, see comments written directly in the respective .m files. Keep in mind that one complete run of the analysis takes about 24hrs to complete on a 40 core shared memory machine.

# Additional analyses
1. ValidateModelCalibration.m: Cross-validation of the calibrated models based on first and second moments.
2. ValidateStateReconstruction.m: Cross-validation of the trajectory inference results based on first and second moments.
3. ValidateTrajectoryExcludion.m: Certain trajectories can cause the SMC scheme to become unstable, which are thus excluded. This script compares the moments of gene expression responses in the presence and absence of trajectory exclusion for each individual promoter and condition.

# External libraries / code
Our repository relies on the following freely available scripts:
-npermutek.m: 
Matt Fig (2020). N_PERMUTE_K (https://www.mathworks.com/matlabcentral/fileexchange/11462-n_permute_k), MATLAB Central File Exchange. 
