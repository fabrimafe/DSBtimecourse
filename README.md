# DSBtimecourse: Inference of kinetic parameters for time-courses sequenced with the UMI-DSBseq protocol 

Repository for "Uncovering the Dynamics of Precise Repair at CRISPR/Cas9-induced Double-Strand Breaks", Ben Tov*, Mafessoni* et al., biorXiv (2023). 
https://doi.org/10.1101/2023.01.10.523377 

For documentation of the R scripts, just use the --help flag.
.sh scripts are wrappers to run the R scripts on a cluster.
README*.txt document the development and testing stage for Ben Tov*, Mafessoni * et al.
The procedure for the modeling follows this general scheme:

### 1) data preparation
Prepare input files following prepare_data.R, and error matrices using calculate_error_matrix.R.
If necessary create input bootstraps files with create_stbootstrap.R. If for your data stratified bootstrapping is not possible due to the absence of repeated measures, alternative bootstrapping procedures are implemented in timeseriesbootstraps.R.

### 2) optimization
This is the core of the procedure, fitting the maximum likelihood parameters for the selected model. Use optimize_model_backbone.v1.R on the original dataset and on the bootstrapped data.

### 3) likelihood confidence intervals
Particularly when boostrapping is not an option, one can calculate likelihood-ratio based confidence intervals with calculate_CI.v1.R. While this would not be necessary per se, the output of this file can be used to generates plots in step 4)

### 4) plotting and summarizing
The program plot_bootstraps.R print flows and plots trajectories and induction curves. It can be run individually on a single time course or bootstrap or on mean-data+boostraps. In the former case likelihood-ratio based confidence intervals are used, in the latter bootstrap based confidence intervals.
