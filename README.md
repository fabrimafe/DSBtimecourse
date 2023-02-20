# DSBtimecourse: Inference of kinetic parameters for time-courses of DNA Double-Strand Breaks sequenced with the UMI-DSBseq protocol 

Repository for "Uncovering the Dynamics of Precise Repair at CRISPR/Cas9-induced Double-Strand Breaks", Ben Tov*, Mafessoni* et al., biorXiv (2023). 
https://doi.org/10.1101/2023.01.10.523377 

For documentation of the R scripts, just use the --help flag.
.sh scripts are wrappers to run the R commands on a cluster.
README*.txt document the development and testing stage for Ben Tov*, Mafessoni * et al.
The procedure for the modeling follows this general scheme:

### 1) data preparation
Prepare input files: a time-course-file with the abundance of the different types of molecules, and an error-matrix file, specifying the expected error rates. Example files can be found in test/. The files should be tab separated files. The time-course file has five columns, the first describing the time and the others the count of molecules for each type. For example:

|time|    y1|      y2|      y3|      y4|
|----|------|--------|--------|--------|
|0   | 3095 |   0    |   62   |   9    |
|0   | 1608 |   0    |   34   |   16   |
|6   | 2289 |  29    |   41   |   21   |
|6   | 3259 |  34    |  162   |   33   |

....

where y1 indicates the intact molecules, y2 the precise DSBs, y3 processed DSBs and y4 indels. Several replicates can be provided with a same time.
Error-matrix files are four-entries tab separated matrices specifying the probability that a molecule (in row) is observed as such or as other types. For example:

|from\observed as|intact           |precise DSB         |processed DSB       |indels             |
|----------------|-----------------|--------------------|--------------------|-------------------|
|intact          |0.985826643022484|0.00901705556457197 |0.000276158814197761|0.00488014259874579|
|precise DSBs    |0                |   0.994843698587056|0.000276158814197761|0.00488014259874579|
|processed DSBs  |0                |0                   |1                   |0                  |
|indels          |0                |0                   |0                   |1                  |

indicates that intact molecules have a ~0.5% chance of being classified as indels. Such files can be prepared following prepare_data.R, and error matrices using calculate_error_matrix.R. If necessary create input bootstraps files with create_stbootstrap.R. If for your data stratified bootstrapping is not possible due to the absence of repeated measures, alternative bootstrapping procedures are implemented in timeseriesbootstraps.R.

### 2) optimization
This is the core of the procedure, fitting the maximum likelihood parameters for the selected model. Use DSBtimecourse_optimizer.R on the original dataset and on the bootstrapped data. To run on the test dataset using the 4 state model in Ben Tov et al.,2023:
```
./DSBtimecourse_optimizer.R -T test/timecourse_RNP_Psy1.txt -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_nok12 -o test/output_4states.tsv -z 3 -n 20
 The 3-state model can be run with:
./DSBtimecourse_optimizer.R -T test/timecourse_RNP_Psy1.txt -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_3states.tsv -z 3 -n 20
```
where the arguments -z sets the shape of the induction curve to have 3 degrees of freedom, and n is the number of independent iterations. The default is 100 iterations, but more are recommended, especially for more complex models like the 4-states. In Ben Tov et al. (2023) we used 50000 iterations, though a smaller number is usually sufficient. 

### 3) likelihood confidence intervals
Particularly when boostrapping is not an option, one can calculate likelihood-ratio based confidence intervals with calculate_CI.v1.R. While this would not be necessary per se, the output of this file can be used to generates plots in step 4). Here, the output of the optimization step is given as an input with the flag -i, together with the original timecourse and the error matrix. 

```
./calculate_CI.R -i test/output_3states.tsv -d test/timecourse_RNP_Psy1.txt -o output.CI -z 3 -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_3states_CI.tsv
```
The previous steps can be iterated for independent bootstraps of the data.

### 4) plotting and summarizing
The program plot_bootstraps.R print flows and plots trajectories and induction curves. It can be run individually on a single time course or bootstrap or on mean-data+boostraps. In the former case likelihood-ratio based confidence intervals are used, in the latter bootstrap based confidence intervals.
