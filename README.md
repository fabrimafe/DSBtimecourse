# DSBtimecourse: Inference of kinetic parameters for time-courses of DNA Double-Strand Breaks sequenced with the UMI-DSBseq protocol 

Repository for "Uncovering the Dynamics of Precise Repair at CRISPR/Cas9-induced Double-Strand Breaks", Ben Tov*, Mafessoni* et al., biorXiv (2023). 
https://doi.org/10.1101/2023.01.10.523377 

For documentation of the R scripts, you can just use the --help flag. To run the code it is adviced to have R v4.0.1 or higher installed. Here you can also find the folders:
* **scripts**: shell wrappers to run the R commands on a IBM LSF cluster
* **logs**: README*.txt files documenting the development and testing stage for Ben Tov*, Mafessoni * et al.
* **test**: files used to simulate data and test the program.

To run estimates of Double Strand Breaks (DSB) repair dynamics on your data follow the procedure described below.

### 1) data preparation
The necessary input files are a time-course-file with the abundance of the different types of molecules, and an error-matrix file, specifying the expected error rates. The latter can be estimated from control data, in which UMIDSB-seq is applied without any induced DSB, as described below. Example files can be found in test/. The files should be tab separated files. The time-course file has five columns, the first describing the time and the others the count of molecules for each type. For example:

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

indicates that intact molecules have a ~0.5% chance of being classified as indels. Such files can be prepared following prepare_data.R, and error matrices using calculate_error_matrix.R. If necessary create input bootstraps files with create_bootstrap.R. If for your data stratified bootstrapping is not possible due to the absence of repeated measures, alternative bootstrapping procedures are implemented in timeseriesbootstraps.R.

A possible way to create example data is to **simulate** them using the script DSBtimecourse_simulate.R. This script takes as input a set of parameters (specified in the format of file test/params_modelDSBs1i1_3x4_k0.05_r0.01_induction.CI and fed as input flag -p) and a set of times specifying total number of reads and time of sampling (example file test/timecourse_n2k_72h.txt), a model (with flag -m, as example a 3-state model "modelDSBs1i1_3x4") and an output file with flag -o.
Example:
```
./DSBtimecourse_simulate.R -T test/timecourse_n2k_72h.txt -p test/params_modelDSBs1i1_3x4_k0.05_r0.01_induction.CI -m modelDSBs1i1_3x4 -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -o test/simulateddata.tsv
```
Here test/simulateddata.tsv can be used for subsequent analyses.

### 2) optimization
This is the core of the procedure, which consists in fitting the maximum likelihood parameters for the selected model. Use DSBtimecourse_optimizer.R on the original dataset and on the bootstrapped data. To run on the test dataset using the 4 state model in Ben Tov et al.,2023:
```
./DSBtimecourse_optimizer.R -T test/timecourse_RNP_Psy1.txt -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_nok12 -o test/output_4states.tsv -z 3 -n 20
```
 The 3-state model can be run with:
```
./DSBtimecourse_optimizer.R -T test/timecourse_RNP_Psy1.txt -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_3states.tsv -z 3 -n 20
```
where the arguments -z sets the shape of the induction curve to have 3 degrees of freedom. Since the algorithm relies on a numerical optimization that does not guarantee to find the global maximum, it is advised to run it with different initial combinations of parameters. This can be set with the flag -n, which specifies the number of independent random starting points to ensure that a global maximum is found. The default is 100 iterations, but more are recommended for more complex models like the 4-states. In Ben Tov et al. (2023) we used 5000 iterations, though a smaller number is usually sufficient. 
To run the optimization on the simulated data for the 3-state model, run
```
./DSBtimecourse_optimizer.R -T test/simulateddata.tsv -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_simulated_3states.tsv -z 3 -n 20
```
while optimization for the 4-state model can be run as:
```
./DSBtimecourse_optimizer.R -T test/simulateddata.tsv -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_nok12 -o test/output_simulated_4states.tsv -z 3 -n 20
```

If an induction curve is estimated or known a priori, for example through FACS or lucifase data, optimization can be constrained to the estimated curve. This can be done using the -u flag. With this flag, one can specify induction parameters in a file with the same format as the file used to simulate time-courses with known rate parameters. For example:
```
./DSBtimecourse_optimizer.R -T test/simulateddata.tsv -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_simulated_3states_constrained.tsv -z 3 -n 20 -u test/params_FACStable.5.txt
```
where the file test/params_FACStable.5.txt describes an induction curve estimated for tomato protoplast one the basis of FACS data.

### 3) likelihood confidence intervals
Two different strategies can be used to estimate confidence intervals. When boostrapping is not an option, one can calculate likelihood-ratio based confidence intervals with calculate_CI.v1.R. While this would not be necessary per se, the output of this file can be used to generates plots in step 4). Here, the output of the optimization step is given as an input with the flag -i, together with the original timecourse and the error matrix. 

```
./calculate_CI.R -i test/output_simulated_3states.tsv -d test/timecourse_RNP_Psy1.txt -z 3 -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -o test/output_simulated_3states_CI.tsv 
```

This provides two output files. An .RData file containing all the temporary object obtained during the confidence intervals computation; and a summary file, which takes the names provided in the -o flag, and which shows a table with 4 fields:
- max: the maximum likelihood estimate
- CIlow: the lowest 95% confidence interval obtained with a Monte-Carlo approximation of the asymptotic likelihood ratio based confidence interval.
- CIhigh: the highest 95% confidence interval obtained with a Monte-Carlo approximation of the asymptotic likelihood ratio based confidence interval.
- rate: the rate for which the estimate is reported. These files CI can also be used to simulate new data (for example to generate parametric bootstraps).

The previous steps can be iterated for independent bootstraps of the data. Bootstraps of the data can be generated using timeseriesbootstraps.R. Use ./timeseriesbootstraps.R --help to see arguments and options, as different bootstrapping strategies are implemented (maximum entropy bootstrap, stationary bootstrap). To generate a stratified bootstrap (i.e. data points are resampled within a single time-point; possible only when repeated measures are available):
```
./timeseriesbootstraps.R -i test/timecourse_RNP_Psy1.txt -o test/bootstraps -n 10
```
Note that calculating uncertainty through a bootstrap is advisable as it accounts additional uncertainty that is not simply due to the random sampling of reads compared to the likelihood-based method, i.e. overdispersion in the binomial sampling of reads. Use --help to get information on additional bootstrapping procedures, which can be specified with flag -m.

### 4) plotting and summarizing
The program plot_bootstraps.R print flows and plots trajectories and induction curves. It can be run individually on a single time course or bootstrap or on mean-data+boostraps. In the former case likelihood-ratio based confidence intervals are used, in the latter bootstrap based confidence intervals. Use
```
./plot_bootstraps.R --help
```
for information on the usage. The flag -i specifies the input data. If this is a .RData, this is assumed to be a single raw confidence interval file computed with calculate_CI.R. Otherwise, a text file containing the path of individual bootstraps is expected, and formatted as described in the --help. In the former case, if run on the test data:
```
./plot_bootstraps.R -i test/output_3states_CI.tsv.RData -d test/timecourse_RNP_Psy1.txt -E test/error_matrix4_Psy1_errorsfromunbroken.tsv -m modelDSBs1i1_3x4 -n 100 -z 3 -o test/plot3states
```

which will generate two files: 

- plot3states_plot.trajectories.pdf: a plot of the fitted trajectories of the timecourse
- plot3states_plot.induction.pdf: a plot of the induction curve
- plot3states_plot.flow.tab: a table showing AIC values and the flow (in terms of proportion of initial molecules) for each rate. Note that proportions can be higher than 1 when molecules go through a given process more than once.

### Notation

Rates are referred to with the following nomenclature, with names in brackets corresponding to the terminology used in Ben Tov et al.,2023:
- k11 (K<sub>cut</sub>): cutting rate from intact to precise DSBs (or any DSBs if no difference is made between precise DSBs and processed DSBs, i.e. 3 state model);
- k12: cutting rate from intact to processed DSBs; this parameter is not used in Ben Tov*, Mafessoni* et al.,2023, but can be implemented if researchers are concerned about off-site CRISPR/Cas9 cutting;
- rr12 (K<sub>processing</sub>): processing rate from precise DSBs to processed DSBs; 
- r11 (P<sub>direct</sub>): repair from precise DSBs to intact molecules (or from any DSBs if no difference is made between precise DSBs and processed DSBs, i.e. 3 state model), i.e. *precise repair* or *precise repair from precise DSBs*;
- r12 (E<sub>direct</sub>): repair from precise DSBs to indels, i.e. *repair-error from precised DSBs*;
- r21 (P<sub>processed</sub>): repair from processed DSBs to intact molecules;
- r22: (E<sub>processed</sub>): repair from processed DSBs to indels;
- r0 (r): speed of induction;
- K (U): maximum level of induction (untransfected fraction);
- r2 (decay): exponential decay in cutting rate over time;

### Other models

The scripts can also be used to test additional models in respect to the 3-states and 4-states models used in Ben Tov et al.,2023, using the -m flag of DSBtimecourse_optimizer. These are still experimental, but feel free to ask for help if you want help implementing a tailored model for your use case. Specifically:
- *modelDSBs1i1_mini*: 4 state model with same degrees of freedom as the 3 states, in which indels are assumed to come from processed DSBs and precise repair through precise DSB, i.e. r12=0 and r21=0  
- *modelDSBs1i1_fullimprecise*: identical to the 4-state model except that also processed DSBs are not assumed to come from processing, but also from imprecise cutting of intact moleculels, i.e. k12 is estimated and not fixed to 0.
- We immplemented and tested other models, which are at the moment not used for Ben Tov*, Mafessoni et al.,2023. **If you have specific requests, let me know** and I'd be happy to help or add alternative models that could be useful for you.


