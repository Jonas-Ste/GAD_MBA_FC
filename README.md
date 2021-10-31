# Analysis Scripts and Data for Steinhäuser et al. "Functional dissection of neural connectivity in generalized anxiety disorder using Bayesian and frequentist inference"
Jonas L. Steinhäuser, Adam R. Teed, Obada Al-Zoubi, René Hurlemann, Gang Chen, Sahib S. Khalsa

These datafiles and scripts reproduce the results reported in the main manuscript and supplement.

## 00_data
This folder contains two csv-files. One provides the (anonymized) demographical and clinical information of the 55 subjects from the analyzed study sample. The second one contains the correlation coefficients for all ROIs analyzed for each subject individually. 

## 01_demographics_and_headmotion
This script produces the demographical table describing the study sample. Further it contains information on average headmotion as evaluated during first-level preprocessing. 

## 02_NHST_analysis
The first script "NHST_test_apriori_hypotheses.R" tests a-priori defined hypotheses of functional correlativity between ROIs as described in the [pre-registration of this study](https://osf.io/j29qv). The second script "correlation_clinical_variables.R" produces the table with correlations between the vmPFC-PMI z-score and clinical variables that is included in the supplementary material.

## 03_MBA_analysis
This folder contains all files necessary to re-run the Bayesian Multilevel model. Please note that an up-to-date installation of AFNI is needed in order for the scripts to run.
Please refer to the [help file of the AFNI MBA function](https://afni.nimh.nih.gov/pub/dist/doc/program_help/MBA.html) on how to set-up MBA correctly.
Then use the "prepare_MBA_analysis.R" script to produce the input file to MBA in the correct structure. For convenience, the generated input file ("MBA_input_full.txt" is included in the folder but can also be reproduced by running the script. 
Lastly, the file "run_MBA_full.txt" contains the AFNI command to call the MBA function with the correct model and the data-table produced by the script above. 
This can be executed by running the following command in your console:
```
nohup tcsh -x run_MBA_full.txt > diary.txt &
```
For details please refer to the MBA help file linked above.
Please note that the final outputs will not exactly match the results reported in the manuscript and supplement due to nature of randomness involved in the simulations run in this script. They will, however, always be very close to each other.
