# Analysis Scripts and Data for Steinhäuser et al. "Functional dissection of neural connectivity in generalized anxiety disorder"
Jonas L. Steinhäuser, Adam R. Teed, Obada Al-Zoubi, René Hurlemann, Gang Chen, Sahib S. Khalsa

These datafiles and scripts reproduce the results reported in the main manuscript.
The root of this project contains a file with the preprocessing parameters used in AFNI for preprocessing raw MRI data ("preprocessing_command_afniproc.sh"). This script is meant to transparently communicate choices during preprocessing of the data but can not be run since raw MRI data of subjects are not provided due to their size and the sensitive nature of this data. 
The files "FC_analysis_mask_resampled+tlrc.HEAD" and "FC_analysis_mask_resampled+tlrc.BRIK" are an AFNI mask in MNI-space that contain the ROIs used in this analysis and are based on the [Brainnetome-Atlas](https://atlas.brainnetome.org/) (see materials & methods section of the manuscript for the details). These can be used in AFNI to display the ROIs we used for the analysis - each ROI is coded separately and will appear in a different color in AFNI.

## 00_data
This folder contains two csv-files. One provides the (anonymized) demographical and clinical information of the 55 subjects from the analyzed study sample. The second one contains the correlation coefficients for all ROIs analyzed for each subject individually. 

## 01_demographics_and_headmotion
The script "compare_groups.R" produces the demographical table describing the study sample (Table 2 of the manuscript). Further it contains information on average headmotion as evaluated during first-level preprocessing. 

## 02_MBA_analysis
This folder contains all files necessary to re-run the Bayesian multilevel model. Please note that an up-to-date installation of AFNI is needed in order for the scripts to run.
Please refer to the [help file of the AFNI MBA function](https://afni.nimh.nih.gov/pub/dist/doc/program_help/MBA.html) on how to set-up AFNI and MBA correctly.
Then use the "prepare_MBA_analysis.R" script to produce the input file to MBA in the correct structure. For convenience, the generated input file ("MBA_input_full.txt") is included in the folder but can also be reproduced by running the script. The file "ROIlist.txt" contains the names for the ROIs as they will appear in the output of the MBA function.
Lastly, the file "run_MBA_full.txt" contains the AFNI command to call the MBA function with the correct model and the data-table produced by the script above. 
This can be executed by running the following command in your console:
```
nohup tcsh -x run_MBA_full.txt > diary.txt &
```
For details please refer to the MBA help file linked above.
The folder "exemplary_output" contains the the results from running the Bayesian multilevel model that are reported in the manuscript. The R-Workspace file containing the results is also used in the "04_Figures"-section to create Figure 1 and 2.
Please note that the final outputs from running the scripts above will not exactly match the results from the "exemplary_output"-folder due to nature of randomness involved in the simulations run in this script. They will, however, always be very close to each other.

## 03_NHST_analysis
The first script "NHST_test_apriori_hypotheses.R" tests a-priori defined hypotheses of functional correlativity (yields Table 3 of the manuscript) between ROIs as described in the [pre-registration of this study](https://osf.io/j29qv) or in Table 1 of the manuscript. The second script "correlation_clinical_variables.R" produces the table with correlations between the vmPFC-PMI z-score and clinical variables (Table 4 of the manuscript).

## 04_Figures
This folder contains all scripts and files necessary to produce Figures 1-4 from the manuscript as well as the "raw_figures"-subfolder that contains the raw figures in high-resolution PDF format generated from the files/scripts in this directory.
The scripts "Figure2_create_ridge_plot_region_pairs.R" and "Figure3_create_ridge_plot_region_effects.R" are courtesy of Gang Chen that have been modified to suit the data structure in this repository. They will create the ridge plots of region pairs and region effects (Figure 2 and 3). Please note that the labels of Figure 2 have been slightly modified (bold font and/or asterisks) to highlight region pairs included in the NHST-analysis.
The script "Figure4_create_raincloudplot.R" recreates Panel B of Figure 4, the raincloud plots of the vmPFC-PMI z-scores for each group individually.
The file "Figure4_Illustrator_File.ai" is a project-file for Adobe Illustrator used to create Figure 4 with panels A and B. Please note that Adobe Illustrator, a commercial software, is needed in order to view and edit this file. The file "Figure1_CONSORT_diagram.pptx" recreates the CONSORT diagram of the study.
