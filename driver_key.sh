#!/bin/bash
#SBATCH --time=3:00:00 --partition=broadwl --mem=14GB

###############################################################################
# Input Data
###############################################################################

# Directory created by "time_step_independent_qtl_pipelines" scripts
# Contains 1 file per sample with information on each test (variant, target region)
# Each file (sample) has the same number of lines (tests)
cht_input_file_dir="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/cht_input_files/"

# File containing all of the target regions we are using. We are using this file to convert from gene positions to ensable id
target_region_input_file="/project2/gilad/bstrober/ipsc_differentiation/time_step_independent_qtl_pipelines/wasp/target_regions/target_regions_cis_distance_50000_maf_cutoff_0.1_min_reads_90_min_as_reads_20_min_het_counts_4_merged.txt"


# Files containing mapping from sample id to Nirmal's pseudotime predictions
# 3 state HMM
pseudotime_predictions_3_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/3_state_output_14.csv"
# 4 state HMM
pseudotime_predictions_4_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/4_state_output_14.csv"
# 5 state HMM
pseudotime_predictions_5_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess_input_data/pseudotime_predictions_14_lines/5_state_output_14.csv"


# cell line specific pcs
cell_line_specific_pc_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/covariates/cell_line_ignore_missing_principal_components_5.txt"

# Gene expression data for all samples
# Expression is RPKM transformed, then quantile normalized.
# Script will standardize each gene
total_expression_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/processed_total_expression/quantile_normalized_no_projection.txt"

# Dosage-based genotypes for all samples
genotype_file="/project2/gilad/bstrober/ipsc_differentiation/preprocess/genotype/YRI_genotype.vcf"

###############################################################################
# Output directories (aasume all of these exist prior to starting analysis)
###############################################################################

# Root directory for this of all ipsc data based results
output_root="/project2/gilad/bstrober/ipsc_differentiation/expression_slope_dynamic_qtl_pipelines//"

# Directory containing necessary input files to qtl tests
input_data_dir=$output_root"input_data/"


# Directory containing text files with results from dynamic qtl analysis
qtl_results_dir=$output_root"qtl_results/"

# Directory containing visualization of results found in qtl_results_dir
qtl_visualization_dir=$output_root"qtl_visualization/"





###############################################################################
# Dynamic QTL Calling
###############################################################################


##########################################
# Step 1: Create joint_test_input_file
##########################################
# joint_test_input_file is a table with 1 row per sample
# each row has 3 columns:
#### 1. Sample id
#### 2. environmental variable
#### 3. Absolute directory to CHT input file for that sample
# NOTE: This script is very specific to our data
# Takes less than a minute to run
# $environmental_variable is a parameter that describes how we parameterize the environmental variable. So far, this is done with:
### 1. 'time_steps': raw time format
### 2. 'pseudotime_predictions_3': pseudotime predictions (from Nirmal) using hmm with 3 latent variables
### 3. 'pseudotime_predictions_4': pseudotime predictions (from Nirmal) using hmm with 4 latent variables
### 4. 'pseudotime_predictions_5': pseudotime predictions (from Nirmal) using hmm with 5 latent variables
### 5. 'uniform_4': 
### 6. 'time_steps_pseudotime'
### 7. 'time_steps_endpoints'

environmental_variable_form="time_steps"
joint_test_input_file=$input_data_dir"joint_test_input_file_"$environmental_variable_form".txt"
if false; then
python create_joint_test_input_file.py $cht_input_file_dir $joint_test_input_file $environmental_variable_form $pseudotime_predictions_3_file $pseudotime_predictions_4_file $pseudotime_predictions_5_file $cell_line_specific_pc_file
fi

##########################################
# Step 2: Prepare input files for dynamic qtl calling
##########################################
# How to hanlde genotypes in the model. Options currently include:
## 1. "round"
## 2. "dosage"
genotype_version="dosage"

# Stem for all output files
input_data_file=$input_data_dir"expression_slope_dynamic_qtl_input_file_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_input.txt"
if false; then
python prepare_expression_slope_dynamic_qtl_input_files.py $joint_test_input_file $target_region_input_file $total_expression_file $genotype_file $genotype_version $input_data_file $environmental_variable_form
fi





##########################################
# Step 3: Run expression slope dynamic qtl calling
##########################################
## PC Options:
##### 1. "none"
##### 2. "pc1"
##### 3. "pc1_2"
##### 4. "pc1_3"
covariate_method="none"
parameter_string="expression_slope_dynamic_qtl_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_cov_method_"$covariate_method
sh expression_slope_dynamic_qtl_shell.sh $input_data_file $covariate_method $parameter_string $qtl_results_dir $qtl_visualization_dir

if false; then
covariate_method="pc1"
parameter_string="expression_slope_dynamic_qtl_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_cov_method_"$covariate_method
sbatch expression_slope_dynamic_qtl_shell.sh $input_data_file $covariate_method $parameter_string $qtl_results_dir $qtl_visualization_dir

covariate_method="pc1_2"
parameter_string="expression_slope_dynamic_qtl_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_cov_method_"$covariate_method
sbatch expression_slope_dynamic_qtl_shell.sh $input_data_file $covariate_method $parameter_string $qtl_results_dir $qtl_visualization_dir

covariate_method="pc1_3"
parameter_string="expression_slope_dynamic_qtl_environmental_variable_"$environmental_variable_form"_genotype_version_"$genotype_version"_cov_method_"$covariate_method
sbatch expression_slope_dynamic_qtl_shell.sh $input_data_file $covariate_method $parameter_string $qtl_results_dir $qtl_visualization_dir
fi






