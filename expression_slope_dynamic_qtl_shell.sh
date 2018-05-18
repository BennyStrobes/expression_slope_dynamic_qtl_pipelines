#!/bin/bash
#SBATCH --time=20:00:00 --partition=broadwl --mem=14GB

input_data_file="$1"
covariate_method="$2"
parameter_string="$3"
qtl_results_dir="$4"
qtl_visualization_dir="$5"
standard_dynamic_qtl_results="$6"

#############################################
# Part 1: Run analysis for real data
#############################################
permute="False"
qtl_real_results_file=$qtl_results_dir$parameter_string"_permute_"$permute"_qtl_results.txt"
Rscript run_expression_slope_dynamic_qtl.R $input_data_file $covariate_method $qtl_real_results_file $permute


#############################################
# Part 2: Run analysis for permuted data
#############################################
permute="True"
qtl_perm_results_file=$qtl_results_dir$parameter_string"_permute_"$permute"_qtl_results.txt"
Rscript run_expression_slope_dynamic_qtl.R $input_data_file $covariate_method $qtl_perm_results_file $permute



#############################################
# Part 3: Visualize distributions of both real and permuted data
#############################################
Rscript visualize_dynamic_qtls.R $qtl_real_results_file $qtl_perm_results_file $parameter_string $qtl_visualization_dir




#############################################
# Part 4: Run emperical FDR correction
#############################################
# output file for eFDR analysis
efdr_file=$qtl_results_dir$parameter_string"_eFDR_results.txt"
# Run eFDR correction
Rscript eFDR_correction.R $qtl_real_results_file $qtl_perm_results_file $efdr_file


fdr_thresh=".05"
# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
significant_efdr_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_efdr_gene_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant_egenes.txt"
python assess_significance_efdr_approach.py $efdr_file $qtl_real_results_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh


fdr_thresh=".1"
# Assess genome wide significance of actual data based on the eFDR approach with FDR <= $fdr_thresh
# Output file for all significant variant gene pairs
significant_efdr_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant.txt"
# Output file for significant egenes and their strongest associated variant
significant_efdr_gene_results=$qtl_results_dir$parameter_string"_efdr_"$fdr_thresh"_significant_egenes.txt"
python assess_significance_efdr_approach.py $efdr_file $qtl_real_results_file $significant_efdr_results $significant_efdr_gene_results $fdr_thresh




#############################################
# Part 5: Compare results to standard dynamic qtl calling
#############################################
dynamic_qtl_results_comparison_file=$qtl_results_dir$parameter_string"_dynamic_qtl_results_comparison.txt"
python organize_dynamic_qtl_comparison.py $standard_dynamic_qtl_results $qtl_real_results_file $dynamic_qtl_results_comparison_file
echo "DONE"

Rscript visualize_dynamic_qtl_comparison.R $dynamic_qtl_results_comparison_file $qtl_visualization_dir$parameter_string
