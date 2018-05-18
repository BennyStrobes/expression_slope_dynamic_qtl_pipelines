import numpy as np
import os
import sys
import pdb

# Create mapping from variant-gene pair to standard dynamic qtl pvalue
def extract_standard_dynamic_pvalues(standard_dynamic_qtl_results):
    f = open(standard_dynamic_qtl_results)
    head_count = 0
    pvalue_dict = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        # Skip header
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Extract relevent fields for this test
        rs_id = data[2]
        ensamble_id = data[5]
        test_name = rs_id + '_' + ensamble_id
        pvalue = float(data[-3])
        # Add test to dictionary
        pvalue_dict[test_name] = pvalue
    return pvalue_dict

# Create mapping from variant-gene pair to expression slope dynamic qtl pvalue
def extract_expression_slope_dynamic_pvalues(expression_slope_dynamic_qtl_results):
    f = open(expression_slope_dynamic_qtl_results)
    pvalue_dict = {}
    for line in f:
        line = line.rstrip()
        data = line.split()
        # Extract relevent fields
        rs_id = data[0]
        ensamble_id = data[1]
        pvalue = float(data[-1])
        test_name = rs_id + '_' + ensamble_id
        # Add test to dictionary
        pvalue_dict[test_name] = pvalue
    return pvalue_dict


standard_dynamic_qtl_results = sys.argv[1]
expression_slope_dynamic_qtl_results = sys.argv[2]
output_file = sys.argv[3]

# Create mapping from variant-gene pair to standard dynamic qtl pvalue
standard_dynamic_pvalues = extract_standard_dynamic_pvalues(standard_dynamic_qtl_results)
# Create mapping from variant-gene pair to expression slope dynamic qtl pvalue
expression_slope_dynamic_qtl_pvalues = extract_expression_slope_dynamic_pvalues(expression_slope_dynamic_qtl_results)


# Print results to output file
t = open(output_file, 'w')
# Print header
t.write('rs_id\tensamble_id\tstandard_dynamic_qtl_pvalue\texpression_slope_dynamic_qtl_pvalue\n')

# Loop through variant gene pairs
for test_name in standard_dynamic_pvalues.keys():
    rs_id = test_name.split('_')[0]
    ensamble_id = test_name.split('_')[1]
    # Extact pvalues for test from two dictionaries
    standard_pvalue = standard_dynamic_pvalues[test_name]
    slope_pvalue = expression_slope_dynamic_qtl_pvalues[test_name]

    t.write(rs_id + '\t' + ensamble_id + '\t' + str(standard_pvalue) + '\t' + str(slope_pvalue) + '\n')
t.close()