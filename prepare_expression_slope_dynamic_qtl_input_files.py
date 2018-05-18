import numpy as np 
import os
import sys
import pdb
import gzip




# Make mapping from gene identifier to ensamble id
def make_gene_mapping(target_region_input_file):
    # Get input for only 1 sample 
    # Only need one because genes are the same in all sample input files
    f = open(target_region_input_file)
    gene_mapping = {}  # Initialize gene mapping
    for line in f:
        line = line.rstrip()
        data = line.split()
        ensamble_id = data[6].split('_')[-1]
        gene_identifier = data[0] + '_' + data[7] + '_' + data[8]
        gene_mapping[gene_identifier] = ensamble_id
    f.close()
    return gene_mapping


# Extract dictionary of test_names (variant_gene pairs), variants, and gene_names
def extract_test_names(joint_test_input_file, gene_mapping):
    # Initialize output variables
    test_names = {}
    gene_names = {}
    variant_names = {}
    # First get input file from 1 of the samples (doesn't matter which sample)
    f = open(joint_test_input_file)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        head_count = head_count + 1
        if head_count != 2:  # only need input file for 1 sample (so we take the first sample)
            continue
        # this is the input file
        input_file_name = data[2]
    f.close()
    # Open input file
    g = gzip.open(input_file_name)
    # Stream input file
    head_count = 0  # to skip header
    for line in g:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue

        gene_identifier = data[0] + '_' + data[7] + '_' + data[8]
        # convert from gene identifier to ensamble id using pre-built dictionary
        ensamble_id = gene_mapping[gene_identifier]
        # extract rs_id
        rs_id = data[2]
        # add test to dictionaries
        test_names[rs_id + '_' + ensamble_id] = 0
        gene_names[ensamble_id] = 0
        variant_names[rs_id] = 0
    return test_names, variant_names, gene_names


# Create vectors of names of samples, and their corresponding time step
def get_sample_info(joint_test_input_file):
    # initialize output vectors
    sample_names = []
    time_steps = []
    pc1 = []
    pc2 = []
    pc3 = []
    # Stream joint test input file
    f = open(joint_test_input_file)
    head_count = 0  # skip header
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count = head_count + 1
            continue
        sample_name = data[0]
        time_step = data[1]
        sample_names.append(sample_name)
        time_steps.append(time_step)
        pc1.append(data[3])
        pc2.append(data[4])
        pc3.append(data[5])
    return np.asarray(sample_names), np.asarray(time_steps), np.asarray(pc1), np.asarray(pc2), np.asarray(pc3)

def check_to_make_sure_no_missing_entries(dicti):
    for key in dicti.keys():
        if len(dicti[key]) == 1:
            print('MISSING KEY ERRORR!!')
            pdb.set_trace()
    return

# Now extract gene expression vector for all genes
def extract_gene_expresssion(genes, total_expression_file, sample_names):
    head_count = 0  # for header
    f = open(total_expression_file)
    for line in f:
        line = line.rstrip()
        data = np.asarray(line.split())
        if head_count == 0:  # header
            # Using header, find ordered indices that correspond to order of sample_names
            head_count = head_count + 1
            reordered_indices = []
            for sample_name in sample_names:
                for i, ele in enumerate(data):
                    if ele == sample_name:
                        reordered_indices.append(i)
            if np.array_equal(data[reordered_indices],sample_names) == False:
                print('ASSUMPTION ERRORRO!')
                pdb.set_trace()
            continue
        # Name of gene for current line
        ensamble_id = data[0]
        # Extract gene expression measurements
        ordered_data = data[reordered_indices].astype(float)
        # Standardize gene expression measurements
        standardized_ordered_data = (ordered_data - np.mean(ordered_data))/(np.std(ordered_data))
        # Check if gene is in our dictionary of genes (used genes)
        if ensamble_id in genes:
            # If it is, add gene expression vector
            genes[ensamble_id] = standardized_ordered_data
    f.close()
    check_to_make_sure_no_missing_entries(genes)
    return genes

# Convert vector of sample names (cellLine_timeStep) to a vector of cell lines
def sample_to_cell_line_names(sample_names):
    cell_lines = []
    for ele in sample_names:
        cell_lines.append(ele.split('_')[0])
    return np.asarray(cell_lines)


def extract_genotype_data(variants, genotype_file, sample_names, genotype_version):
    cell_line_names = sample_to_cell_line_names(sample_names)

    f = open(genotype_file)
    for line in f:
        line = line.rstrip()
        data = np.asarray(line.split())
        if line.startswith('#CHROM'):  # header
            # Using header, find ordered indices that correspond to order of sample_names
            reordered_indices = []
            for cell_line in cell_line_names:
                for i, ele in enumerate(data):
                    if ele == cell_line:
                        reordered_indices.append(i)
            if np.array_equal(data[reordered_indices],cell_line_names) == False:
                print('ASSUMPTION ERRORRO!')
                pdb.set_trace()
            continue
        if line.startswith('#'):
            continue
        # Name of variant for current line
        rs_id = data[2]
        # Get genotype data in same order as sample_names
        ordered_data = data[reordered_indices].astype(float)
        # If we want to round genotypes
        if genotype_version == 'round':
            ordered_data = np.round(ordered_data)
        
        # Check if variant is in our dictionary of variants (used variants)
        if rs_id in variants:
            # If it is, add gene expression vector
            variants[rs_id] = ordered_data
    f.close()
    check_to_make_sure_no_missing_entries(variants)
    return variants

def time_steps_extraction(arr):
    time_steps = []
    pseudotime = []
    for ele in arr:
        aaa = ele.split('_')[0]
        bbb = ele.split('_')[1]
        time_steps.append(aaa)
        pseudotime.append(bbb)
    return np.asarray(time_steps), np.asarray(pseudotime)

def get_cell_line_data_structure(unique_lines, sample_names):
    data_struct = {}
    for cell_line in unique_lines:
        data_struct[cell_line] = {}
        for time_step in range(16):
            sample_id = cell_line + '_' + str(time_step)
            if sample_id in sample_names:
                positions = np.where(sample_names==sample_id)[0]
                if len(positions) != 1:
                    print("ASSUMPTIONER OEROROE")
                    pdb.set_trace()
                position = positions[0]
                data_struct[cell_line][time_step] = position
            else:
                data_struct[cell_line][time_step] = -10  # encodes missing values
    return data_struct

def extract_variable_in_cell_line_structured_string(cell_line_data_structure, unique_lines, variable):
    stringers = []
    for cell_line in unique_lines:
        cell_line_arr = []
        for time_step in range(16):
            position = cell_line_data_structure[cell_line][time_step]
            if position != -10:
                cell_line_arr.append(variable[position])
        cell_line_string = ','.join(np.asarray(cell_line_arr).astype(str))
        stringers.append(cell_line_string)
    return ';'.join(np.asarray(stringers))

def extract_variable_in_cell_line_structured_string_for_first_time_step(cell_line_data_structure, unique_lines, variable):
    stringers = []
    for cell_line in unique_lines:
        time_step = 0
        position = cell_line_data_structure[cell_line][time_step]
        if position != -10:
            stringers.append(str(variable[position]))
    return ';'.join(np.asarray(stringers))


def print_helper(output_file, test_names, variants, genes, sample_names, time_steps, environmental_variable_form, pc1, pc2, pc3):
    unique_lines = np.unique(sample_to_cell_line_names(sample_names))

    cell_line_data_structure = get_cell_line_data_structure(unique_lines, sample_names)

    time_step_string = extract_variable_in_cell_line_structured_string(cell_line_data_structure, unique_lines, time_steps)

    pc1_string = extract_variable_in_cell_line_structured_string_for_first_time_step(cell_line_data_structure, unique_lines, pc1)
    pc2_string = extract_variable_in_cell_line_structured_string_for_first_time_step(cell_line_data_structure, unique_lines, pc2)
    pc3_string = extract_variable_in_cell_line_structured_string_for_first_time_step(cell_line_data_structure, unique_lines, pc3)


    # Open file handle to output file
    t = open(output_file, 'w')
    for test_name in sorted(test_names.keys()):
        # Extract rs_id and ensamble_id from test_name
        rs_id = test_name.split('_')[0]
        ensamble_id = test_name.split('_')[1]
        # extract genotype vector
        geno_string = extract_variable_in_cell_line_structured_string_for_first_time_step(cell_line_data_structure, unique_lines, variants[rs_id])
        # extract total expression vector
        te_string = extract_variable_in_cell_line_structured_string(cell_line_data_structure, unique_lines, genes[ensamble_id])
        # Print
        t.write(rs_id + '\t' + ensamble_id + '\t' + ';'.join(unique_lines) + '\t' + pc1_string + '\t' + pc2_string + '\t' + pc3_string + '\t' + geno_string + '\t' + time_step_string + '\t' + te_string + '\n')
    t.close()

##########################
# Command line args
##########################
# File containing all samples to be tested for dynamic qtl analysis (along with their time step)
joint_test_input_file = sys.argv[1]
# File containing info on all target regions (ie the names of the genes and their corresponding genomic locations)
target_region_input_file = sys.argv[2]
# Total expression file
total_expression_file = sys.argv[3]
# Genotype file
genotype_file = sys.argv[4]
# Whether to use rounded or imputed genotypes
genotype_version = sys.argv[5]
# Output file
output_file = sys.argv[6]
# Form of environmental variables
environmental_variable_form = sys.argv[7]



## Create mapping from gene location to gene name
gene_mapping = make_gene_mapping(target_region_input_file)


# Create vectors of names of samples, and their corresponding time step
sample_names, time_steps, pc1, pc2, pc3 = get_sample_info(joint_test_input_file)


# Extract dictionary of test_names (variant_gene pairs), variants, and gene_names
test_names, variants, genes = extract_test_names(joint_test_input_file, gene_mapping)


# Now extract gene expression vector for all genes
genes = extract_gene_expresssion(genes, total_expression_file, sample_names)



# Now extract genotype vector for all variants
variants = extract_genotype_data(variants, genotype_file, sample_names, genotype_version)
print('start')

print_helper(output_file, test_names, variants, genes, sample_names, time_steps, environmental_variable_form, pc1, pc2, pc3)