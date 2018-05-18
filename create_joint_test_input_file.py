import numpy as np 
import os
import sys
import pdb

# Function to extract environmental_variable when environmental_variable_form=='time_steps'
def extract_environmental_variable_time_step_form(sample_name):
    return sample_name.split('_')[1]

def extract_environmental_variable_pseudotime_form(sample_name, mapping_file):
    f = open(mapping_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        line_sample = data[0]
        pseudotime = data[1]
        if line_sample == sample_name:
            return pseudotime
    f.close()
    # ERROR, sample name not in nirmals file
    print('ASSUMPTION ERROROR in nirmals pseudotime file')
    pdb.set_trace()
    return

def extract_max_environmental_variable_pseudotime_form(sample_name, mapping_file):
    max_pseudotime = -1
    f = open(mapping_file)
    for line in f:
        line = line.rstrip()
        data = line.split(',')
        line_sample = data[0]
        pseudotime = int(data[1])
        if sample_name.split('_')[0] == line_sample.split('_')[0]:
            if pseudotime > max_pseudotime:
                max_pseudotime = pseudotime
    f.close()

    return str(max_pseudotime)


def extract_environmental_variable_uniform_form(sample_name, num_states):
    true_time_step = int(sample_name.split('_')[1])
    return str(int(np.floor(true_time_step/num_states)))

def get_cell_line_specific_pc(sample_name, file_name, pc_num):
    cell_liner = sample_name.split('_')[0]
    f = open(file_name)
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        if cell_liner == data[0]:
            return data[pc_num]

input_directory = sys.argv[1]  # Directory containing all CHT input files
output_file = sys.argv[2]  # Output file
environmental_variable_form = sys.argv[3]  # String option describing how to parameterize the environmetal variable
pseudotime_predictions_3_file = sys.argv[4]
pseudotime_predictions_4_file = sys.argv[5]
pseudotime_predictions_5_file = sys.argv[6]
cell_line_specific_pc_file = sys.argv[7]

t = open(output_file, 'w')  # Open output file handle
t.write('sample_id\tenvironmental_variable\tcht_input_file\tcell_line_pc1\tcell_line_pc2\tcell_line_pc3\n')

for file_name in sorted(os.listdir(input_directory)):
    if file_name.startswith('haplotype_read_counts_rand_hap_cis_') == False:
        continue
    if file_name.endswith('.txt.gz') == False:
        continue
    # Extract sample id from filename
    sample_name = file_name.split('.')[2]
    # Extract environmental variable (depends on environmental_variable_form)
    if environmental_variable_form == 'time_steps':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)
    elif environmental_variable_form == 'pseudotime_predictions_3':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_3_file)
    elif environmental_variable_form == 'pseudotime_predictions_4':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_4_file)
    elif environmental_variable_form == 'pseudotime_predictions_5':
        environmental_variable = extract_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_5_file)
    elif environmental_variable_form == 'uniform_4':
        environmental_variable = extract_environmental_variable_uniform_form(sample_name, 4)
    elif environmental_variable_form == 'time_steps_pseudotime':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)
        pseudotime = extract_max_environmental_variable_pseudotime_form(sample_name, pseudotime_predictions_5_file)
        environmental_variable = environmental_variable + '_' + pseudotime
    elif environmental_variable_form == 'time_steps_endpoints':
        environmental_variable = extract_environmental_variable_time_step_form(sample_name)
        if environmental_variable != '0' and environmental_variable != '15':
            continue
    pc1 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 1)
    pc2 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 2)
    pc3 = get_cell_line_specific_pc(sample_name, cell_line_specific_pc_file, 3)
    # Print information to output file
    t.write(sample_name + '\t' + environmental_variable + '\t' + input_directory + file_name + '\t' + pc1 + '\t' + pc2 + '\t' + pc3 + '\n')
t.close()