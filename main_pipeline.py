# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import Functions as func
import Bio.PDB as pdb

# Parse the input arguments:

input_f = './TEMPLATES/'  # "./3kuy.pdb"
output = './output.pdb'

stoichiometry_file = './stoch.tbl'  # a file containing information about the stoichiometry. This is mandatory
number_subunits_file = './subunits_num.tbl'  # a file containing the number of subunits, if known
subunits_seq_file = './subunits_seq.fa'  # a file containing the sequences of the subunits to include. This is mandatory for any chain that is not a random DNA seuqnce.

# process the input

if os.path.isfile(input_f):

    # when the input is a file you have to generate all the interacting pairs, rotated and translated

    if input_f.split('.')[-1] != 'pdb':
        raise EnvironmentError('The provided complex has to be a PDB file')

    Templates_dir = './TEMPLATES/'

    func.Generate_pairwise_subunits_from_pdb(input_f, Templates_dir)

elif os.path.isdir(input_f):

    Templates_dir = input_f

else:
    raise EnvironmentError('You have to provide a valid input path')

# generate info about the Templates
# PDB_info information about the unique chains: Keys: filename, Values: {Chain: unique ID}
# Uniq_seqs is a set with all the unique IDs
PDB_info, Uniq_seqs = func.Generate_PDB_info(Templates_dir, subunits_seq_file)

# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# parse file 1 and 2 and get a structure for each
for filename1 in os.listdir(Templates_dir):

    # we start with the structure of the first pairwise interaction, this is now the current model
    current_structure = p.get_structure("pr1", Templates_dir + filename1)
    for chain in current_structure.get_chains():
        curr_id=chain.id
        chain.id=[x for x in PDB_info[filename1] if x[0]==curr_id][0]


    final_model=func.build_complex(current_structure, Templates_dir, PDB_info)
