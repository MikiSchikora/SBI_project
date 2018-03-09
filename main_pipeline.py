
# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import Functions_miki as func

# Parse the input arguments:

input = "./TEMPLATES/"
output = './output.pdb'

stoichiometry_file = './stoch.tbl'  # a file containing information about the stoichiometry. This is mandatory
number_subunits_file = './subunits_num.tbl' # a file containing the number of subunits, if known
subunits_seq_file = './subunits_seq.fa'  # a file containing the sequences of the subunits to include. This is mandatory for any chain that is not a random DNA seuqnce.

# process the input

if os.path.isfile(input):

   # when the input is a file you have to generate all the interacting pairs, rotated and translated

   if input.split('.')[-1] != 'pdb':
      raise EnvironmentError('The provided complex has to be a PDB file')

   func.Generate_pairwise_subunits_from_pdb(input)

   Templates_dir = './TEMPLATES/'

elif os.path.isdir(input):

   Templates_dir = input

else:
   raise EnvironmentError('You have to provide a valid input path')


# generate info about the Templates
PDB_info, Uniq_seqs = func.Generate_PDB_info(Templates_dir,subunits_seq_file)





