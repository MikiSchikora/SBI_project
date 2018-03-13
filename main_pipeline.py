
# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import Functions_miki as func
import Bio.PDB as pdb

# Parse the input arguments:

input_f = './TEMPLATES/' #"./3kuy.pdb"
output = './output.pdb'

stoichiometry_file = './stoch.tbl'  # a file containing information about the stoichiometry. This is mandatory
number_subunits_file = './subunits_num.tbl' # a file containing the number of subunits, if known
subunits_seq_file = './subunits_seq.fa'  # a file containing the sequences of the subunits to include. This is mandatory for any chain that is not a random DNA seuqnce.

# process the input

if os.path.isfile(input_f):

   # when the input is a file you have to generate all the interacting pairs, rotated and translated

   if input_f.split('.')[-1] != 'pdb':
      raise EnvironmentError('The provided complex has to be a PDB file')

   Templates_dir = './TEMPLATES/'

   func.Generate_pairwise_subunits_from_pdb(input_f,Templates_dir)

elif os.path.isdir(input_f):

   Templates_dir = input_f

else:
   raise EnvironmentError('You have to provide a valid input path')


# generate info about the Templates
# PDB_info information about the unique chains: Keys: filename, Values: {Chain: unique ID}
# Uniq_seqs is a set with all the unique IDs
PDB_info, Uniq_seqs = func.Generate_PDB_info(Templates_dir,subunits_seq_file)


# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# parse file 1 and 2 and get a structure for each
for filename1 in os.listdir(Templates_dir):
   if filename1.startswith("A"):
      # we start with the structure of the first pairwise interaction, this is now the current model
      current_structure = p.get_structure("pr1", Templates_dir+filename1)
      current_chains=PDB_info[filename1]

      # iterate through the chains of the current model
      for common_chain1, common_id in current_chains.items():

         # check if another file has this chain
         for filename2 in os.listdir(Templates_dir):
            # this file needs to be different from the first and also have a subunit in common with it
            rotating_chain=None
            common_chain2=None
            structure2=None
            if filename2 != filename1 and common_id in PDB_info[filename2].values():

               structure2 = p.get_structure("pr2", Templates_dir+filename2)

               # label the two chains as common (if equal id as previous chain) or rotating
               for other_chain, other_id in PDB_info[filename2].items():
                  if other_id != common_id:
                     rotating_chain = other_chain
                  else:
                     common_chain2 = other_chain

               # if this is an homodimer, the first chain is set to be the rotating
               if list(PDB_info[filename2].values())[0] == list(PDB_info[filename2].values())[1]:
                  rotating_chain=list(PDB_info[filename2])[0]

               current_structure=func.superimpose_and_rotate(common_chain1, common_chain2, rotating_chain, current_structure, structure2)

               # save in a pdb file
               io = pdb.PDBIO()
               io.set_structure(current_structure)
               io.save('out1.pdb')
