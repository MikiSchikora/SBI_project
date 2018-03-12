
# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import Functions_miki as func
import Bio.PDB as pdb

# Parse the input arguments:

input_f = "./3kuy.pdb"
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

m
# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# parse file 1 and 2 and get a structure for each
filename1 = "TEMPLATES/A_and_B.pdb"
structure1 = p.get_structure("pr1", filename1)

common_chain='A'
for filename2 in os.listdir('TEMPLATES'):
   if filename2.startswith(common_chain) and filename2[6] != 'B':
      structure2 = p.get_structure("pr2", "TEMPLATES/"+filename2)

      rotating_chain = filename2[6]

      # same chain is retrieved from the 2 structures. Example: chain A
      common_chain_s1 = structure1[0][common_chain]
      common_chain_s2 = structure2[0][common_chain]

      # get the atoms of the common chain in a list
      common_chain_atoms_s1 = list(common_chain_s1.get_atoms())
      common_chain_atoms_s2 = list(common_chain_s2.get_atoms())

      # use the Superimposer
      sup = pdb.Superimposer()

      # first argument is fixed, second is moving. both are lists of Atom objects
      sup.set_atoms(common_chain_atoms_s1, common_chain_atoms_s2)
      print(sup.rotran)
      print(sup.rms)

      # rotate moving atoms
      sup.apply(list(structure2[0][rotating_chain].get_atoms()))

      # add to the fixed structure, the moved chain
      structure1[0].add(structure2[0][rotating_chain])

      # save in a pdb file
      io = pdb.PDBIO()
      io.set_structure(structure1)
      io.save('out1.pdb')