import Bio.PDB as pdb
from Bio.PDB import Structure as pdb_struct
from Bio.PDB import Model as pdb_model
import numpy as np
import os
import copy
import Bio.SeqIO as Seq_IO
from Bio import pairwise2
import random
import string

def Generate_pairwise_subunits_from_pdb(pdb_file_path,Templates_dir):

   """This function takes an existing complex and fragments it into each of the pairwise interactions between subunits.

   pdb_file_path is the path where the complex PDB is

   It does not consider nucleic acid sequences, as this is only for testing the program on different complexes"""


   TEMPLATES_path = Templates_dir

   parser = pdb.PDBParser(PERMISSIVE=1)

   structure = parser.get_structure('pdb_name', pdb_file_path)

   model = structure[0]

   # free the ./TEMPLATES_path/
   os.system('rm -rf '+TEMPLATES_path+'*')

   # initialize the saved pairs
   saved_pairs = set()

   for chain1 in model.get_chains():

      for chain2 in model.get_chains():

         comb = chain1.id+'.'+chain2.id
         comb_rev = chain2.id+'.'+chain1.id

         if chain1 is not chain2 and comb not in saved_pairs:

            # save the combination
            saved_pairs.add(comb)
            saved_pairs.add(comb_rev)

            # ask if any of the residues is interacting, if so save the PDB

            chains_interacting = 0

            for residue1 in chain1:
               if chains_interacting==1:
                  break
               for residue2 in chain2:
                  if residue1 != residue2:
                     # compute distance between CA atoms
                     try:
                        distance = residue1['CA'] - residue2['CA']
                     except KeyError:
                        ## no CA atom, e.g. for H_NAG
                        continue
                     if distance < 10:
                        chains_interacting = 1

            if chains_interacting==1:

               # create a structure object
               ID = chain1.id+'_and_'+chain2.id
               new_structure = pdb_struct.Structure(ID)

               new_model = pdb_model.Model(0)
               new_model.add(copy.deepcopy(chain1))
               new_model.add(copy.deepcopy(chain2))
               new_structure.add(new_model)

               # move the coordinates of the structure to simulate what would happen if they were coming from different files

               #rotation = np.array([[0.01,0.01,0.03],[0.01,0.03,0.02],[0.05,0.01,0.02]])
               rotation = np.array([[1,0,0],[0,1,0],[0,0,1]])
               translation = np.array((0, 0, 1), 'f')

               for atom in new_structure.get_atoms():
                  atom.transform(rotation, translation)

               # write to new pdb:

               io = pdb.PDBIO()
               io.set_structure(new_structure)
               io.save(TEMPLATES_path+ID+'.pdb')

input_f = './3kuy.pdb'
Generate_pairwise_subunits_from_pdb(input_f,'./TEMPLATES/')