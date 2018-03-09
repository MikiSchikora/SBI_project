
# This script inputs a PDB complex and splits it into all possible interacting pairs of chains.


import Bio.PDB as pdb
from Bio.PDB import Structure as pdb_struct
from Bio.PDB import Model as pdb_model
import numpy as np
import os
import copy

parser = pdb.PDBParser(PERMISSIVE=1)

structure = parser.get_structure('3kuy', '3kuy.pdb')

model = structure[0]

# free the ./TEMPLATES_DIR/
os.system('rm -rf ./TEMPLATES/*')

# initialize the saved pairs
saved_pairs = set()

for chain1 in model.get_chains():

   #new_structure = chain1

   for chain2 in model.get_chains():

      #comb = chain1.id+'.'+chain2.id
      #comb_rev = chain2.id+'.'+chain1.id

      if chain1 is not chain2:

         # save the combination
         #saved_pairs.add(comb)
         #saved_pairs.add(comb_rev)

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
                  if distance < 15:
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

            rotation = np.random.rand(3,3)
            translation = np.array((0, 0, 1), 'f')

            for atom in new_structure.get_atoms():
               #
               #
               #atom.transform(rotation, translation)
               pass

            # write to new pdb:

            io = pdb.PDBIO()
            io.set_structure(new_structure)
            io.save('TEMPLATES/'+ID+'.pdb')








         #new_structure.add(chain2.get_residues())



