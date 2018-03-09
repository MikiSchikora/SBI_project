
# This script inputs a PDB complex and splits it into all possible interacting pairs of chains.


import Bio.PDB as pdb
from Bio.PDB import Structure as pdb_struct
from Bio.PDB import Model as pdb_model
import os

parser = pdb.PDBParser(PERMISSIVE=1)

structure = parser.get_structure('3kuy', '3kuy.pdb')

model = structure[0]

# free the ./TEMPLATES_DIR/
os.system('rm -rf ./TEMPLATES/*')

for chain1 in model.get_chains():

   #new_structure = chain1

   for chain2 in model.get_chains():

      if chain1 is not chain2:

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
                  if distance < 6:
                     chains_interacting = 1

         if chains_interacting==1:

            # create a structure object
            ID = chain1.id+'_and_'+chain2.id
            new_structure = pdb_struct.Structure(ID)
            new_model = pdb_model.Model(0)
            new_model.add(chain1)
            new_model.add(chain2)
            new_structure.add(new_model)

            # move the coordinates of the structure to simulate what would happen if they were coming from different files

            # write to new pdb:

            io = pdb.PDBIO()
            io.set_structure(new_structure)
            io.save('TEMPLATES/'+ID+'.pdb')








         #new_structure.add(chain2.get_residues())



