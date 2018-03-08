
# This script inputs a PDB complex and splits it into all possible interacting pairs of chains.


from Bio.PDB import *

parser = PDBParser(PERMISSIVE=1)

structure = parser.get_structure('3kuy', '3kuy.pdb')

model = structure[0]

for chain1 in model.get_chains():

   new_structure = chain1

   for chain2 in model.get_chains():

      if chain1 is not chain2:

         new_structure.add(chain2.get_residues())



