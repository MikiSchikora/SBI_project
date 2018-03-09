# this script takes the pdbs in TEMPLATES/ and creates an information dictionary

# We assume that each of the provided PDBs will have always two chains (some of them being the same molecule (a homodimer))


import Bio.PDB as pdb
import os

List_PDBs = os.listdir('./TEMPLATES/')
print('hi')

PDB_info = {} # This will contain information about the unique chains: Keys: filename, Values: {Chain: unique ID}
Unique_Chains = {} # This has Keys: sequence, Value: Integrer, that is the ID of the unique chain

ppb = pdb.Polypeptide.PPBuilder() # polypeptide parser
parser = pdb.PDBParser(PERMISSIVE=1) # pdb parser

ID = 0

for file in List_PDBs:

   path = './TEMPLATES/'+file
   chains = list(parser.get_structure('structure', path).get_chains())
   Chains_dict = {}

   for chain in chains:
      aaSeq = str(ppb.build_peptides(chain)[0].get_sequence())

      if chain.id=='B' or chain.id=='F':
         print(chain.id, aaSeq)

      # add to Unique_Chains if not there
      if aaSeq not in Unique_Chains.keys():
         Unique_Chains[aaSeq] = ID
         ID += 1

      Chains_dict[chain.id] = Unique_Chains[aaSeq]

   PDB_info[file] = Chains_dict








   print(file)