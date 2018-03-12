
# This is a python module that includes the functions and classes used in the main_pipeline.py


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

               rotation = np.array([[0.01,0.01,0.03],[0.01,0.03,0.02],[0.05,0.01,0.02]])
               translation = np.array((0, 0, 1), 'f')

               for atom in new_structure.get_atoms():
                  atom.transform(rotation, translation)

               # write to new pdb:

               io = pdb.PDBIO()
               io.set_structure(new_structure)
               io.save(TEMPLATES_path+ID+'.pdb')

def Generate_PDB_info(Templates_dir,subunits_seq_file, min_identity_between_chains = 30):

   """This function takes the Templates_dir and creates a dictionary with information about each template

   Each template get's incorporated into a dictionary with a unique identifiers for it's chains, in a way that you save all files that are useful

   min_identity_between_chains refers to the minimum identity required between two molecules to state that they are the same

   subunits_seq_file is a fasta file containing all the molecules that form the complex"""

   List_PDBs = os.listdir(Templates_dir)

   PDB_info = {} # This will contain information about the unique chains: Keys: filename, Values: {Chain: unique ID}

   # create a dictionary that has the known sequences
   Seqs ={}
   for record in Seq_IO.parse(subunits_seq_file, "fasta"):
      Seqs[record.id] = record.seq._data


   ppb = pdb.Polypeptide.PPBuilder() # polypeptide parser
   parser = pdb.PDBParser(PERMISSIVE=1) # pdb parser

   for file in List_PDBs:

      path = Templates_dir+file
      chains = list(parser.get_structure('structure', path).get_chains())
      Chains_dict = {} # wll contain the unique ID for each chain

      for chain in chains:
         aaSeq = "".join([str(x.get_sequence()) for x in ppb.build_peptides(chain)])

         # Assign an ID (id of the fasta or a new for the molecule) to the sequence:

         Seq_id = None
         Best_identity = min_identity_between_chains

         for my_id in Seqs.keys():

            alignment = pairwise2.align.globalxx(aaSeq, Seqs[my_id])
            percent_match = (alignment[0][2] / len(min(aaSeq, Seqs[my_id]))) * 100

            if percent_match > Best_identity:
               Seq_id = my_id
               Best_identity = percent_match

         # if it doesn't fit anything, it may be DNA or RNA, noty included in the provided fasta file

         if Seq_id is None and all(i in 'ACUTG' for i in aaSeq):
            Seq_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6)) #random ID
            Seqs[Seq_id] = aaSeq

         if Seq_id is None:
            print(file+" contains non provided and non acid nucleic sequences!!")
            break

         else:
            Chains_dict[chain.id] = Seq_id

      PDB_info[file] = Chains_dict

   return PDB_info, set(Seqs.keys())


