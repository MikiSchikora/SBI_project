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

               #rotation = np.array([[0.01,0.01,0.03],[0.01,0.03,0.02],[0.05,0.01,0.02]])
               rotation = np.array([[1,0,0],[0,1,0],[0,0,1]])
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

   subunits_seq_file is a fasta file containing all the molecules that form the complex

   It also writes a new set of PDBs in ./NEW_TEMPLATES so that the chains are organized accordingly"""

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
      PDB_structure = parser.get_structure('structure', path)
      chains = list(PDB_structure.get_chains())
      Chain_IDs = [] # wll contain the unique ID for each chain

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
            ID_new = str(chain.id+"|||"+Seq_id)
            Chain_IDs.append(ID_new)

      PDB_info[file] = Chain_IDs


   return PDB_info, set(Seqs.keys())


def is_Steric_clash(structure, rotating_chain):

   """This function returns a boolean that indicates if a rotating_chain clashes against
    a structure.

    The clash crteria is that at least one of the atoms are at a lower distance than
    min_VDW_distance"""

   Van_der_Waals_radius = {'H':1.2,'C':1.7,'N':1.55,'O':1.52,'F':1.47,'P':1.8,'S':1.8}

   NS = pdb.NeighborSearch(list(structure.get_atoms()))

   Number_clashes = 0

   for at in rotating_chain.get_atoms():
      neighbors = NS.search(at.get_coord(), 1.5)



      # if you find a close atom, ask if it is a clash
      if len(neighbors)>0:

         for neigh in neighbors:

            sum_VDW_radius = Van_der_Waals_radius[neigh.element] + Van_der_Waals_radius[at.element]
            Distance = neigh - at

            if Distance < sum_VDW_radius:
               Number_clashes += 1

            if Number_clashes==1:
               return True

   return False


def get_all_res_list(struct, chain):
    """Get a list with all residues of a particular chain"""
    print(struct, chain)
    curr_chain = struct[0][chain.id]
    chain_res = [x for x in curr_chain.get_residues()]
    return chain_res


def get_list_of_common_res(reslist1, reslist2):
    """Takes 2 lists of residues and returns a list of common residues (equal number)"""
    common_res = [x for x in reslist1 if x.get_id()[1] in [y.get_id()[1] for y in reslist2]]
    return common_res


def get_atom_list_from_res_list(reslist):
    """Get a list of atom objects from a list of residue objects"""
    atomlist = []
    for res in reslist:
        for atom in res.get_atoms():
            atomlist.append(atom)
    return atomlist


def superimpose_and_rotate(eq_chain1, eq_chain2, moving_chain, curr_struct, struct2):
    # all residues from same chain (common chain) are retrieved from the 2 structures. Example: chain A
    # print(eq_chain1, eq_chain2, moving_chain, curr_struct, struct2)
    #res_chain1 = get_all_res_list(curr_struct, eq_chain1)
    #res_chain2 = get_all_res_list(struct2, eq_chain2)

    #print(eq_chain1.id)

    res_chain1 = list(curr_struct[0][eq_chain1.id].get_residues())
    res_chain2 = list(struct2[0][eq_chain2.id].get_residues())
    #print('residues 'res_chain1)
    #print(res_chain2)

    # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
    # so first we obtain a list of the common residues
    common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
    common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

    # then we obtain a list of atom objects to use it later
    common_atoms_s1 = get_atom_list_from_res_list(common_res_s1)
    common_atoms_s2 = get_atom_list_from_res_list(common_res_s2)

    # use the Superimposer
    sup = pdb.Superimposer()

    # first argument is fixed, second is moving. both are lists of Atom objects
    sup.set_atoms(common_atoms_s1, common_atoms_s2)
    # print(sup.rotran)
    # print(sup.rms)

    # rotate moving atoms
    sup.apply(list(struct2[0][moving_chain.id].get_atoms()))

    # add to the fixed structure, the moved chain
    added=0
    if not is_Steric_clash(curr_struct,moving_chain):

        my_id=moving_chain.id
        chain_names = [x.id for x in curr_struct[0].get_chains()]
        #print('chain_names',chain_names)
        while(added==0):
            rand = '|||'+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))  # random ID
            if my_id+rand not in chain_names:
                moving_chain.id=my_id+rand

                print('id if the moving chain',moving_chain.id)

                curr_struct[0].add(struct2[0][moving_chain.id])
                added=1
    else:
        print("problem:clash")

    return curr_struct, added



def build_complex(current_str, mydir, PDB_dict):
    """This function builds a complex from a set of templates"""
    sth_added=0
    p = pdb.PDBParser(PERMISSIVE=1)
    for chain1 in current_str[0].get_chains():
        id_chain1 = chain1.id.split('|||')[1]

        for filename2 in os.listdir(mydir):
            rotating_chain = None
            common_chain2 = None
            if id_chain1 in [x.split('|||')[1] for x in PDB_dict[filename2]]:
                structure2 = p.get_structure("pr2", mydir + filename2)
                for chain in structure2.get_chains():
                    curr_id=chain.id
                    chain.id=[x for x in PDB_dict[filename2] if x[0]==curr_id][0]

                #print(PDB_dict)

                #print('id of the first chain', PDB_dict[filename2][0].split('|||')[1])
                #print('id of the second chain ', PDB_dict[filename2][1].split('|||')[1])

                for chain2 in structure2[0].get_chains():
                    id_chain2= chain2.id.split('|||')[1]

                    if PDB_dict[filename2][0].split('|||')[1]==PDB_dict[filename2][1].split('|||')[1]:
                        rotating_chain= structure2[0][PDB_dict[filename2][0]]
                        common_chain2= structure2[0][PDB_dict[filename2][1]]

                    elif id_chain1 == id_chain2:
                        common_chain2 = chain2
                    else:
                        rotating_chain = chain2

                #print('rotating_chain:',rotating_chain,'common_chain:',common_chain2)

                current_str, sth_added = superimpose_and_rotate(chain1, common_chain2, rotating_chain, current_str, structure2)

    if sth_added==1:
        build_complex(current_str,dir,PDB_dict)
    else:
        print(list(current_str.get_chains()))
        #return current_str


        return

        io = pdb.PDBIO()

        hey=current_str[0].get_list()
        print(hey)

        io.set_structure(current_str)
        io.save('out1.pdb')



