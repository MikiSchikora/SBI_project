# this is the structure in created structures

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
import itertools as iter
import rmsd

def get_list_of_common_res(reslist1, reslist2):
    """Takes 2 lists of residues and returns a list of common residues (equal number)"""
    common_res = [x for x in reslist1 if x.get_id()[1] in [y.get_id()[1] for y in reslist2]]
    return common_res


def get_atom_list_from_res_list(reslist):
    """Get a list of atom objects from a list of residue objects"""
    atomlist = []
    for res in reslist:
        for atom in res.get_atoms():
            if atom.id == 'CA':
                atomlist.append(atom)
    return atomlist

def hello():

    res_chain1 = list(curr_struct[0][eq_chain1.id].get_residues())
    res_chain2 = list(struct2[0][eq_chain2.id].get_residues())

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
    rms = sup.rms
    # print(sup.rotran)
    # print(sup.rms)

    # rotate moving atoms
    sup.apply(list(struct2[0][moving_chain.id].get_atoms()))


def structure_in_created_structures(structure, created_structures):

    """ This function asks if the structure (with many chains) is already in  created_structures (a list of structures). returning a boolean if so"""

    # make a deepcopy of these objects
    structure = copy.deepcopy(structure)
    created_structures = copy.deepcopy(created_structures)

    # Get the ids of the chains in structure
    #chain_ids_structure = tuple(sorted([x.id.split('|||')[1] for x in structure.get_chains()]))
    chain_ids_structure = tuple(sorted([x.id for x in structure.get_chains()]))

    # loop through each of the contents of created_structures:
    for created_structure in created_structures:

        # get the chains of created_structure
        #chain_ids_created_structure = tuple(sorted([x.id.split('|||')[1] for x in created_structure.get_chains()]))
        chain_ids_created_structure = tuple(sorted([x.id for x in created_structure.get_chains()]))

        # a boolean that indicates structure matches this created_structures
        structure_is_matching = False

        # ask if the number of each and ids of the chains are the same:
        if chain_ids_structure == chain_ids_created_structure:


            # create a list of the possible partners
            possible_partners = list(chain_ids_created_structure)

            # loop through all the chains in structure
            for chain_str in structure.get_chains():

                # a boolean that indicates if this chain in structure have a partner in created_strcuture
                chain_has_a_partner = False

                #id_str = chain_str.id.split('|||')[1]
                id_str = chain_str.id

                # try to find a partner in created_structure:
                for chain_created_str in created_structure.get_chains():

                    # id_created_str = chain_created_str.id.split('|||')[1]
                    id_created_str = chain_created_str.id

                    # if they have the same id they are potential partners. The id_created_str has also to be avaliable in possible_partners
                    if id_str==id_created_str and id_created_str in possible_partners:

                        # align them:
                        res_chain1 = list(chain_str.get_residues())
                        res_chain2 = list(chain_created_str.get_residues())

                        # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
                        # so first we obtain a list of the common residues
                        common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
                        common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

                        # then we obtain a list of atom objects to use it later
                        common_atoms_s1 = get_atom_list_from_res_list(common_res_s1)
                        common_atoms_s2 = get_atom_list_from_res_list(common_res_s2)

                        rms = 100 # initialize the rms

                        # some molecules don't have CA, like DNA or water molecules
                        if len(common_atoms_s1)>0:

                            # use the Superimposer
                            sup = pdb.Superimposer()

                            # first argument is fixed, second is moving. both are lists of Atom objects
                            sup.set_atoms(common_atoms_s1, common_atoms_s2)
                            rms = sup.rms

                        # if it is a partner
                        if rms <= 2.0:

                            chain_has_a_partner = True

                            # delete this chain_created_str as it already has a partner
                            possible_partners.remove(id_created_str)

                            # exit the partner searching
                            break

                # go to the next  created_structure if a chain has no partner:
                if not chain_has_a_partner:
                    break

                # determine that structure_is_matching if you have an empty possible_partners
                if len(possible_partners)==0:
                    structure_is_matching = True # this is only going to be true at the end if this is the last iteration
                    break

            # return if the structure is matching:
            if structure_is_matching:
                return True

    # if you didn't find any match return false:
    return False

# check:
parser = pdb.PDBParser(PERMISSIVE=1)
structure1 = parser.get_structure('pdb_name', './TEMPLATES/A_and_B.pdb')
cretaed_structures = [copy.deepcopy(structure1)]


a = structure_in_created_structures(structure1,cretaed_structures)
print(a)
