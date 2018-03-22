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
import rmsd


def Generate_pairwise_subunits_from_pdb(pdb_file_path, TEMPLATES_path, file_type):
    """This function takes an existing complex and fragments it into each of the pairwise interactions between subunits.
    pdb_file_path is the path where the complex PDB is
    It does not consider nucleic acid sequences, as this is only for testing the program on different complexes"""

    num_file = 0

    if file_type == 'PDB':
        parser = pdb.PDBParser(PERMISSIVE=1)
    else:
        parser = pdb.MMCIFParser()

    structure = parser.get_structure('pdb_name', pdb_file_path)

    model = structure[0]

    # free the ./TEMPLATES_path/
    os.system('rm -rf ' + TEMPLATES_path + '*')

    # initialize the saved pairs
    saved_pairs = set()

    for chain1 in model.get_chains():

        for chain2 in model.get_chains():

            comb = chain1.id + '.' + chain2.id
            comb_rev = chain2.id + '.' + chain1.id

            if chain1 is not chain2 and comb not in saved_pairs:

                # save the combination
                saved_pairs.add(comb)
                saved_pairs.add(comb_rev)

                # ask if any of the residues is interacting, if so save the PDB

                chains_interacting = 0

                for residue1 in chain1:
                    if chains_interacting == 1:
                        break
                    for residue2 in chain2:
                        if residue1 != residue2:
                            # compute distance between CA atoms
                            try:
                                distance = residue1['CA'] - residue2['CA']
                            except KeyError:
                                # no CA atom, e.g. for H_NAG
                                continue
                            if distance < 5:
                                chains_interacting = 1

                if chains_interacting == 1:

                    # create a structure object
                    ID = str(num_file)
                    num_file += 1
                    new_structure = pdb_struct.Structure(ID)

                    new_model = pdb_model.Model(0)
                    new_model.add(copy.deepcopy(chain1))
                    new_model.add(copy.deepcopy(chain2))

                    # check that the chain.id is correct (only 1 char)
                    # if it has 2 chars, change them to lowercase ('a' and 'b')
                    for new_chain in new_model.get_chains():
                        if len(new_chain.id) == 2:
                            if "a" in new_model:
                                new_chain.id = "b"
                            else:
                                new_chain.id = "a"

                    new_structure.add(new_model)

                    # move the coordinates of the structure to simulate what would happen if they were coming from different files

                    # rotation = np.array([[0.01,0.01,0.03],[0.01,0.03,0.02],[0.05,0.01,0.02]])
                    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                    translation = np.array((0, 0, 1), 'f')

                    for atom in new_structure.get_atoms():
                        atom.transform(rotation, translation)

                    # write to new pdb:

                    io = pdb.PDBIO()
                    io.set_structure(new_structure)
                    io.save(TEMPLATES_path + ID + '.pdb')


def create_random_chars_id(n):
    """Creates an identifier with n random characters."""
    rand_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))
    return rand_id


def add_new_seqid_to_dicts(seq, dict1, dict2, n):
    """Adds a random identifier of n characters to 2 dictionaries, the first having a sequence, the other initialising and empty set"""
    seq_id = create_random_chars_id(n)  # create random ID
    while seq_id in dict1:
        seq_id = create_random_chars_id(n)  # random ID is already in Seqs dict
    dict1[seq_id] = seq
    dict2[seq_id] = set()

    return seq_id, dict1[seq_id], dict2[seq_id]


def Generate_PDB_info(Templates_dir, subunits_seq_file, min_identity_between_chains=95):
    """This function takes the Templates_dir and creates a dictionary with information about each template
    Each template gets incorporated into a dictionary with a unique identifiers for it's chains, in a way that you save all files that are useful
    min_identity_between_chains refers to the minimum identity required between two molecules to state that they are the same
    subunits_seq_file is a fasta file containing all the molecules that form the complex
    it also returns a dictionary with keys: chain ids, values: array of filenames that contain that key"""

    List_PDBs = os.listdir(Templates_dir)
    PDB_info = {}  # This will contain information about the unique chains: Keys: filename, Values: {Chain: unique ID}
    Seq_to_filenames = {}  # dictionary that contains the filenames that contain each seq ID

    ppb = pdb.Polypeptide.PPBuilder()  # polypeptide parser
    parser = pdb.PDBParser(PERMISSIVE=1)  # pdb parser

    # create a dictionary that has the known sequences if a multifasta is given
    Seqs = {}
    if subunits_seq_file:
        for record in Seq_IO.parse(subunits_seq_file, "fasta"):
            Seqs[record.id] = record.seq._data
            Seq_to_filenames[record.id] = set()

    for file in List_PDBs:
        path = Templates_dir + file
        PDB_structure = parser.get_structure('structure', path)
        chains = list(PDB_structure.get_chains())
        Chain_IDs = []  # wll contain the unique ID for each chain

        for chain in chains:
            # obtain the sequence of the chain
            aaSeq = "".join([str(x.get_sequence()) for x in ppb.build_peptides(chain)])
            if aaSeq == '':
                aaSeq = create_random_chars_id(len(chain))

            # assign an ID (id of the fasta or a new for the molecule) to the sequence:
            Seq_id = None
            Best_identity = min_identity_between_chains

            # compare the current chain sequence with the ones present in Seqs
            if Seqs:
                for my_id in Seqs:
                    # align my current sequence with the sequences stored in Seqs dictionary
                    alignment = pairwise2.align.globalxx(aaSeq, Seqs[my_id])
                    # compute percentage identity
                    percent_match = (alignment[0][2] / min(len(aaSeq), len(Seqs[my_id]))) * 100

                    # if percentage identity is lar
                    if percent_match > Best_identity:
                        Seq_id = my_id
                        Best_identity = percent_match

                # my current sequence didn't match anything in Seqs, so it is added to the dictionary
                if Seq_id is None and subunits_seq_file is None:
                    Seq_id, Seqs[Seq_id], Seq_to_filenames[Seq_id] = add_new_seqid_to_dicts(aaSeq, Seqs, Seq_to_filenames, 6)

            # add the current chain sequence to Seqs if it is empty
            else:
                Seq_id, Seqs[Seq_id], Seq_to_filenames[Seq_id] = add_new_seqid_to_dicts(aaSeq, Seqs, Seq_to_filenames, 6)

            # if my current seq doesn't fit anything, it may be DNA or RNA, not included in the provided fasta file
            if Seq_id is None and all(i in 'ACUTG' for i in aaSeq): # only if a multifasta file is provided
                Seq_id, Seqs[Seq_id], Seq_to_filenames[Seq_id] = add_new_seqid_to_dicts(aaSeq, Seqs, Seq_to_filenames, 6)

            if Seq_id is None:  # only if a multifasta file is provided
                print(file + " contains non provided and non acid nucleic sequences!!")
                break
            # save info info dictionary Seq_to_filenames {Seq_id:filename}
            else:
                ID_new = str(chain.id + "|||" + Seq_id)
                Chain_IDs.append(ID_new)
                Seq_to_filenames[Seq_id].add(file)

        # save info into dictionary PDB_info {filename:chain_ids_list}
        if len(Chain_IDs) == 2:
            PDB_info[file] = Chain_IDs

    return PDB_info, Seq_to_filenames


def is_Steric_clash(structure, rotating_chain, distance_for_clash=1):

    """This function returns False if there's no clash between a rotating chain and the current structure.
    If there's a clash it can return 1 (clash between two different chains) or 2 (you are trying to superimpose something in the place it was already)

    it also returns the ids of the chains in structure that are clashing

    The clash criteria is that at least one of the atoms are at a lower distance than distance_for_clash A"""

    NS = pdb.NeighborSearch(list(structure.get_atoms()))

    clashing_chains = set()
    n_clashes = 0

    for at in rotating_chain.get_atoms():
        neighbors = NS.search(at.get_coord(), distance_for_clash)

        if len(neighbors) > 0:
            for neigh in neighbors:
                clashing_chains.add(neigh.get_parent().get_parent().id)
                n_clashes += 1

    if len(clashing_chains) > 1:
        # a clash against different chains:
        val_to_return = 1

    elif len(clashing_chains) == 1 and n_clashes > 100 and rotating_chain.id.split('|||')[1] == list(clashing_chains)[0].split('|||')[1]:

        # a clash MAYBE because you are trying to superimpose something in the place it was already
        # define the clashing chain:
        clash_chain = structure[0][list(clashing_chains)[0]]
        res_chain1 = list(clash_chain.get_residues())
        res_chain2 = list(rotating_chain.get_residues())

        # so first we obtain a list of the common residues
        common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
        common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

        # then we obtain a list of atom objects to use it later
        common_atoms_s1 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s1) if x.id == 'CA'])
        common_atoms_s2 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s2) if x.id == 'CA'])

        RMSD = rmsd.kabsch_rmsd(common_atoms_s1, common_atoms_s2)
        # print("clash chain: ", clash_chain, " rotating chain: ", rotating_chain, RMSD)

        if RMSD <= 2.0:
            # it is the same chain
            val_to_return = 2

        else:
            # it is another chain or the same with different structure
            val_to_return = 1

    elif n_clashes > 0:
        # it is another chain
        val_to_return = 1

    else:
        # no clash
        val_to_return = False

    if val_to_return:
        return val_to_return, clashing_chains
    else:
        return val_to_return, None


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


def superimpose_and_rotate(eq_chain1, eq_chain2, moving_chain, curr_struct, struct2):
    # all residues from same chain (common chain) are retrieved from the 2 structures. Example: chain A

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

    # add to the fixed structure, the moved chain
    added = 0
    clash, clashing_chains = is_Steric_clash(curr_struct, moving_chain)

    if not clash and rms <= 3.0:
        my_id = moving_chain.id
        chain_names = [x.id for x in curr_struct[0].get_chains()]
        while added == 0:
            rand = '|||' + create_random_chars_id(6)  # random ID
            if my_id + rand not in chain_names:
                moving_chain.id = my_id + rand
                curr_struct[0].add(struct2[0][moving_chain.id])
                added = 1

    return curr_struct, added, clash, clashing_chains, moving_chain


def generate_new_permutations(all_files, filename):

    """ This function inputs a list of files (all_files) and a filename of interest
    outputs a list of tuples, in which each of them has filename in a different position """

    new_permutations = []

    for idx in range(0, len(all_files)):

        all_files.remove(filename)
        all_files.insert(idx, filename)

        new_permutations.append(tuple(all_files))

    return new_permutations


def structure_in_created_structures(structure, created_structures):

    """ This function asks if the structure (with many chains) is already in  created_structures (a list of structures). returning a boolean if so"""

    # make a deepcopy of these objects
    structure = copy.deepcopy(structure)
    created_structures = copy.deepcopy(created_structures)

    # Get the ids of the chains in structure
    chain_ids_structure = tuple(sorted([x.id.split('|||')[1] for x in structure.get_chains()]))
    #chain_ids_structure = tuple(sorted([x.id for x in structure.get_chains()]))

    # loop through each of the contents of created_structures:
    for created_structure in created_structures:

        # get the chains of created_structure
        chain_ids_created_structure = tuple(sorted([x.id.split('|||')[1] for x in created_structure.get_chains()]))
        #chain_ids_created_structure = tuple(sorted([x.id for x in created_structure.get_chains()]))

        # ask if the number of each and ids of the chains are the same:
        if chain_ids_structure == chain_ids_created_structure:

            # create a list of the possible partners
            possible_partners = list(chain_ids_created_structure)

            # loop through all the chains in structure
            for chain_str in structure.get_chains():

                # a boolean that indicates if this chain in structure have a partner in created_strcuture
                chain_has_a_partner = False

                id_str = chain_str.id.split('|||')[1]
                #id_str = chain_str.id

                # try to find a partner in created_structure:
                for chain_created_str in created_structure.get_chains():

                    id_created_str = chain_created_str.id.split('|||')[1]
                    #id_created_str = chain_created_str.id

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
                        if rms <= 1.0:

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
                return True # this is only going to be true at the end if this is the last iteration

    # if you didn't find any match return false:
    return False

class TwoChainException(Exception):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "Input file %s does not have 2 chains" % self.file


def build_complex(saved_models, current_str, mydir, PDB_dict, Seq_to_filenames, this_is_a_branch=False, this_is_a_complex_recursion=False, non_brancheable_clashes=set(), tried_operations=set(), rec_level_branch=0, rec_level_complex=0, tried_branch_structures=list(), stoich=None):

    """This function builds a complex from a set of templates """

    # update the rec_level:

    if this_is_a_branch:
        rec_level_branch += 1

    if this_is_a_complex_recursion:
        rec_level_complex += 1

    print('BRANCH LEVEL: ', rec_level_branch, 'COMPLEX LEVEL: ', rec_level_complex)

    # a parser that is going to be used many times:
    p = pdb.PDBParser(PERMISSIVE=1)

    # files that have to be used
    all_files = list(PDB_dict.keys())

    if len(all_files) < 2:
        raise Exception("Only one file provided: no complex can be built.")

    # a boolean that indicates if there's something added at this level
    something_added = False

    # iterate through the chains of the current structure
    for chain1 in current_str[0].get_chains():

        id_chain1 = chain1.id.split('|||')[1]

        # iterate through the PDB files of the directory to compare them to the current struct
        for filename2 in all_files:
            rotating_chain = None
            common_chain2 = None

            # if the chain id of the current structure is found in the second PDB file:
            if id_chain1 in [x.split('|||')[1] for x in PDB_dict[filename2]]:
                # get the structure
                structure2 = p.get_structure("pr2", mydir + filename2)

                if len(list(structure2.get_chains())) != 2:
                    raise TwoChainException(filename2)

                # change the chain id names
                for chain in structure2.get_chains():
                    curr_id = chain.id
                    chain.id = [x for x in PDB_dict[filename2] if x[0] == curr_id][0]

                for chain2 in structure2.get_chains():
                    id_chain2 = chain2.id.split('|||')[1]

                    # if it is a homodimer
                    if PDB_dict[filename2][0].split('|||')[1] == PDB_dict[filename2][1].split('|||')[1]:
                        rotating_chain = structure2[0][PDB_dict[filename2][0]]
                        common_chain2 = structure2[0][PDB_dict[filename2][1]]

                    elif id_chain1 == id_chain2:
                        common_chain2 = chain2

                    else:
                        rotating_chain = chain2

                # check if this operation has already been tried before, and skip if so
                operation = (chain1.id, filename2, rotating_chain.id[0])

                if operation not in tried_operations:

                    # tried_operations.add(operation)
                    current_str, sth_added, clash, clashing_chains, added_chain = superimpose_and_rotate(chain1, common_chain2, rotating_chain, current_str, structure2)

                    # when something is added
                    if sth_added == 1:

                        # record if there has been any chain added:
                        something_added = True

                    # when there's a aberrant clash
                    if clash == 1:

                        # a branch complex will be created if the rotating chain is not one of the previously branch-opening clashing chains
                        # or if the clashes happen against the chain that opened this branch

                        open_branch = False

                        if this_is_a_branch:

                            clash_ids = set([(x, added_chain.id.split('|||')[1]) for x in clashing_chains])

                            # check if the clashes are not an exception
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            # set the id of the chain you were trying to add:
                            added_chain = copy.deepcopy(added_chain)
                            added_chain.id = added_chain.id + '|||' + create_random_chars_id(6)

                            # make a new branch:
                            branch_new_str = copy.deepcopy(current_str)

                            # remove the chains that chains that are clashing:
                            for clashing_chain in clashing_chains:
                                branch_new_str[0].detach_child(clashing_chain)

                            # create a set of tuples (added_chain, clashing_chains_ids) of branch opening exceptions
                            for x in clashing_chains:
                                non_brancheable_clashes.add((added_chain.id, x.split('|||')[1]))

                            # add the chain that was clashing:
                            branch_new_str[0].add(added_chain)

                            # check if the branch_new_str is not in tried_branch_structures:
                            if not structure_in_created_structures(branch_new_str, tried_branch_structures):

                                # add this branch to the ones already tested:
                                tried_branch_structures.append(branch_new_str)

                                # create a new structure based on this branch:
                                build_complex(saved_models, branch_new_str, mydir, PDB_dict, Seq_to_filenames, this_is_a_branch=True, non_brancheable_clashes=non_brancheable_clashes,
                                              rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich)

    if something_added:

        if this_is_a_branch:
            set_this_is_a_branch = True
        else:
            set_this_is_a_branch = False

        build_complex(saved_models, current_str, mydir, PDB_dict, Seq_to_filenames, this_is_a_complex_recursion=True, this_is_a_branch=set_this_is_a_branch, non_brancheable_clashes=non_brancheable_clashes,
                      rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich)

    else:
        print("trying to save model")

        if stoich:
            final_stoich = {}
            for chain in current_str.get_chains():
                chain_id = chain.id.split('|||')[1]

                if chain_id not in final_stoich:
                    final_stoich[chain_id] = 1
                else:
                    final_stoich[chain_id] += 1

            divisors = set()
            for key, value in final_stoich.items():
                divisor = value / stoich[key]
                divisors.add(divisor)

            if len(divisors) == 1 and not structure_in_created_structures(current_str, saved_models):
                saved_models.append(current_str)
                print("saved models: ", saved_models)

        elif not structure_in_created_structures(current_str, saved_models):
            saved_models.append(current_str)
            print("saved models: ", saved_models)

    return saved_models

