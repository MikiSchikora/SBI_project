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
import itertools as iter
import rmsd


def Generate_pairwise_subunits_from_pdb(pdb_file_path, Templates_dir):
    """This function takes an existing complex and fragments it into each of the pairwise interactions between subunits.
    pdb_file_path is the path where the complex PDB is
    It does not consider nucleic acid sequences, as this is only for testing the program on different complexes"""

    TEMPLATES_path = Templates_dir

    parser = pdb.PDBParser(PERMISSIVE=1)

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
                    ID = chain1.id + '_and_' + chain2.id
                    new_structure = pdb_struct.Structure(ID)

                    new_model = pdb_model.Model(0)
                    new_model.add(copy.deepcopy(chain1))
                    new_model.add(copy.deepcopy(chain2))
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
        PDB_info[file] = Chain_IDs

    return PDB_info, Seq_to_filenames


def is_Steric_clash(structure, rotating_chain, distance_for_clash=1):

    """This function returns False if there's no clash between a rotating chain and the current structure.
    If there's a clash it can return 1 (clash between two different chains) or 2 (you are trying to superimpose something in the place it was already)

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
        return 1

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
        common_atoms_s1 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s1) if x.id=='CA'])
        common_atoms_s2 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s2) if x.id=='CA'])

        RMSD = rmsd.kabsch_rmsd(common_atoms_s1, common_atoms_s2)
        print("clash chain: ", clash_chain, " rotating chain: ", rotating_chain, RMSD)

        if RMSD <= 1:
            # it is the same chain
            return 2

        else:
            # it is another chain or the same with different structure
            return 1

    elif n_clashes > 0:
        # it is another chain
        return 1

    else:
        # no clash
        return False


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
    # print(sup.rotran)
    # print(sup.rms)

    # rotate moving atoms
    sup.apply(list(struct2[0][moving_chain.id].get_atoms()))

    # add to the fixed structure, the moved chain
    added = 0
    clash = is_Steric_clash(curr_struct, moving_chain)

    if not clash:
        my_id = moving_chain.id
        chain_names = [x.id for x in curr_struct[0].get_chains()]
        while added == 0:
            rand = '|||' + create_random_chars_id(6)  # random ID
            if my_id + rand not in chain_names:
                moving_chain.id = my_id + rand
                curr_struct[0].add(struct2[0][moving_chain.id])
                added = 1
    elif clash == 1:
        pass
        print("clash=1")

    return curr_struct, added


def build_complex(current_str, mydir, PDB_dict, Seq_to_filenames):

    """This function builds a complex from a set of templates"""

    # call the global variables:

    # global tried_operations
    # global rec_level

    # rec_level += 1

    # print(Seq_to_filenames)

    sth_added = 0
    p = pdb.PDBParser(PERMISSIVE=1)

    # define a list of the files containing the chains in current_str
    filenames_all = []
    for chain1 in current_str[0].get_chains():

        id_chain1 = chain1.id.split('|||')[1]
        filenames_all += list(Seq_to_filenames[id_chain1])

    filenames_all = list(set(filenames_all))  # unique

    # print(filenames_all)
    n_perm = 0
    permutations = []
    for permutation in iter.permutations(filenames_all):
        permutations.append(permutation)
        n_perm += 1
        # print(permutation)
        if n_perm == 100:
            break

    # print(permutations)

    # iterate through the chains of the current structure
    for chain1 in current_str[0].get_chains():

        id_chain1 = chain1.id.split('|||')[1]

        # iterate through the PDB files of the directory to compare them to the current struct
        for filename2 in os.listdir(mydir):
            rotating_chain = None
            common_chain2 = None

            # if the chain id of the current structure is found in the second PDB file:
            if id_chain1 in [x.split('|||')[1] for x in PDB_dict[filename2]]:
                # get the structure
                structure2 = p.get_structure("pr2", mydir + filename2)

                # change the chain id names
                for chain in structure2.get_chains():
                    curr_id = chain.id
                    chain.id = [x for x in PDB_dict[filename2] if x[0] == curr_id][0]

                for chain2 in structure2.get_chains():
                    id_chain2 = chain2.id.split('|||')[1]

                    if PDB_dict[filename2][0].split('|||')[1] == PDB_dict[filename2][1].split('|||')[1]:
                        rotating_chain = structure2[0][PDB_dict[filename2][0]]
                        common_chain2 = structure2[0][PDB_dict[filename2][1]]

                    elif id_chain1 == id_chain2:
                        common_chain2 = chain2
                    else:
                        rotating_chain = chain2

                # check if this operation has already been tried before, and skip if so

                # operation = (chain1.id, filename2, rotating_chain.id[0])

                current_str, sth_added = superimpose_and_rotate(chain1, common_chain2, rotating_chain, current_str, structure2)

    if sth_added == 1:
        build_complex(current_str, mydir, PDB_dict, Seq_to_filenames)

    else:
        # we go through all the chains of the structure and rename them alphabetically
        final_chains = list(current_str.get_chains())  # list of chain objects
        chain_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                          "T", "U", "V", "W", "X", "Y", "Z"]
        alphabet_pos = 0
        for final_chain in final_chains:
            final_chain.id = chain_alphabet[alphabet_pos]
            alphabet_pos += 1
            # print("final chain id")
            # print(final_chain.id)

        # then we can finally save the obtained structure object into a pdb file
        io = pdb.PDBIO()
        io.set_structure(current_str)
        io.save('out1.pdb')

        return

    return
