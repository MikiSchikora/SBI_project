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
import itertools

# some variables:

chain_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
            "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"
            , "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "0", "1", "2", "3", "4"
            , "5", "6", "7", "8", "9"]

#chain_alphabet = ["A", "B", "C", "D"]


repeat2 = [''.join(p) for p in itertools.product(chain_alphabet, repeat=2)]
repeat3 = [''.join(p) for p in itertools.product(chain_alphabet, repeat=3)]

complete_chain_alphabet = chain_alphabet + repeat2 + repeat3

def Generate_pairwise_subunits_from_pdb(pdb_file_path, TEMPLATES_path, file_type, verbose):



    """This function takes an existing complex and fragments it into each of the pairwise interactions between subunits.

    pdb_file_path is the path where the complex PDB is

    TEMPLATES_path is the folder where the resulting folders have to be saved

    file_type indicates which is the type of file

    verbose indicates if a log of the program execution has to be saved

    It does not consider nucleic acid sequences, as this is only for testing the program on different complexes"""

    num_file = 0

    if file_type == 'PDB':
        parser = pdb.PDBParser(PERMISSIVE=1)
    else:
        parser = pdb.MMCIFParser()

    structure = parser.get_structure('pdb_name', pdb_file_path)

    # give unique chain identifiers to a structure, it has to be similar to the ids of the chains used in build_complex, to be able to use further the structure_in_created_structures() function
    id_nch = 0
    for chain in structure.get_chains():
        actual_id = chain.id
        chain.id = (complete_chain_alphabet[id_nch]+'_',actual_id)
        id_nch += 1

    # free the ./TEMPLATES_path/
    os.system('rm -rf ' + TEMPLATES_path + '*')

    # initialize the saved pairs and structures
    saved_pairs = set()
    saved_structures = []

    # loop through all possible pairwise files

    for chain1 in structure.get_chains():

        for chain2 in structure.get_chains():

            # the following strings define the pairs already saved
            comb = tuple(list(chain1.id) + list(chain2.id))
            comb_rev = tuple(list(chain2.id) + list(chain1.id))

            if chain1 is not chain2 and comb not in saved_pairs:

                # save the combination
                saved_pairs.add(comb)
                saved_pairs.add(comb_rev)

                # ask if any of the residues is interacting, if so save the PDB

                chains_interacting = False

                for residue1 in chain1:
                    if chains_interacting is True:
                        break
                    for residue2 in chain2:
                        if residue1 != residue2:

                            # define which is the important residue of each chain:
                            atoms1 = [x.id for x in residue1.get_atoms()]
                            atoms2 = [x.id for x in residue2.get_atoms()]

                            important_atom1 = None
                            if 'CA' in atoms1:
                                important_atom1 = residue1['CA']
                            elif 'P' in atoms1:
                                important_atom1 = residue1['P']

                            important_atom2 = None
                            if 'CA' in atoms2:
                                important_atom2 = residue2['CA']
                            elif 'P' in atoms2:
                                important_atom2 = residue2['P']

                            # compute the distance:
                            if important_atom1 is not None and important_atom2 is not None:
                                distance = important_atom1 - important_atom2
                            else:
                                continue

                            if distance < 7:
                                chains_interacting = True
                                break

                if chains_interacting is True:

                    # create a structure object
                    ID = str(num_file)
                    num_file += 1
                    new_structure = pdb_struct.Structure(ID)

                    new_model = pdb_model.Model(0)
                    new_model.add(copy.deepcopy(chain1))
                    new_model.add(copy.deepcopy(chain2))

                    new_structure.add(new_model)

                    # move the coordinates of the structure to simulate what would happen if they were coming from different files
                    rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
                    translation = np.array((0, 0, 1), 'f')
                    for atom in new_structure.get_atoms():
                        atom.transform(rotation, translation)

                    # write to new pdb:

                    if structure_in_created_structures(new_structure,saved_structures) is False:

                        # record as a saved structure:
                        saved_structures.append(new_structure)

                        # give unique chains to a structure (A and B)
                        id_nch = 0
                        for chain in new_structure.get_chains():
                            chain.id = chain_alphabet[id_nch]
                            id_nch += 1

                        if verbose:
                            print('writing PDB file with the interaction of %s and %s into %s.pdb'%(chain1.id[1],chain2.id[1],ID))

                        # write using our customized writer
                        io = pdb.PDBIO()
                        io.set_structure(new_structure)
                        io.save(TEMPLATES_path + ID + '.pdb')


def create_random_chars_id(n):

    """Creates an identifier with n random characters."""

    rand_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))
    return rand_id


def add_new_seqid_to_dicts(seq, dict, n):

    """Adds a random identifier of n characters to 2 dictionaries, the first having a sequence, the other initialising and empty set"""

    seq_id = create_random_chars_id(n)  # create random ID
    while seq_id in dict:
        seq_id = create_random_chars_id(n)  # random ID is already in Seqs dict
    dict[seq_id] = seq

    return seq_id, dict[seq_id]


def Generate_PDB_info(Templates_dir, subunits_seq_file, verbose, min_identity_between_chains=95):

    """This function takes the Templates_dir and creates a dictionary with information about each template
    Each template gets incorporated into a dictionary with a unique identifiers for it's chains, in a way that you save all files that are useful

    min_identity_between_chains refers to the minimum identity required between two molecules to state that they are the same

    subunits_seq_file is a fasta file containing all the molecules that form the complex """

    # each chain in the template will have a unique identifier, composed by a tuple of (chain accession , molecule identifier). Accession is the provided chain id, and molecule identifier indicates which is the molecule of that chain

    List_PDBs = os.listdir(Templates_dir)
    PDB_info = {}  # This will contain information about the unique chains: Keys: filename, Values: [(chain accession , molecule identifier) ... ]

    ppb = pdb.Polypeptide.PPBuilder()  # polypeptide parser
    parser = pdb.PDBParser(PERMISSIVE=1)  # pdb parser

    # create a dictionary that has the known sequences if a multifasta is given
    Seqs = {}
    if subunits_seq_file:
        for record in Seq_IO.parse(subunits_seq_file, "fasta"):
            Seqs[record.id] = record.seq._data

    for file in List_PDBs:

        if verbose:
            print("Processing the content of file %s"%(file))


        # generate the structure of this file:
        path = Templates_dir + file
        PDB_structure = parser.get_structure('structure', path)
        chains = list(PDB_structure.get_chains())
        Chain_IDs = []  # wll contain the unique ID for each chain

        for chain in chains:

            # obtain the sequence of the chain
            molSeq = "".join([str(x.get_sequence()) for x in ppb.build_peptides(chain)])

            # some chains are not properly characterized and have unassignable residues (UNK). we drop them if so
            if molSeq == '':

                # it may be nucleic acids
                molSeq = ''.join([x.resname[2] for x in chain.get_residues() if x.resname[2] in 'ACGTU' and len(x.resname)==3])

                if not len(molSeq)>0:
                    break

            # assign an ID (id of the fasta or a new for the molecule) to the sequence:
            Seq_id = None
            Best_identity = min_identity_between_chains

            # compare the current chain sequence with the ones present in Seqs
            if Seqs:
                for my_id in Seqs:
                    # align my current sequence with the sequences stored in Seqs dictionary
                    alignment = pairwise2.align.globalxx(molSeq, Seqs[my_id])
                    # compute percentage identity
                    percent_match = (alignment[0][2] / min(len(molSeq), len(Seqs[my_id]))) * 100

                    # if percentage identity is lar
                    if percent_match > Best_identity:
                        Seq_id = my_id
                        Best_identity = percent_match

                # my current sequence didn't match anything in Seqs, so it is added to the dictionary
                if Seq_id is None and subunits_seq_file is None:
                    Seq_id, Seqs[Seq_id] = add_new_seqid_to_dicts(molSeq, Seqs, 6)

            # add the current chain sequence to Seqs if it is empty
            else:
                Seq_id, Seqs[Seq_id] = add_new_seqid_to_dicts(molSeq, Seqs, 6)

            if Seq_id is None:  # only if a multifasta file is provided and this chain is not there
                print("WARNING: %s contains non expected sequences!! It will not be used for building the model."%(file))
                break

            else:
                ID_new = (chain.id,Seq_id)
                Chain_IDs.append(ID_new)

        # save info into dictionary PDB_info
        if len(Chain_IDs) == 2:
            PDB_info[file] = Chain_IDs

    return PDB_info, Seqs


def is_Steric_clash(structure, rotating_chain, distance_for_clash=1.4):

    """This function returns False if there's no clash between a rotating chain and the current structure.
    If there's a clash it can return 1 (clash between two different chains) or 2 (you are trying to superimpose something in the place it was already)

    it also returns the ids of the chains in structure that are clashing

    The clash criteria is that at least one 10 the atoms are at a lower distance than distance_for_clash Amstrongs.
    The logic behind is that few atoms clashing indicate that the resulting structure could be viable after optimization with the rosetta software.

    It is considered that two chains are in the same place if the RMSD between them is lower or equal to 3.0. """

    # initialize the neighbor search
    NS = pdb.NeighborSearch(list(structure.get_atoms()))

    clashing_chains = set() # the set of clashing chains (in structure)
    n_clashes = 0 # the number of clashes

    for at in rotating_chain.get_atoms():
        neighbors = NS.search(at.get_coord(), distance_for_clash)

        if len(neighbors) > 0:
            for neigh in neighbors:
                clashing_chains.add(neigh.get_parent().get_parent().id)
                n_clashes += 1

    if len(clashing_chains) > 1 and n_clashes > 20:
        # a clash against different chains:
        val_to_return = 1

    elif len(clashing_chains) == 1 and n_clashes > 20 and rotating_chain.id[1] == list(clashing_chains)[0][1]:

        # a clash MAYBE because you are trying to superimpose something in the place it was already

        # define the clashing chain:
        clash_chain = structure[0][list(clashing_chains)[0]]
        res_chain1 = list(clash_chain.get_residues())
        res_chain2 = list(rotating_chain.get_residues())

        # so first we obtain a list of the common residues
        common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
        common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

        # then we obtain a list of atom objects to use it later
        common_atoms_s1 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s1)])
        common_atoms_s2 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_s2)])

        RMSD = rmsd.kabsch_rmsd(common_atoms_s1, common_atoms_s2)

        if RMSD <= 3.0:
            # it is the same chain
            val_to_return = 2

        else:
            # it is another chain or the same with different structure
            val_to_return = 1

    elif n_clashes > 20:
        # it is ine chain and the previous conditions are not fullfilled
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

    # define the set of ids of the valid residues in reslist2, they have to include 'CA' and 'P' in both
    ids_res_2 = set([x.id[1] for x in reslist2 if ('CA' in x or 'P' in x)])

    common_res = [x for x in reslist1 if (x.id[1] in ids_res_2 and ('CA' in x or 'P' in x))]

    return common_res


def get_atom_list_from_res_list(reslist):

    """Get a list of atom objects from a list of residue objects"""

    atomlist = []
    for res in reslist:
        for atom in res.get_atoms():
            if atom.id == 'CA' or atom.id == 'P':
                atomlist.append(atom)
                break

    return atomlist


def superimpose_and_rotate(eq_chain1, eq_chain2, moving_chain, curr_struct, rec_level_complex, filename2=None):

    """This function inputs several objects for performing the adding of chains to the current complex.

    eq_chain1 is the common chain in the current structure (curr_struct)

    eq_chain2 is the common chain in the structure from which you want to add a chain

    moving_chain is the chain that may be added to the current complex

    rec_level_complex is the recursion level of building the complex

    filename2 is the name of the file that contained the moving_chain

    It returns the current structure, information about the clashes and a flag for having added something"""

    # all residues from same chain (common chain) are retrieved from the 2 structures. Example: chain A

    res_chain1 = list(eq_chain1.get_residues())
    res_chain2 = list(eq_chain2.get_residues())

    # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
    # so first we obtain a list of the common residues
    common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
    common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

    # then we obtain a list of atom objects to use it later
    common_atoms_s1 = get_atom_list_from_res_list(common_res_s1)
    common_atoms_s2 = get_atom_list_from_res_list(common_res_s2)


    # debug
    if len(common_atoms_s1)!=len(common_atoms_s2):
        return curr_struct, 0, False, set(), moving_chain
    

    # use the Superimposer
    sup = pdb.Superimposer()

    # first argument is fixed, second is moving. both are lists of Atom objects
    sup.set_atoms(common_atoms_s1, common_atoms_s2)
    rms = sup.rms

    # rotate moving atoms
    sup.apply(list(moving_chain.get_atoms()))

    # add to the fixed structure, the moved chain
    added = 0
    clash, clashing_chains = is_Steric_clash(curr_struct, moving_chain)

    # something is added if there's no clashes and the RMSD is very low, indicating that the two chains are actually the same
    if not clash and rms <= 3.0:
        my_id = moving_chain.id
        chain_names = [x.id for x in curr_struct[0].get_chains()]
        while added == 0:
            rand = (create_random_chars_id(6),str(rec_level_complex)) # random ID + a number that indicates the recursion level at which this chain has been added
            if my_id + rand not in chain_names:
                moving_chain.id = tuple(list(my_id) + list(rand))
                curr_struct[0].add(moving_chain)
                added = 1

    return curr_struct, added, clash, clashing_chains, moving_chain


def structure_in_created_structures(structure, created_structures):

    """ This function asks if the structure (with many chains) is already in  created_structures (a list of structures). returning a boolean if so.

    This is considered to be true if all of the chains in structure have a chain in one of the structures in created_structures that have an RMSD of 3.0 or less.

    This is compared after superposing structure and one of the structures in created_structures by a common chain (below called )"""

    # make a deepcopy of these objects
    structure = copy.deepcopy(structure)
    created_structures = copy.deepcopy(created_structures)

    # Get the ids of the chains in structure
    chain_ids_structure = tuple(sorted([x.id[1] for x in structure.get_chains()]))

    # loop through each of the contents of created_structures:
    for created_structure in created_structures:

        # get the chains of created_structure
        chain_ids_created_structure = tuple(sorted([x.id[1] for x in created_structure.get_chains()]))

        # ask if the number of each and ids of the chains are the same:
        if chain_ids_structure == chain_ids_created_structure:

            # pick one chain in structure to compare with the chains in created
            chain_str = list(structure.get_chains())[0]
            id_str = chain_str.id[1]

            # try to find a partner in created_structure:
            for chain_created_str in created_structure.get_chains():

                id_created_str = chain_created_str.id[1]

                # if they have the same id they are potential partners. Superimpose these next. The id_created_str has also to be avaliable in possible_partners
                if id_str == id_created_str:

                    # get list of residues:
                    res_chain1 = list(chain_str.get_residues())
                    res_chain2 = list(chain_created_str.get_residues())

                    # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
                    # so first we obtain a list of the common residues
                    common_res_s1 = get_list_of_common_res(res_chain1, res_chain2)
                    common_res_s2 = get_list_of_common_res(res_chain2, res_chain1)

                    # then we obtain a list of atom objects to use it later
                    common_atoms_s1 = get_atom_list_from_res_list(common_res_s1)
                    common_atoms_s2 = get_atom_list_from_res_list(common_res_s2)

                    # continue if the common atoms is full
                    if len(common_atoms_s1) > 0:

                        # use the Superimposer
                        sup = pdb.Superimposer()

                        # first argument is fixed, second is moving. both are lists of Atom objects
                        sup.set_atoms(common_atoms_s2, common_atoms_s1)

                        # if I have superimposed same ID but different structure, try another chain
                        if sup.rms > 3.0:
                            continue

                        # apply rotation to whole common structure
                        sup.apply(list(structure.get_atoms()))

                        # if the previous chain_str and chain_created_str are real partners they should also result in haveing all they cross-superimposed chains with partners
                        partners = set()

                        for searching_partner in created_structure.get_chains():
                            partner_found = False

                            for possible_partner in structure.get_chains():

                                if partner_found is True:
                                    break

                                if possible_partner.id[1] == searching_partner.id[1] and possible_partner not in partners:
                                    # get list of residues:
                                    res_partner1 = list(searching_partner.get_residues())
                                    res_partner2 = list(possible_partner.get_residues())

                                    # get the atoms of the previous list, ONLY belonging to common RESIDUES! to be then able to superimpose
                                    # so first we obtain a list of the common residues
                                    common_res_p1 = get_list_of_common_res(res_partner1, res_partner2)
                                    common_res_p2 = get_list_of_common_res(res_partner2, res_partner1)

                                    # then we obtain a list of coordinates
                                    common_coords_p1 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_p1)])
                                    common_coords_p2 = np.array([list(x.get_coord()) for x in get_atom_list_from_res_list(common_res_p2)])

                                    rms = rmsd.kabsch_rmsd(common_coords_p2, common_coords_p1)
                                    if rms <= 3.0:
                                        partners.add(possible_partner)
                                        partner_found = True

                        if len(partners) == len(list(created_structure.get_chains())):
                            return True  # all chains have a partner, which means that the structure is in the created_structures

    # if you didn't find any match return false:
    return False


class TwoChainException(Exception):
    def __init__(self, file):
        self.file = file

    def __str__(self):
        return "Input file %s does not have 2 chains. We reccomend to remove it from the templates directory, as this pipeline requires templates to be PAIRWISE interactions." % self.file


def print_topology_of_complex(structure):

    """This function prints to STDOUT the topology of a structure"""

    print('The current complex has %i chains, with the molecules:'%(len(list(structure.get_chains()))))

    # create a dictionary of the molecules
    chains_dict = {}
    for chain_in_comp in structure.get_chains():

        id_ch = chain_in_comp.id[1]
        if id_ch not in chains_dict:
            chains_dict[id_ch] = 1
        else:
            chains_dict[id_ch] += 1

    # print
    for key, val in chains_dict.items():
        print(key, ': ', str(val))


def build_complex(saved_models, current_str, mydir, PDB_dict, num_models, exhaustive, n_chains, this_is_a_branch=False, this_is_a_complex_recursion=False, non_brancheable_clashes=set(), rec_level_branch=0, rec_level_complex=0, tried_branch_structures=list(), stoich=None, verbose=False):

    """This function builds the complex. It has many recursion-opening possibilities, returning saved_models (a list of the currently saved models) in each of them.

    current_str is the structure being built

    mydir contains the files with the template pairwise interactions

    PDB_dict is the dictionary with the information about the chain contents of each file

    num_models is the desired number of different models to generate with this function

    exhaustive is a boolean that indicates if all possible complexes have to be generated. This is dangerous for large complexes.

    this_is_a_branch indicates if the current function is running in a branch of complex building.

    this_is_a_complex_recursion indicates if the current function is trying to add recursively new chains to the model previosuly created

    non_brancheable_clashes is a set of tuples that correspond to clashes that should not open a branch.

    rec_level_branch indicates which is the level of the branch opened.

    rec_level_complex indicates which is in which recursion level, in terms of adding new chains, is the current function.

    tried_branch_structures indicates which are the structures that opened a branch.

    stoich is a dictionary that includes information about the expected stoichiometry of the complex.

    verbose is a boolean that indicates if the progression of the program has to be considered.

    NOTE THAT if exhaustive is true,  stoich is provided or n_models exceeds 1 this function will open branches for building more than one complex."""

    # return as soon as possible:
    if exhaustive is False and num_models==len(saved_models):
        return saved_models

    # update the rec_level:

    if this_is_a_branch: # level of the recursions opened by the clash of a chain against another
        rec_level_branch += 1

    if this_is_a_complex_recursion:  # level of the recursions opened by adding subunits
        rec_level_complex += 1

    if verbose:
        print('The current branch level is: ',rec_level_branch, 'The current complex building level: ',rec_level_complex)

        if exhaustive:
            print('The number of complexes already created is %i and you want to generate all possible complexes'%(len(saved_models)))
        else:
            print('The number of complexes already created is %i and you want to generate %i'%(len(saved_models),num_models))

    # a parser that is going to be used many times:
    p = pdb.PDBParser(PERMISSIVE=1)

    # files that have to be used. The shuffle is for generating alternative paths, so that running the program many times may generate slightly different resukts
    all_files = list(PDB_dict.keys())
    random.shuffle(all_files)

    if len(all_files) < 2:
        raise Exception("Only one file provided: no complex can be built.")

    # a boolean that indicates if there's something added at this level
    something_added = False

    # print the conformation of the current complex:
    if verbose:
        print_topology_of_complex(current_str)

    current_chains = list(current_str.get_chains())

    # iterate through the chains of the current structure
    for chain1 in current_chains:

        # when you are in deeper recursion levels check that the chains are the ones of the previous level, for avoiding many comparisons
        interesting_complex_rec_level = None

        if this_is_a_complex_recursion:
            # we are interested in analyzing chains that have been added JUST in the previous recursion levels
            interesting_complex_rec_level = rec_level_complex - 1

        elif this_is_a_branch:
            # we are interested in analyzing chains of the current rec_level_complex
            interesting_complex_rec_level = rec_level_complex

        # remember that chain1.id will be atuple of ONLY two elements if rec_level_complex=0
        if rec_level_complex>0 and chain1.id[-1]!=str(interesting_complex_rec_level):
            continue

        if verbose:
            print('Trying to add chains through the common chain: ', chain1.id)

        id_chain1 = chain1.id[1]

        # iterate through the PDB files of the directory to compare them to the current struct
        for filename2 in all_files:

            # initialize the common and rotating chain
            rotating_chain = None
            common_chain2 = None


            # if the chain id of the current structure is found in the second PDB file:
            if id_chain1 in [x[1] for x in PDB_dict[filename2]]:

                # get the structure
                structure2 = p.get_structure("pr2", mydir + filename2)

                if len(list(structure2.get_chains())) != 2:
                    raise TwoChainException(filename2)

                # change the chain id names
                for chain in structure2.get_chains():
                    curr_id = chain.id # the so-called chain accession, usually a letter (A,B,C,D ...)
                    chain.id = [x for x in PDB_dict[filename2] if x[0] == curr_id][0]


                # if it is a homodimer set both pssibilities of rotating / common
                if PDB_dict[filename2][0][1] == PDB_dict[filename2][1][1]:

                    # generate an array with the two possibilities

                    rotating_chains = [structure2[0][PDB_dict[filename2][0]],structure2[0][PDB_dict[filename2][1]]]
                    common_chains2 = [structure2[0][PDB_dict[filename2][1]],structure2[0][PDB_dict[filename2][0]]]


                else:
                    # decide which is (rotating / common)
                    for chain2 in structure2.get_chains():
                        id_chain2 = chain2.id[1]

                        if id_chain1 == id_chain2:
                            common_chain2 = chain2
                        else:
                            rotating_chain = chain2

                    # define the array
                    rotating_chains = [rotating_chain]
                    common_chains2 = [common_chain2]

                #go through each possible rotating and common chains

                for I in range(0,len(rotating_chains)):

                    # define each of the rotating/common, without perturbing the original ones
                    if len(rotating_chains)>1:
                        rotating_chain = copy.deepcopy(rotating_chains[I])
                        common_chain2 = copy.deepcopy(common_chains2[I])
                    else:
                        rotating_chain = rotating_chains[I]
                        common_chain2 = common_chains2[I]

                    # try to add the new chain to the complex
                    current_str, sth_added, clash, clashing_chains, added_chain = superimpose_and_rotate(chain1, common_chain2, rotating_chain, current_str, rec_level_complex, filename2=filename2)

                    # when something is added, in this try
                    if sth_added == 1:
                        # record if there has been any chain added:
                        something_added = True

                    # when there's a aberrant clash and you fulfill one of the branch-opening  conditions
                    if clash == 1 and (exhaustive or num_models>1 or stoich):

                        # a branch complex will be created if the rotating chain is not one of the previously branch-opening clashing chains
                        # or if the clashes happen against the chain that opened this branch

                        open_branch = False

                        if this_is_a_branch:

                            clash_ids = set([(x, added_chain.id[1]) for x in clashing_chains])

                            # check if the clashes are not an exception
                            if len(non_brancheable_clashes.intersection(clash_ids)) == 0:
                                open_branch = True

                        else:
                            open_branch = True

                        if open_branch:

                            # set the id of the chain you were trying to add:
                            added_chain = copy.deepcopy(added_chain)
                            added_chain.id = tuple(list(added_chain.id) + [create_random_chars_id(6), str(rec_level_complex)])

                            # make a new branch:
                            branch_new_str = copy.deepcopy(current_str)

                            # remove the chains that chains that are clashing:
                            for clashing_chain in clashing_chains:
                                branch_new_str[0].detach_child(clashing_chain)

                            # create a set of tuples (added_chain, clashing_chains_ids) of branch opening exceptions
                            for x in clashing_chains:
                                non_brancheable_clashes.add((added_chain.id, x[1]))

                            # add the chain that was clashing:
                            branch_new_str[0].add(added_chain)

                            # check if the branch_new_str is not in tried_branch_structures:
                            if structure_in_created_structures(branch_new_str, tried_branch_structures) is False:

                                # indicate that a branch is opening
                                if verbose:
                                    print('%s of %s  clashes against the complex. Opening a new branch...' %(added_chain.id[1],filename2))

                                # add this branch to the ones already tested:
                                tried_branch_structures.append(branch_new_str)

                                # create a new structure based on this branch:
                                saved_models = build_complex(saved_models, branch_new_str, mydir, PDB_dict, num_models, exhaustive, n_chains, this_is_a_branch=True, non_brancheable_clashes=non_brancheable_clashes,
                                              rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

                                # return as soon as possible:
                                if exhaustive is False and num_models == len(saved_models):
                                    return saved_models

    if something_added and len(list(current_str.get_chains()))<=n_chains:

        if this_is_a_branch:
            set_this_is_a_branch = True
        else:
            set_this_is_a_branch = False

        saved_models = build_complex(saved_models, current_str, mydir, PDB_dict, num_models, exhaustive, n_chains, this_is_a_complex_recursion=True, this_is_a_branch=set_this_is_a_branch, non_brancheable_clashes=non_brancheable_clashes,
            rec_level_complex=rec_level_complex, rec_level_branch=rec_level_branch, tried_branch_structures=tried_branch_structures, stoich=stoich, verbose=verbose)

        # return as soon as possible:
        if exhaustive is False and num_models == len(saved_models):
            return saved_models

    else:
        if verbose:
            print("trying to save model")

        # check stoichiomatry
        if stoich:
            final_stoich = {}
            for chain in current_str.get_chains():
                chain_id = chain.id[1]

                if chain_id not in final_stoich:
                    final_stoich[chain_id] = 1
                else:
                    final_stoich[chain_id] += 1

            divisors = set()
            for key, value in final_stoich.items():
                divisor = value / stoich[key]
                divisors.add(divisor)

            if len(divisors) == 1 and structure_in_created_structures(current_str, saved_models) is False:
                saved_models.append(current_str)

                if verbose:
                    print("saving model")
                    print_topology_of_complex(current_str)
                    print('\n\n')

        elif structure_in_created_structures(current_str, saved_models) is False:

            saved_models.append(current_str)
            if verbose:
                print("saving model")
                print_topology_of_complex(current_str)
                print('\n')

    return saved_models


