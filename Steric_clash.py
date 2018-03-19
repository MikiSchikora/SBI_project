
import Bio.PDB as pdb
import copy as cp
import rmsd
import numpy as np

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

def is_Steric_clash(structure, rotating_chain , distance_for_clash = 1.3):

    """This function returns False if there's no clash between a rotating chain and the current structure.
    If there's a clash it can return 1 (clash between two different chains) or 2 (you are trying to superimpose something in the place it was already)

    The clash crteria is that at least one of the atoms are at a lower distance than distance_for_clash A"""

    NS = pdb.NeighborSearch(list(structure.get_atoms()))

    clashing_chains = set()
    n_clashes = 0

    for at in rotating_chain.get_atoms():
        neighbors = NS.search(at.get_coord(), distance_for_clash)

        if len(neighbors)>0:

            for neigh in neighbors:

                clashing_chains.add(neigh.get_parent().get_parent().id)
                n_clashes += 1


    print(clashing_chains,n_clashes)

    if len(clashing_chains)>1:

        # a clash against different chains:
        return 1

    elif len(clashing_chains)==1 and n_clashes>100:

        # a clash MAYBE because you are trying to superimpose something in the place it was already

        #define the clashing chain:
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

        if RMSD <= 0.01:

            # it is the same chain
            return 2

        else:

            # it is anoother chain
            return 1

    else:

        # no clash
        return False




# check:
parser = pdb.PDBParser(PERMISSIVE=1)
structure1 = parser.get_structure('pdb_name', '3kuy.pdb')

rotating_chain = structure1[0]['C']

structure1[0]['C'].id = "C|||GAGSH|||HDHDJSSHJS"
rotating_chain.id = "C|||GAGSH|||HDHDHJS"

a = is_Steric_clash(structure1,rotating_chain)
print(a)

