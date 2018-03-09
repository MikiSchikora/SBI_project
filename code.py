import Bio.PDB as pdb
import os


# argument parser
# input
# output
# verbose

# if input is a single pdb and we want to split it
# use Generate_pairwise_interactions (creates a file with 2 subunits
# if they are at less than a given distance
# at the end we have a path (directory with all these pdb pairs


# generate PDB info (it is a dictionary of a dictionary)
# primary key: filename, secondary key: chain id, value: chain unique accession number (it is a number)




# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# parse file 1 and 2 and get a structure for each
filename1 = "TEMPLATES/A_and_B.pdb"
structure1 = p.get_structure("pr1", filename1)

common_chain='A'
for filename2 in os.listdir('TEMPLATES'):
    if filename2.startswith(common_chain) and filename2[6] != 'B':
        structure2 = p.get_structure("pr2", "TEMPLATES/"+filename2)

        rotating_chain = filename2[6]

        # same chain is retrieved from the 2 structures. Example: chain A
        common_chain_s1 = structure1[0][common_chain]
        common_chain_s2 = structure2[0][common_chain]

        # get the atoms of the common chain in a list
        common_chain_atoms_s1 = list(common_chain_s1.get_atoms())
        common_chain_atoms_s2 = list(common_chain_s2.get_atoms())

        # use the Superimposer
        sup = pdb.Superimposer()

        # first argument is fixed, second is moving. both are lists of Atom objects
        sup.set_atoms(common_chain_atoms_s1, common_chain_atoms_s2)
        print(sup.rotran)
        print(sup.rms)

        # rotate moving atoms
        sup.apply(list(structure2[0][rotating_chain].get_atoms()))

        # add to the fixed structure, the moved chain
        structure1[0].add(structure2[0][rotating_chain])

        # save in a pdb file
        io = pdb.PDBIO()
        io.set_structure(structure1)
        io.save('out1.pdb')
