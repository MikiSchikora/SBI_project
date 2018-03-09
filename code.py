import Bio.PDB as pdb


# argument parser
# input
# output
# verbose

#




# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# parse file 1 and 2 and get a structure for each
filename1 = "TEMPLATES/A_and_B.pdb"
structure1 = p.get_structure("pr1", filename1)
filename2 = "PAIR_AC.pdb"
structure2 = p.get_structure("pr2", filename2)

# same chain is retrieved from the 2 structures. Example: chain A
chainAfromAB=structure1[0]['A']
chainAfromAC=structure2[0]['A']

# get the atoms of the common chain in a list
atomsAfromAB=list(chainAfromAB.get_atoms())
atomsAfromAC=list(chainAfromAC.get_atoms())

# use the Superimposer
sup = pdb.Superimposer()

# first argument is fixed, second is moving. both are lists of Atom objects
sup.set_atoms(atomsAfromAB,atomsAfromAC)
print(sup.rotran)
print(sup.rms)

# rotate moving atoms
sup.apply(list(structure2[0]['C'].get_atoms()))

# add to the fixed structure, the moved chain
structure1[0].add(structure2[0]['C'])

# save in a pdb file
io = pdb.PDBIO()
io.set_structure(structure1)
io.save('out1.pdb')

