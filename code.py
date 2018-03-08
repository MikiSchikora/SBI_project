import Bio.PDB as pdb

p = pdb.PDBParser(PERMISSIVE=1)
filename1 = "prova1.pdb"
structure1 = p.get_structure("pr1", filename1)
filename2 = "prova2.pdb"
structure2 = p.get_structure("pr2", filename2)

structureAB = pdb.Model.Model('AB')
structureAB.add(structure1[0]['A'])

print(type(structureAB))
for chain in structureAB.get_chains():
   print(chain)

structureAC = pdb.Model.Model('AC')
structureAC.add(structure2[0]['A'])


for chain in structureAB.get_chains():
   if chain.get_id() == 'A':
      AB = list(chain.get_atoms())

for chain in structureAC.get_chains():
   if chain.get_id() == 'A':
      AC = list(chain.get_atoms())

print(AB)

for atom in AB:
   #print(dir(atom))
   print(atom.get_coord())
   break
   #print(atom)


sup = pdb.Superimposer()

# first argument is fixed, second is moving. both are lists of Atom objects
sup.set_atoms(AC,AB)
print(sup.rotran)
print(sup.rms)

# rotate moving atoms
sup.apply(AB)


io = pdb.PDBIO()
for atom in AB:
   print(atom.get_coord())
   a=atom.get_parent()
   #print(list(a.get_atoms()))
   # io.save('out.pdb')
for atom in AC:
   print(atom.get_coord())
   a=atom.get_parent()
