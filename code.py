from Bio.PDB import *

p = PDBParser(PERMISSIVE=1)
filename = "3kuy.pdb"
structure = p.get_structure("3kuy", filename)

print(dir(structure))
print(type(structure))


chains = structure.get_chains()
print(chains)
print(dir(chains))

for chain in chains:
   print(dir(chain))
   if chain.get_id() == 'B':
      AB = list(chain.get_atoms())
   if chain.get_id() == 'B':
      BC = list(chain.get_atoms())

#print(AB)

for atom in AB:
   #print(dir(atom))
   print(atom.get_coord())
   break
   #print(atom)


sup = Superimposer()

# first argument is fixed, second is moving. both are lists of Atom objects
sup.set_atoms(AB,BC)
print(sup.rotran)
print(sup.rms)

# rotate moving atoms
sup.apply(BC)


io = PDBIO()
for atom in BC:
   print(type(atom))
   io.set_structure(atom)
io.save('out.pdb')

#print(residuesBC)
