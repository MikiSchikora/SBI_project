
import Bio.PDB as pdb
import copy as cp

def is_Steric_clash(structure, rotating_chain):

   """This function returns a boolean that indicates if a rotating_chain clashes against
    a structure.

    The clash crteria is that at least one of the atoms are at a lower distance than
    min_VDW_distance"""

   Van_der_Waals_radius = {'H':1.2,'C':1.7,'N':1.55,'O':1.52,'F':1.47,'P':1.8,'S':1.8}

   NS = pdb.NeighborSearch(list(structure.get_atoms()))

   Number_clashes = 0

   for at in rotating_chain.get_atoms():
      neighbors = NS.search(at.get_coord(), 2.0)



      # if you find a close atom, ask if it is a clash
      if len(neighbors)>0:

         for neigh in neighbors:

            sum_VDW_radius = Van_der_Waals_radius[neigh.element] + Van_der_Waals_radius[at.element]
            Distance = neigh - at

            if Distance < sum_VDW_radius:
               Number_clashes += 1
               print(at,neigh)
               print(sum_VDW_radius)
               print(Distance)
               print(Number_clashes)

            if Number_clashes==100:
               print(Number_clashes)
               return True

   return False



# check:
parser = pdb.PDBParser(PERMISSIVE=1)
structure = parser.get_structure('pdb_name', '3kuy.pdb')
structure1 = structure[0]['B']
rotating_chain = structure[0]['B']

a = is_Steric_clash(structure1,rotating_chain)
print(a)
if a is True:
   print(a)
   print ('theres a clash')

