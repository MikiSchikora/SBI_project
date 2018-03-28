
#this is an example for refining the structure of the
from pyrosetta.toolbox import *

import pyrosetta

pyrosetta.init()

from pyrosetta.toolbox import cleanATOM

cleanATOM("./example_pair.pdb")

print(dir(pyrosetta.rosetta))

print(dir(pyrosetta))

pose = pyrosetta.pose_from_pdb('example_pair.clean.pdb')

from pyrosetta.teaching import *

scorefxn = get_fa_scorefxn()

print(scorefxn)

print(pose.phi(5))

relax = pyrosetta.rosetta.protocols.relax.ClassicRelax()
relax.set_scorefxn(scorefxn)
relax.apply(pose)

pose.dump_pdb('example_pair_relaxed.pdb')

relax.set_scorefxn(scorefxn)
relax.apply(pose)

print(pose)

pose = pose_from_pdb("small.clean.pdb")

print(pose)

pose = pose_from_pdb("1YY8.clean.pdb")





# stuff associated to the energy calculation

import pyrosetta
pyrosetta.init()

from pyrosetta.toolbox import cleanATOM

cleanATOM('example_pair_relaxed.pdb')
pose = pyrosetta.pose_from_pdb('example_pair_relaxed.clean.pdb')

from pyrosetta.teaching import *
scorefxn = get_fa_scorefxn()

#print(scorefxn(pose))

print(scorefxn.show(pose))




print(scorefxn(pose))