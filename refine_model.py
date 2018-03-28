# this script inputs a pdb file, refines it and outputs the refined model

import argparse

from pyrosetta.toolbox import *
import pyrosetta
pyrosetta.init()

from pyrosetta.toolbox import cleanATOM
from pyrosetta.teaching import *

from contextlib import redirect_stdout

# arguments

parser = argparse.ArgumentParser(description="This script optimizes a pdb and prints information aboy the energy calculation of the optimized structure") 

parser.add_argument('-o', '--output-file',
   dest="outputfile",
   action="store",
   default="./refined_model.pdb",
   help="File name of the optimized model")

parser.add_argument('-i', '--input',
   dest="input",
   action="store",
   default=None,
   help="Path of the PDB file to optimize. It has to have .pdb extension")

options = parser.parse_args()

# debug the input
if options.input is None:
    raise Exception("No input provided") 	

if options.input.split('.')[-1] != 'pdb':
	raise Exception("The input has to be a .pdb file")

# clean the input file to be correct for pyrosetta
cleanATOM(options.input)
clean_path = options.input.split('.')[0]+'.clean.pdb'

# initialize the pose object that will be optimized
pose = pyrosetta.pose_from_pdb(clean_path)

# define a general energy score function for the energy optimization
scorefxn = get_fa_scorefxn()

# apply a classic relax pipeline(from P. Bradley, K. M. S. Misura & D. Baker, “Toward high-resolution de novo structure prediction for small proteins.” Science 309, 1868-1871 (2005)) 
relax = pyrosetta.rosetta.protocols.relax.ClassicRelax()
relax.set_scorefxn(scorefxn)
relax.apply(pose)

# save into a pdb file
pose.dump_pdb(options.outputfile)

# save the info about the energy of the refined model

print("This file contains information abot the energy components of the refined pdb structure %s.\n\n"%(options.outputfile))
print("Below is the calculation of the energy scoring function terms:\n\n")

# energy scoring function terms
scorefxn.show(pose)



