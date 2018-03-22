# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import Functions as func
import Bio.PDB as pdb
import argparse


parser = argparse.ArgumentParser(description="This program does BLA BLA BLA")  # WRITE DESCRIPTION!!!!!

parser.add_argument('-o', '--output-dir',
   dest="outputdir",
   action="store",
   default="./output_models",
   help="output directory")


parser.add_argument('-i', '--input',
   dest="input",
   action="store",
   default=None,
   help="Input a directory with several PDB files or a PDB complex to test the program")


parser.add_argument('-v', '--verbose',
   dest="verbose",
   action="store_true",
   default=False,
   help="Print log in stderr")


parser.add_argument('-seq', '--sequences',
   dest="sequences",
   action="store",
   default=None,
   help="Add a MULTIFASTA with the sequences of the different subunits of the complex to be modelled")


parser.add_argument('-exh', '--exhaustive',
   dest="exhaustive",
   action="store_true",
   default=False,
   help="Try all the possibilities. It may take a long time!")


options = parser.parse_args()

# Parse the input arguments:
if options.input:
    # process the input
    if os.path.isfile(options.input):

        # when the input is a file you have to generate all the interacting pairs, rotated and translated
        if options.input.split('.')[-1] == 'pdb':  # OTHER EXTENSIONS??????????
            filetype = 'PDB'
        elif options.input.split('.')[-1] == 'cif':
            filetype = 'CIF'
        else:
            raise EnvironmentError('The provided complex has to be a PDB or mmCIF file')

        Templates_dir = './TEMPLATES/'
        func.Generate_pairwise_subunits_from_pdb(options.input, Templates_dir, filetype)

    elif os.path.isdir(options.input):
        Templates_dir = options.input

    else:
        raise EnvironmentError('You have to provide a valid input path')
else:
    raise IOError("No input provided")


# define the output, if it is not specified, then it is the default
output=options.outputdir

# handle multifasta file input
# a file containing the sequences of the subunits to include. This is mandatory for any chain that is not a random DNA seuqnce.
if options.sequences:
    # if a multifasta is provided
    subunits_seq_file = options.sequences
else:
    # if a multifasta is not provided
    subunits_seq_file = None
    print("NO SUBUNITS SEQ FILE")


stoichiometry_file = './stoch.tbl'  # a file containing information about the stoichiometry. This is mandatory
number_subunits_file = './subunits_num.tbl'  # a file containing the number of subunits, if known


# generate info about the Templates
# PDB_info information about the unique chains: Keys: filename, Values: {Chain: unique ID}
# Uniq_seqs is a set with all the unique IDs
PDB_info, Seq_to_filenames = func.Generate_PDB_info(Templates_dir, subunits_seq_file)

# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

final_models = []

# parse file 1 and 2 and get a structure for each
for filename1 in os.listdir(Templates_dir):

    # we start with the structure of the first pairwise interaction, this is now the current model
    current_structure = p.get_structure("pr1", Templates_dir + filename1)

    # change the chain id names
    for chain in current_structure.get_chains():
        curr_id=chain.id
        chain.id=[x for x in PDB_info[filename1] if x[0]==curr_id][0]

    # NOTE: all the chains, but the initial 2, in the current_structure will have this naming for a proper working of the code:
        # (A|||H2A3_B|||AGS6G), that corresponds to (chain_accession|||chain_id|||random_id)

    # initialize some global vars

    # tried_operations = set() # a set containing the tried operations (chain1_full_id, filename2, rotating_chain_accession)
    rec_level = 0  # the recursivity level

    # then we call the function 'build_complex'
    final_models = func.build_complex(final_models, current_structure, Templates_dir, PDB_info, Seq_to_filenames)
    break


for final_model in final_models:
    # we go through all the chains of the structure and rename them alphabetically
    final_model.id = '_'.join([str(x.id) for x in final_model.get_chains()])

    print(final_model.id)

    final_chains = list(final_model.get_chains())  # list of chain objects
    chain_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S",
                      "T", "U", "V", "W", "X", "Y", "Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o",
                      "p","q","r","s","t","u","v","w","x","y","z","1","2","3","4","5","6","7","8","9","0"]
    alphabet_pos = 0
    for final_chain in final_chains:
        final_chain.id = chain_alphabet[alphabet_pos]
        alphabet_pos += 1

    # then we can finally save the obtained structure object into a pdb file
    written_id = 0
    PDB_name = 'model_' + str(written_id) + '.pdb'
    while PDB_name in os.listdir('./Output_models/'):
        written_id += 1
        PDB_name = 'model_'+str(written_id)+'.pdb'


    io = pdb.PDBIO()
    io.set_structure(final_model)
    io.save('./Output_models/'+PDB_name)
