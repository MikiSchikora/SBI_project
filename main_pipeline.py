# This is the main script for modelling a complex from a set of Pairwise interactions

import os
import sys
sys.path.append(os.curdir)
import Functions as func
import Bio.PDB as pdb
import argparse
import random
import copy


parser = argparse.ArgumentParser(description="This program builds a complex from the interacting pairwise subunits")  # WRITE DESCRIPTION!!!!!

parser.add_argument('-o', '--output-dir',
   dest="outputdir",
   action="store",
   default="./output_models/",
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


parser.add_argument('-sto', '--stoichiometry',
   dest="stoich",
   action="store",
   default=None,
   help="Add a tabular file with sequence ID and stoichiometry")

parser.add_argument('-n_models', '--number-of-models',
   dest="n_models",
   type=int,
   action="store",
   default=1,
   help="This pipeline can generate several models out of one input. With this you can indicate the number of models (an integrer) that have to be generated")

parser.add_argument('-n_chains', '--number-of-chains',
   dest="n_chains",
   type=int,
   action="store",
   default=1000,
   help="Specify the maximum number of chains that a model should have")


options = parser.parse_args()

# Parse the input arguments:
if options.input:

    # process the input
    if os.path.isfile(options.input):

        # when the input is a file you have to generate all the interacting pairs, rotated and translated
        if options.input.split('.')[-1] == 'pdb' or options.input.split('.')[-1] == 'ent':
            filetype = 'PDB'
        elif options.input.split('.')[-1] == 'cif':
            filetype = 'CIF'
        else:
            raise Exception('The provided complex has to be a PDB or mmCIF file')

        if options.verbose:
            print("Your input is a complex already. This is the program testing mode.\n Splitting input into pairwise subunits...")

        templates_dir = './templates/'
        func.generate_pairwise_subunits_from_pdb(options.input, templates_dir, filetype, options.verbose)

    elif os.path.isdir(options.input):
        templates_dir = options.input

    else:
        raise Exception('You have to provide a valid input path')
else:
    raise Exception("No input provided")

if options.exhaustive:
    sys.stderr.write("!!!!!!!!! WARNING: you specified the exhaustive option.  !!!!!!!! \n"
          "This means that the program will undergo a tree-like recursive approach for building any possible complex with your input files, taking a long time for large structures.\n"
          "We recommend not to use this unless you did not succeed with any other options. \n"
          "If you suspect that different isoforms could result of your input we recommend to preset the number of isoforms expected (with the option n_models) and/or specify the desired stoichiometry (with option -sto).\n")


# define the output, if it is not specified, then it is the default
output = options.outputdir
if not os.path.exists(output):
    os.makedirs(output)

# handle multifasta file input
# a file containing the sequences of the subunits to include. This is mandatory for any chain that is not a random DNA seuqnce.
if options.sequences:
    # if a multifasta is provided
    subunits_seq_file = options.sequences
else:
    # if a multifasta is not provided
    subunits_seq_file = None
    sys.stderr.write("!!!!!!!!! WARNING: you did not specify the expected molecules of the complex. !!!!!!!! \n"
          " We recommend to provide a file with the expected molecules that are part of the complex with the -seq option.\n "
          "This is important for you to understand which are the molecules that form each chain in the resulting complex. If not specified, each molecule gets an arbitrary identifier.\n ")

if options.stoich:
    if options.sequences:
        stoich_file = options.stoich
        fd = open(stoich_file)
        stoich_dict = {}
        for line in fd:
            line = line.strip()
            try:
                seq_id, sto = line.split("\t")
            except ValueError:
                raise Exception("Your stoichiometry file does not have the correct format")

            stoich_dict[seq_id] = int(sto)
    else:
        raise Exception("You have to provide a MULTIFASTA to link ID and sequence if you want stoichiometry to be considered.")
else:
    stoich_dict = None
    sys.stderr.write("!!!!!!!!! WARNING: you did not specify the expected stoichiometry of the complex. !!!!!!!! \n"
          " We recommend to provide a file with the expected stoichiometry with the -sto option. This is important for filtering out non-expected complexes. Use this option only if you have prior evidence about the stoichiometry and you expect many possible stoichiometries.\n")

number_subunits_file = './subunits_num.tbl'  # a file containing the number of subunits, if known


# generate info about the templates
if options.verbose:
    print("Generating information about your input files ... ")

# PDB_info information about the unique chains: Keys: filename, Values: [(chain accession , molecule identifier) ... ]
# Seqs_info has the sequence id and sequence

PDB_info, Seqs_info = func.generate_PDB_info(templates_dir, subunits_seq_file, options.verbose)


# initialise PDB files parser
p = pdb.PDBParser(PERMISSIVE=1)

# Pick initial filename randomly:
filename = random.choice(list(PDB_info.keys()))

# we start with the structure of the first pairwise interaction, this is now the current model
current_structure = p.get_structure("pr1", templates_dir + filename)

if len(list(current_structure.get_chains())) != 2:
    raise func.TwoChainException(filename)


# change the chain id names
for chain in current_structure.get_chains():
    curr_id = chain.id  # the so-called chain accession, usually a letter (A,B,C,D ...)
    chain.id = [x for x in PDB_info[filename] if x[0] == curr_id][0]

# NOTE: all the chains, but the initial 2, in the current_structure will have this naming for a proper working of the code:
    # (A,H2A3_B,AGS6G,3), that corresponds to (chain_accession, chain_id, random_id, recursivity level at which it has been added)

if options.verbose:
    print(" BUILDING THE COMPLEX...")

# initialize vars
final_models = []

# then we call the function 'build_complex'
final_models = func.build_complex(final_models, current_structure, templates_dir, PDB_info, options.n_models, options.exhaustive, options.n_chains, stoich=stoich_dict, verbose=options.verbose)

if len(final_models) == 0:
    if options.stoich:
        sys.stderr.write("No complex could be obtained with your specified stoichiometry. We recommend specifying a number of models to explore (-n_models) or run the default pipeline (as sometimes the obtained stoichiometry is very similar to the expected). \n")
    else:
        sys.stderr.write("No complex could be obtained with the provided files. \n")

# printing the model
for final_model in final_models:

    # if final model has more than 62 chains, split into different models
    structure = func.pdb_struct.Structure('id')  # empty structure:

    # create a new model for each group of 62 chains, and add to structure and change the id of the chains and record the info about which is the molecule of each chain
    chain_counter = 0
    model_counter = 0
    legend = ""

    # change the ids of the chains to completely random, for further
    chains = list(final_model.get_chains())

    for chain in chains:

        # define the new id:
        new_id = func.chain_alphabet[chain_counter]

        # update what has to be printed:
        legend += 'CHAIN    %s   %s   %s\n' % (new_id, chain.id[1], Seqs_info[chain.id[1]])

        # change to the new id
        chain.id = new_id

        # improve the chain counter
        chain_counter += 1

        if options.verbose:
            print('printing chain: ', chain_counter, chain.id)

        # initialize model
        if chain_counter == 1:
            model_counter += 1
            model = func.pdb_model.Model(model_counter)

        # add chain to the model
        model.add(chain)

        # reset counter and add model after 62 rounds
        if chain_counter == len(func.chain_alphabet):
            chain_counter = 0

            model = copy.deepcopy(model)
            structure.add(model)

            # remove the chains in the current model of  of the final_model, to avoid further problems with the id
            for chain_m in model.get_chains():
                final_model[0].detach_child(chain_m.id)

    # add the last model:
    current_models_strcuture = [x.id for x in structure.get_models()]
    if model.id not in current_models_strcuture:
        structure.add(model)

    # define the name of the model, not to overwrite previous existing files
    written_id = 0
    PDB_name = 'model_' + str(written_id) + '.pdb'
    while PDB_name in os.listdir(output):
        written_id += 1
        PDB_name = 'model_'+str(written_id)+'.pdb'

    if options.verbose:
        print('%s created as a model for the complex structure.' % PDB_name)

    # open the pdb and put information about what is each chain
    model_path = output + PDB_name

    # Write the ATOM lines:
    io = pdb.PDBIO()  # using our own PDB writer
    io.set_structure(structure)
    io.save(model_path)

    # write the info lines
    fd = open(model_path, 'a')
    fd.write('CHAIN HEADER    current id   molecule id    sequence\n')
    fd.write(legend)
    fd.write('\n')
    fd.close()
