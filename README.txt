
######### INTRODUCTION #########

This is a usage tutorial of the BuiltYourComplex command-line application. It aims to build a multi-molecular complex from the pairwise interacting subunits, provided by the user.

Developers: Miquel Àngel Schikora Tamarit and Marina Lleal Custey, from the UPF MSc of Bioinformatics for Health Sciences.



######### SYSTEM REQUIREMENTS #########

This application was written and debuged using python 3.5 and 3.6 versions, so we recommend using any of these for a proper usage.

Other application requirements (python modules and packages):

- Biopython package (version 1.7 was the used by the developers, so that any posterior release should also work), with the modules Bio.PDB, Bio.SeqIO, Bio.pairwise2, Bio._py3k, Bio.Data

- PyRosetta4 package (release 173 or later, we have tried it for python 3.6)

- python numpy 

- python os 

- python copy 

- python random 

- python rmsd

- python string

- python argparse

- python sys

######### INSTALLATION #########



######### USAGE #########

The BuildYourComplex may use two scripts for building the complex:

- main_pipeline.py builds the complex from a provided directory of pairwise interacting subunits. It can also split a provided pdb into pairwise interacting subunits and rebuild it (testing mode of the pipeline).


	ARGUMENTS:

	-i specifies the input directory or file. If -i is a directory the program will assume that the contents are pdb files with each of the pairwise interactions. If it is a file it will split into pairwise interacting chain, and create a directory, TEMPLATES/, that will store each of them. This argument is obligatory.

	-o specifies the directory where the created models have to be stored. If the directory is not empty the new models created will be added to this. The default is output_models/

	-v indicates if the progression log has to be printed. You can redirect it to a file with: 
		python main_pipeline.py -i INPUTDIR/ -v > log.txt

	-seq comes together with a multifasta file with the expected molecule (protein, RNA or DNA) sequences in the complex. Provoding this file is necessary for a proper identification of the molecules in the input PDB, and their corresponding identification in the output models. When absent, each of the molecules (often shared between several input files) will get a unique random identifier. An example of such a file for the nucleosome complex is in example_sequences.fa.

	-exh indicates that all the possibilities have to be explored. This pipeline can assume that the provided input files may generate different complexes, so it implements a branch-opening recursive function that can generate many different complexes. The -exh argument forces to explore all the possible paths of this "building-tree", and it may be very time-consuming for large complexes. It is not recommended.

	-sto comes together with a tab-delimited file (see example_stoichiometry.tbl for an example of the nucleosome complex) that indicates which is the expected  stoichiometry of each of the molecules. It has to be provided together with the -seq option. We recommend to use this argument if there's previous evidence indicateing that the complex may have different stoichiometry, and you are interested in a specific one. When provided, the pipeline also implements a branch-opening recursive function for ensuring that the expected stoichiometry will be found if possible.

	-n_models comes together with an integrer that indicates which is the maximum number of models (with different structure) you want to create (default is 1). We recommend to use it only when there's evidence that many possible conformations can be found, as this option also implements branch-opening recursive function that can be slow for  large complexes.

    -n_chains comes together with the maximum number of chains that the model will have. Default is 1000. Use it if you expect you complex to be a

	GENERAL COMMENTS:

	Using the -n_models, -exh and/or -sto arguments will tell our program to assume that many possible complexes can be built with the provided input. This calls a branch-opening recursive function that can be very slow, because it has to try all the possibilities, so be careful with using these options. -sto and -n_models are better than -exh, because the execution will stop after having generated the expected number of models and/or found the expected stoichiometry, so it will be faster.

	The general recommendation is not to use any of the -n_models or -exh options if you obtain models, so simply call this script like:

	python3 main_pipeline.py -i INPUT_DIR/ -o OUTPUT_DIR/


	EXAMPLES:

	We have used the "testing mode" of our pipeline to rebuild several complexes, all of them stored in example_inputs/. This folder contains, for each complex, the original pdb file and a folder with the split pairwise subunits. We have built all these with the default call (python3 main_pipeline.py -i INPUT_DIR/ -o OUTPUT_DIR/), because we expect only one solution for each. They include:

	- Nucleosome

	- Phosphate dehydratase

	- Proteasome

	- ATP synthase

	- Ribosome



	You can use any of them (input or ouput files) for testing the functionality of the program, and try different options.

	The output models are stored in example_outputs/ with the corresponding naming.


	WHAT IS PRINTED TO THE SCREEN DURING THE EXECUTION PROGRESSION?

	When the -v argument is indicated, this script prints to the standard output the progression of the execution, that includes:

		- A recall of the input files that are being written if the input is a pdb file. 

		- A recall of the moment in which the input files' information is being processed and the start of building the complex.

		- A recall of which input file is being processed at the corresponding time.

		- Each time that build_complex(), the function that adds new chains to the complex, is called it indicates which is the current recursivity level of two distinct features:

			Branch level indicates the node in the branch-building tree in which the current build_complex() is found. Higher branch levels indicate that the complex-building process is immersed in deep branches of the tree-building (see the theroretical background for more information).

			Complex building level indicates how many times build_complex() has been called recursively for adding new chains to a previously existing feature. 

		- Each time that build_complex() the composition of the complex (how many chains of each unique molecule are present in the current model) is recorded

		- Each time a chain of an input file is trying to be added to the complex the unique identifier of this chain (a tuple containing the original chain ID, the unique molecule ID, a random identifier and the Complex building recursion level at which it was added) is printed.

		- The reason for a branch-opening event (chains clashing against the complex) is always reported while opening.

	    - While saving the model created the pipeline prints which is the chain molecule id and the random id.

	If you are interested in a deeper understanding of this information we suggest you to read carefully the code and the theoretical explanation. 


- refine_model.py runs an optimization process (using the pyrosetta4 package) on a given structure (it should be one of the models generated by main_pipeline.py), to minimize the energy of the final model and turn it into a "relaxed state".

	ARGUMENTS:

	-i is mandatory and indicates the input .pdb file, that should be one of the models generated with main_pipeline-py

	-o is the file-path of the refined model (it will have been optimized with rosetta)

	GENERAL COMMENTS:

	Our pipeline runs a classical relax optimization pipeline (described in P. Bradley, K. M. S. Misura & D. Baker, “Toward high-resolution de novo structure prediction for small proteins.” Science 309, 1868-1871 (2005)) on the provided structure. 

	This may take ~40 min for a structure of 100 residues, so you should be aware of using it for large complexes.

	The progression prints information to the terminal about the refinement and energy calculation process of the refined model. Here we describe which is the meaning of each of the energy terms reported at the end of running this script (they refer to the energy scoring functions of the refined structure):

		fa_atr                                     Lennard-Jones attractive between atoms in different residues.  Supports canonical and noncanonical residue types.
		fa_rep                                     Lennard-Jones repulsive between atoms in different residues.  Supports canonical and noncanonical residue types.
		fa_sol                                     Lazaridis-Karplus solvation energy.  Supports canonical and noncanonical residue types.
		fa_intra_rep                               Lennard-Jones repulsive between atoms in the same residue.  Supports canonical and noncanonical residue types.
		fa_intra_sol_xover4                        Intra-residue LK solvation, counted for the atom-pairs beyond torsion-relationship.  Supports arbitrary residues types.                 
		lk_ball_wtd                                Weighted sum of lk_ball & lk_ball_iso (w1*lk_ball + w2*lk_ball_iso); w2 is negative so that anisotropic contribution(lk_ball) replaces some portion of isotropic contribution (fa_sol=lk_ball_iso).  Supports arbitrary residue types.
		fa_elec                                    Coulombic electrostatic potential with a distance-dependent dielectric.  Supports canonical and noncanonical residue types.
		pro_close                                  Proline ring closure energy and energy of psi angle of preceding residue.  Supports D- or L-proline, plus D- or L-oligourea-proline.
		hbond_sr_bb                                Backbone-backbone hbonds close in primary sequence.  All hydrogen bonding terms support canonical and noncanonical types.
		hbond_lr_bb                                Backbone-backbone hbonds distant in primary sequence.
		hbond_bb_sc                                Sidechain-backbone hydrogen bond energy.
		hbond_sc                                   Sidechain-sidechain hydrogen bond energy.
		dslf_fa13                                  Disulfide geometry potential.  Supports D- and L-cysteine disulfides, plus homocysteine disulfides or disulfides involving beta-3-cysteine.
		rama                                       Ramachandran preferences.  Supports only the 20 canonical alpha-amino acids and their mirror images.
		omega                                      Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.  Supports alpha-amino acids, beta-amino acids, and oligoureas.  In the case of oligoureas, both amide bonds (called 'mu' and 'omega' in Rosetta) are constarined to planarity.
		fa_dun                                     Internal energy of sidechain rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).  Supports any residue type for which a rotamer library is avalable.
		p_aa_pp                                    Probability of amino acid at Φ/Ψ.  Supports only the 20 canonical alpha-amino acids and their mirror images.
		yhh_planarity                              Sidechain hydroxyl group torsion preference for Tyr residues.
		ref                                        Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.  Supports only the 20 canonical alpha-amino acids and their mirror images.
		rama_prepro                                Backbone torsion preference term that takes into account of whether preceding amono acid is Proline or not.  Currently supports the 20 canonical alpha-amino acids, their mirror-image D-amino acids, oligoureas, and N-methyl amino acids.  Arbitrary new building-blocks can also be supported provided that an N-dimensional mainchain potential can be generated somehow.



	You may consider to redirect the progression log (for future usage) like this:

		python3 refine_model.py -i input_model.pdb -o refined_model.pdb  > log.txt

