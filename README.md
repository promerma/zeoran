# Generate Zeolite Frameworks

## Author's information
	Pablo Romero Marimon
	Eindhoven University of Technology
	Department of Applied Physics
	Materials Simulation and Modelling group
	February 10, 2023

## Description
This repository was done as part of my Master's thesis at the MS&M group. Its main goal is to generate multiple .cif files with zeolite structures having different aluminum substitutions. The aluminum substitutions are randomly generated from 4 different probability distributions, implemented by 4 different algorithms:

(1) chains: Al atoms are introduced forming chains of consecutive Al-O-Al bonds. The number and lengths of the chains are given as an input.

(2) clusters: A given number of Al atoms are introduced in a small spatial region.

(3) merw: A given number of Al atoms are introduced "as spread as possibly" in the structure, i.e., maximizing the entropy of the framework.

(4) random: A given number of Al atoms are introduced by sampling a uniform distribution.

For a detailed description of the algorithms used to generate the zeolite frameworks and the program itself, please check the .pdf file available in this repository. The software is entirely written in C/C++ and contains a bash script to execute it. 

## Run the program
The repository contains an executable shell script that compiles and executes the program. To execute the program open the terminal in the directory containing this file and type:

	./execute
	
The first time, it might be necessary to give execution permission to the file execute. This can be done by right-click on properties or through the terminal by typing:

	chmod +x execute

#### Requirements:
	
	1. The g++ needs to be installed.
	2. The data directory needs to be in the current directory.
	3. The files functions.cpp, generate.input, global.h, heaaders.h and libraries.h need to be in the current directory.

## Content
Here the content of the repository is described. 

(1) data/atom_sites: Contains the atomic positions of all the atoms forming a single unit cell of each zeolite. The structure of a line of this file is: atom_id atom_type x y z charge.

(2) data/unit_cell: Contains the basic information of the unit cell of each framework.

(3) execute: shell script to run the program.

(4) functions.cpp: contains all the functions implemented to run the program, including the main function.

(5) generate.input: input file.

(6) generate_zeolite_frameworks.pdf: detailed explanation of the program.

(7) global.h: contains the definition of the global variables and structures used in the program.

(8) headers.h: contains the headers of all the functions written in the file "functions.cpp".

(9) libraries.h: contains all the libraries used.

(10) log.dat: contains the current state of the program and changes and improvements that still need to be done.


## Structure of generate.input
The generate.input file needs to have a particular structure. Please, note that he number of parameters and their order cannot be changed. This strcture is as follows:

	Zeolite (MFI/MOR/FAU/RHO/MEL/TON/new)
	Algorithm (chains/clusters/merw/random)
	Name_of_output_directory
	Number_of_structures_to_be_generated
	Extra_parameters

The amount and type of Extra_parameters depend on the algorithm we are using (chains/clusters/merw/random). In each case, these parameters are:

#### chains:

	Number_of_chains
	Number_of_substitutions_in_first_chain
	Number_of_substitutions_in_second_chain
	...
	Number_of_substitutions_in_last_chain

#### clusters:

	Number_of_substitutions
	
#### merw:
	
	Number_of_substitutions
	Number_of_equilibration_steps
	Number_of_visits

NOTE: The last two parameters are optional. If they are missing the default values are 100 and 20 respectively. Note that if changed, both of them need to be specified.

#### random:

	Number_of_substitutions

## How to work with a new zeolite

The software contains the information to generate structures for the MFI, MOR, FAU, RHO, TON and MEL frameworks. However, it is implemented to work with all the zeolites. Therefore, to generate structures of other frameworks, only the corresponding files need to be added in the directories "data/atom_sites" and "data/unit_cell". Also, the new zeolite needs to be specified in the first line of generate.input.
