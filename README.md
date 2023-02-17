# zeoran

ZEOlite RANdom generation

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

For a detailed description of the algorithms used to generate the zeolite frameworks and the program itself, please check the .pdf file available in this repository. The software is entirely written in C/C++ and contains a Shell script to execute it. 


## Content
Here the content of the repository is described. 

(1) zeoran_data/atom_sites: Contains the atomic positions of all the atoms forming a single unit cell of each zeolite. The structure of a line of this file is: atom_id atom_type x y z charge.

(2) zeoran_data/unit_cell: Contains the basic information of the unit cell of each framework.

(3) makefile: script to compile and install the program.

(4) functions.cpp: contains all the functions implemented to run the program, including the main function.

(5) generate.input: input file.

(6) zeoran.pdf: detailed explanation of the program.

(7) global.h: contains the definition of the global variables and structures used in the program.

(8) headers.h: contains the headers of all the functions written in the file "functions.cpp".

(9) libraries.h: contains all the libraries used.

(10) log.dat: contains the current state of the program and changes and improvements that still need to be done.


## Run the program


#### Requirements
	
	1. The compiler `g++` and `make` need to be installed. (`sudo apt-get install build-essential`)
	2. The standard Eigen library needs to be installed. (`sudo apt-get install libeigen3-dev`)
	
#### Compilation and installation

The program can be compiled and installed in the system using the makefile. Open a terminal in the directory containing all the files and introduce the commands:

    	```
	make
	sudo make install
    	```

After installing, the source folder may be removed from your user space.

#### Execution

    	```
	zeoran
    	```

#### Uninstalling

    	```
        sudo make uninstall
    	```

make uninstallation is only available if you did not delete the source folder.
If you did remove it, you can uninstall with `sudo rm -rf /usr/local/bin/zeoran /user/local/share/zeoran`

## Structure of generate.input
The generate.input file needs to have a particular structure. Please, note that he number of parameters and their order cannot be changed. This structure has a fixed part, and a part that depends on the algorithm we aim to use (chains/clusters/merw/random), since they require different parameters. Next, we specify the format of the input file in each case:

#### chains:

	Zeolite (MFI/MOR/FAU/RHO/MEL/TON/new)
	Algorithm (chains/clusters/merw/random)
	Name_of_output_directory
	Number_of_structures_to_be_generated
	Number_of_chains
	Number_of_substitutions_in_first_chain
	Number_of_substitutions_in_second_chain
	...
	Number_of_substitutions_in_last_chain

#### clusters:

	Zeolite (MFI/MOR/FAU/RHO/MEL/TON/new)
	Algorithm (chains/clusters/merw/random)
	Name_of_output_directory
	Number_of_structures_to_be_generated
	Number_of_substitutions
	
#### merw:

	Zeolite (MFI/MOR/FAU/RHO/MEL/TON/new)
	Algorithm (chains/clusters/merw/random)
	Name_of_output_directory
	Number_of_structures_to_be_generated
	Number_of_substitutions
	Number_of_equilibration_steps
	Number_of_visits

NOTE: The last two parameters are optional. If they are missing the default values are 100 and 20 respectively. Note that if changed, both of them need to be specified. More information about these parameters can be found in the file generate_zeolite_frameworks.pdf.

#### random:

	Zeolite (MFI/MOR/FAU/RHO/MEL/TON/new)
	Algorithm (chains/clusters/merw/random)
	Name_of_output_directory
	Number_of_structures_to_be_generated
	Number_of_substitutions

## How to work with a new zeolite

The software contains the information to generate structures for the MFI, MOR, FAU, RHO, TON and MEL frameworks. However, it is implemented to work with all the zeolites. Therefore, to generate structures of other frameworks, only the corresponding files need to be added in the directories "data/atom_sites" and "data/unit_cell". Also, the new zeolite needs to be specified in the first line of generate.input.
