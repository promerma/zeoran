#include "libraries.h"
#include "global.h"
#include "headers.h"


int main (void) {

	string name_zeo, name_alg, file_zeo, file_ucell, out_name;
	int Nstruct, Nchains, Nals;
	vector<int> chains, Als;
	vector<vector<int> > neigbrs;
	atom *list;
	ifstream fin;

	//Use current time as seed for random generator
    srand(time(NULL));

	//Check if input exists and if so read
	fin.open("generate.input");
    if(fin.fail()) {
        cerr << "unable to open file generate.input for reading" << endl;
        exit( 1 );
    }
	fin >> name_zeo >> name_alg >> out_name >> Nstruct;


    // Structure which would store the metadata
    struct stat sb;

    //Check if output directory already exists
    if (stat(out_name.c_str(), &sb) == 0) {
    	string yn;
        cout << "Output directory already exists. Do you want to continue anyway (this might overwrite existing files)? (yes/no)" << endl;
    	cin >> yn;
    	if(yn.compare("no") == 0) {
    		cout << "The program has been cancelled." << endl;
    		exit(1);
    	}

    //If directory does not exist, create it
    } else {
    	mkdir(out_name.c_str(), 0777);
    }


	if(name_alg.compare("chains") == 0) {
		fin >> Nchains;
		int aux;
		for(int i=0; i<Nchains; i++) {
			fin >> aux;
			chains.push_back(aux);
		}

	} else if(name_alg.compare("clusters") == 0) {
		fin >> Nals;

	} else if (name_alg.compare("merw") == 0) {
		fin >> Nals;
		if (!fin.eof()) {
			fin >> Neqsteps;
			fin >> Nvisits;
		}

	} else if (name_alg.compare("random") == 0) {
		fin >> Nals;

	} else {
		cerr << "Wrong algorithm name introduced." << endl;
		exit(1);
	}

	fin.close();

	//Read zeolite and unit cell
	file_zeo = "/usr/local/share/zeoran/atom_sites/" + name_zeo + ".txt";
	file_ucell = "/usr/local/share/zeoran/unit_cell/" + name_zeo + ".txt";

	read_unit_cell(file_ucell);
	list = read_atom_sites(file_zeo);

	//Select the algorithm
	if(name_alg.compare("chains") == 0) {

		int **M_T, *Tids;

		//Generate all the structures
		for(int i=0; i<Nstruct; i++) {

			//Output current state
			cout << "Progress: " << i+1 << "/" << Nstruct << endl; 

			//Save indexs of T atoms in list
			Tids=(int*)malloc(Tatoms*sizeof(int));

			//Generate adjacency matrix
			M_T=get_amatrix_Tatoms(list, Tids);

			//For each atom find its neighbours
			neigbrs=get_neighbours(M_T);

			//Generate Al chains
			Als=generate_chains(M_T, neigbrs, chains);

			//Print structure
			print_structure(list, Als, i, name_zeo, name_alg, out_name);

			free(Tids);
		}		

	} else if(name_alg.compare("clusters") == 0) {

		int **M_T, *Tids;

		//Generate all the structures
		for(int i=0; i<Nstruct; i++) {

			//Output current state
			cout << "Progress: " << i+1 << "/" << Nstruct << endl; 			

			//Save indexs of T atoms in list
			Tids=(int*)malloc(Tatoms*sizeof(int));

			//Generate adjacency matrix
			M_T=get_amatrix_Tatoms(list, Tids);

			//For each atom find its neighbours
			neigbrs=get_neighbours(M_T);

			//Generate Al chains
			Als=clusters_substitutions(M_T, Nals);

			//Print structure
			print_structure(list, Als, i, name_zeo, name_alg, out_name);

			free(Tids);
		}	

	} else if (name_alg.compare("merw") == 0) {

		int **M_T, *Tids;

		//Generate all the structures
		for(int i=0; i<Nstruct; i++) {

			//Output current state
			cout << "Progress: " << i+1 << "/" << Nstruct << endl; 

			//Save indexs of T atoms in list
			Tids=(int*)malloc(Tatoms*sizeof(int));

			//Generate adjacency matrix
			M_T=get_amatrix_Tatoms(list, Tids);

			//For each atom find its neighbours
			neigbrs=get_neighbours(M_T);

			//Generate Al chains
			Als=merw_substitutions(M_T, Nals);

			//Print structure
			print_structure(list, Als, i, name_zeo, name_alg, out_name);

			free(Tids);
		}			

	} else if (name_alg.compare("random") == 0) {

		//Generate all the structures
		for(int i=0; i<Nstruct; i++) {

			//Output current state
			cout << "Progress: " << i+1 << "/" << Nstruct << endl; 			

			vector<int> Als = generate_random(Nals);
			
			//Generate structure
			print_structure(list, Als, i, name_zeo, name_alg, out_name);
		}
	}


	return 0;
}



atom* read_atom_sites (string fname) {

	int i, aux;
	atom *list;
	ifstream fin;

	list=(atom*)malloc(Natoms*sizeof(atom));

	fin.open(fname.c_str());

	for(i=0; i<Natoms; i++) {
		fin >> list[i].at >> list[i].id >> list[i].x >> list[i].y >> list[i].z >>list[i].q;
	}

	fin.close();

	return list;
}

void read_unit_cell (string fname) {

	string aux;
	ifstream fin;
	fin.open(fname.c_str());

	fin >> aux >> aux >> aux >> Natoms;
	fin >> aux >> aux >> aux >> Tatoms;
	fin >> aux >> a >> aux >> b >> aux >> c;
	fin >> aux >> alpha >> aux >> beta >> aux >> gama;
	fin >> aux >> setting;
	
	return;
}


int connectivity (atom a1, atom a2) {

	double d, dx, dy, dz, dpbc;

	//Compute distance in all directions
	dx=fabs(a1.x-a2.x);
	dy=fabs(a1.y-a2.y);
	dz=fabs(a1.z-a2.z);

	//Check boundaries
	if(a1.x <= a2.x) {
		dpbc=fabs(a1.x + 1 - a2.x);
	} else {
		dpbc=fabs(a2.x + 1 - a1.x);
	}
	if(dpbc < dx) {
		dx=dpbc;
	}

	if(a1.y <= a2.y) {
		dpbc=fabs(a1.y + 1 - a2.y);
	} else {
		dpbc=fabs(a2.y + 1 - a1.y);
	}
	if(dpbc < dy) {
		dy=dpbc;
	}

	if(a1.z <= a2.z) {
		dpbc=fabs(a1.z + 1 - a2.z);
	} else {
		dpbc=fabs(a2.z + 1 - a1.z);
	}
	if(dpbc < dz) {
		dz=dpbc;
	}

	//Convert distance from u.c to amstrongs
	dx*=a;
	dy*=b;
	dz*=c;

	//Compute distance
	d=sqrt(dx*dx + dy*dy + dz*dz);

	//Return whether or not connected
	if(d < 2.56) {
		return 1;
	} else {
		return 0;
	}

}

int** get_amatrix_Tatoms (atom *list, int *Tids) {

	int i, j, k, ctl, index;
	int **M_all, **M_T;

	//Find T ids
	index=0;
	for(i=0; i<Natoms; i++) {
		if(strcmp(list[i].id,"Si") ==0 ) {
			Tids[index]=i;
			index++;
		}
	}


	//Save memory for adjacency matries
	M_all=(int**)calloc(Natoms, sizeof(int*));
	for(i=0; i<Natoms; i++) {
		M_all[i]=(int*)calloc(Natoms, sizeof(int));
	}

	M_T=(int**)calloc(Tatoms, sizeof(int*));
	for(i=0; i<Tatoms; i++) {
		M_T[i]=(int*)calloc(Tatoms, sizeof(int));
	}

	//Compute adjacency matrix all atoms
	for(i=0; i<Natoms-1; i++) {
		for(j=i+1; j<Natoms; j++) {
			M_all[i][j]=connectivity(list[i], list[j]);
			M_all[j][i]=M_all[i][j];
		}
	}

	//Compute adjacency matrix T atoms
	for(i=0; i<Tatoms-1; i++) {
		for(j=i+1; j<Tatoms; j++) {
			ctl=1;
			for(k=0; k<Natoms && ctl==1; k++) {
				if(M_all[k][Tids[i]]==1 && M_all[k][Tids[j]]==1) {
					M_T[i][j]=1;
					M_T[j][i]=1;
					ctl=0;
				}
			}
		}
	}

	//Free memory
	for(i=0; i<Natoms; i++){
		free(M_all[i]);
	}
	free(M_all);

	return M_T;
}



vector<vector<int> > get_neighbours(int **M_T) {

	//For each atom find its neighbours
	vector<vector<int> > neigbrs(Tatoms, vector<int>(0));
	for(int i=0; i<Tatoms-1; i++) {
		for(int j=i+1; j<Tatoms; j++) {
			if(M_T[i][j]==1) {
				neigbrs[i].push_back(j);
				neigbrs[j].push_back(i);
			}
		}
	}

	return neigbrs;
}


int find_neigb_random(int n, int **M_T) {

	int nei;

	vector<int> possible_neis(0);

	//Find all neighbours
	for(int i=0; i<Tatoms; i++) {
		if(M_T[n][i]==1) {
			possible_neis.push_back(i);
		}
	}

	//Check that there are available neighbours
	if(possible_neis.size() == 0) {
		cerr << "No more neighbours available" << endl;
		cerr << "Structure could not be generated" << endl;
		exit(1);
	}

	//Select randomly one neighbour
	int index=rand()%possible_neis.size();
	nei=possible_neis[index];

	return nei;
}


void delete_neigb(int n, int **M_T) {

	for(int i=0; i<Tatoms; i++) {
		M_T[i][n]=0;
		M_T[n][i]=0;
	}

	return;
}


vector<int> generate_chains (int **M_T, vector<vector<int> > neigbrs, vector<int> chains) {


    int index, next, last;
	vector<int> Als(0), Tpos(Tatoms);
	vector<vector<int> > neigbrs_available=neigbrs;

	//Initialize vector containing all possible atoms to be substituted
	for(int i=0; i<Tatoms; i++) {
		Tpos[i]=i;
	}

	for(unsigned int i=0; i<chains.size(); i++) {

		//Check that there are more root atoms
		if(Tpos.size() == 0) {
			cerr << "No more root atoms available" << endl;
			cerr << "Structure could not be generated" << endl;
			exit(1);
		}

		//Find next root
		index=rand()%Tpos.size();
		Als.push_back(Tpos[index]);

		//Remove from possible roots
		int save=Tpos[index];
		Tpos.erase(Tpos.begin()+index);


		//Remove neighbours from possible roots
		for(unsigned int j=0; j<neigbrs_available[save].size(); j++) {
			int ctl=1;
			for(int k=0; k<Tpos.size() && ctl==1; k++) {
				if(neigbrs_available[save][j] == Tpos[k]) {
					Tpos.erase(Tpos.begin()+k);
					ctl=0;
				}
			}
		}


		//Complete the chain
		last=save;
		for(int j=1; j<chains[i]; j++) {

			//Find a random neighbour of the root
			index=find_neigb_random(last, M_T);
			//Remove from possible root
			int ctl=1;
			for(int k=0; k<Tpos.size() && ctl==1; k++) {
				if(index == Tpos[k]) {
					Tpos.erase(Tpos.begin()+k);
					ctl=0;
				}
			}

			//Remove neighbours from possible roots
			for(unsigned int r=0; r<neigbrs_available[index].size(); r++) {
				ctl=1;
				for(int k=0; k<Tpos.size() && ctl==1; k++) {
					if(neigbrs_available[index][r] == Tpos[k]) {
						Tpos.erase(Tpos.begin()+k);
						ctl=0;
					}
				}
			}

			//Remove all the other neighbours from neigbrs_available
			//Remove links from neighbours
			for(unsigned int r=0; r<neigbrs_available[last].size(); r++) {
				if(neigbrs_available[last][r] != index) {
					//delete rows and columns in adjacency matrix
					delete_neigb(neigbrs_available[last][r], M_T);
				}
			}

			//Remove links from 0
			M_T[last][index]=0;
			M_T[index][last]=0;

			//Add new Al
			Als.push_back(index);
			last=index;
			
			//If it is the end of the chain
			if(j==chains[i]-1) {
				//Remove neighbours of the index from neigbrs_available
				for(unsigned int r=0; r<neigbrs_available[index].size(); r++) {
					//delete rows and columns in adjacency matrix
					delete_neigb(neigbrs_available[index][r], M_T);
				}
			}

		}
	}

	return Als;
}




//Clusters
vector<int> clusters_substitutions (int **M_T, int Nsubst) {

	int index, ini_index, ini, next, ctl;
	vector<int> Als(0);
	vector<vector<int> > neigbrs;

	//Compute neighbours of all atoms
	neigbrs=get_neighbours(M_T);

	//Find First substitution
	index=rand()%Tatoms;
	//Save Al
	Als.push_back(index);

	//Do the rest of substitutions
	int currentAl=0, currentNei=0, countNsubst=1;
	while(countNsubst < Nsubst) {
		//Check if current Al atom has neighbours available
		if(in(neigbrs[Als[currentAl]][currentNei], Als) == false) {
			//Add Al atom
			Als.push_back(neigbrs[Als[currentAl]][currentNei]);
			countNsubst++;
			//Update neighbour
			currentNei++;
			//Check if there are more neighbours
			if(neigbrs[Als[currentAl]].size() <= currentNei) {
				currentAl++;
				currentNei=0;
			}
		} else {
			//Update neighbour
			currentNei++;
			//Check if there are more neighbours
			if(neigbrs[Als[currentAl]].size() <= currentNei) {
				currentAl++;
				currentNei=0;
			}			
		}
	}


	return Als;

}

bool in(int index, vector<int> Als) {

	bool res;

	for(unsigned int i=0; i<Als.size(); i++) {
		if(Als[i] == index) {
			return true;
		}
	}

	return false;
}


//MERW
vector<int> merw_substitutions (int **M_T, int Nsubst) {

	int index, ini_index, ini, next, ctl;
	double **S;
	vector<int> Als(0), Tpos(Tatoms);
	vector<vector<int> > neigbrs;
	int permute=0;

	if(Nsubst > (Tatoms/2)) {
		Nsubst = Tatoms - Nsubst; 
		permute=1;
	}

	//Save memory for the probability matrix (S will change every time we add an Al)
	S=(double**)calloc(Tatoms, sizeof(double*));
	for(int i=0; i<Tatoms; i++) {
		S[i]=(double*)calloc(Tatoms, sizeof(double));
	}

	//Initialize vector containing all possible atoms to be substituted
	for(int i=0; i<Tatoms; i++) {
		Tpos[i]=i;
	}

	//Find First substitution
	index=rand()%Tpos.size();
	//Save Al
	Als.push_back(Tpos[index]);
	//Remove atom bonds from adjacency matrix
	delete_neigb(Tpos[index], M_T);
	//Compute new neighbours of all atoms
	neigbrs=get_neighbours(M_T);
	//Remove atom from possible future substitutions
	Tpos.erase(Tpos.begin()+index);	

	//Do the rest of substitutions
	for(int i=1; i<Nsubst; i++) {
		//Check that there are more available atoms
		if(Tpos.size() == 0) {
			cerr << "No more atoms available" << endl;
			cerr << "Structure could not be generated" << endl;
			exit(1);
		}

		//Compute transition probability matrix
		compute_S(M_T, S);

		//Select randomly the initial position of the MERW (not issolated)
		int count=0;
		do {
			ini_index=rand()%Tpos.size();  //index of the initial point in Tpos
			ini=Tpos[ini_index];		   //id of the initial position
			count++;
		} while((neigbrs[ini].size() == 0) && (count < 20));

		if(count==20) {
			//Save next Al
			Als.push_back(Tpos[ini_index]);
			//Change adjacency matrix
			delete_neigb(Tpos[ini_index], M_T);
			//Compute new neighbours
			neigbrs=get_neighbours(M_T);
			//Remove atom from possible substitution
			Tpos.erase(Tpos.begin()+ini_index);	
		} else {

			//Perform MERW: equilibration and finding next substitution (next)
			next=merw(S, neigbrs, ini);

			//Find the next atom in the Tpos vector
			int l;   //index of next in Tpos
			ctl=0;
			for(unsigned int k=0; k<Tpos.size() && ctl==0; k++) {
				if(Tpos[k] == next) {
					l=k;
					ctl=1;
				}
			}

			//Error check
			if(ctl==0) {
				cerr << "The atom next was not found in the remaining atoms:" << endl;
				cerr << "That means our RW is walking through Al atoms, and that shouldn't happen." << endl;
				exit(1);			
			}

			//Save next Al
			Als.push_back(Tpos[l]);
			//Change adjacency matrix
			delete_neigb(Tpos[l], M_T);
			//Compute new neighbours
			neigbrs=get_neighbours(M_T);
			//Remove atom from possible substitution
			Tpos.erase(Tpos.begin()+l);	
		}

	}

	//Free memory
	for(int i=0; i<Tatoms; i++) {
		free(S[i]);
	}
	free(S);


	//Check if we need to permute Si for Al
	if(permute==1) {

		vector<int> Als2(0);
		for(int i=0; i<Tatoms; i++) {
			int ex=0;
			for(int j=0; j<Nsubst && ex==0; j++) {
				if(Als[j] == i) {
					ex=1;
				}
			}
			if(ex==0) {
				Als2.push_back(i);
			}
		}

		return Als2;
	} else {

		return Als;
	}

}

void compute_S(int **M_T, double **S) {

	int max_ind;
	double max_eval;
	Eigen::EigenSolver<Eigen::MatrixXf> eigensolver;
	Eigen::MatrixXf A = Eigen::MatrixXf(Tatoms,Tatoms);

	//Convert adjacency matrix
	for(int i=0; i<Tatoms; i++) {
		for(int j=0; j<Tatoms; j++) {
			double p=((double) rand() / (RAND_MAX));
			A(i,j)=M_T[i][j] + p*0.01;
		}
	}
	

	//Diagonalize A
	eigensolver.compute(A, /* computeEigenvectors = */ true);
	Eigen::VectorXf eigen_values = eigensolver.eigenvalues().real();
	Eigen::MatrixXf eigen_vectors = eigensolver.eigenvectors().real();

	//Get maximum eigenvalue
	max_eval=eigen_values(0);
	max_ind=0;
	for(int i=1; i<Tatoms; i++) {
		if(eigen_values(i) > max_eval) {
			max_eval=eigen_values(i);
			max_ind=i;
		}
	}

	//Compute probability matrix
	for(int i=0; i<Tatoms; i++) {
		for(int j=0; j<Tatoms; j++) {
			S[i][j]= (1.0*M_T[i][j]/(max_eval)) * (eigen_vectors(j,max_ind)/eigen_vectors(i,max_ind));
		}
	}

	return;
}


int merw(double **S, vector<vector<int> > neigbrs, int ini) {

	int next=-1, i, j; 
	bool done=false;
	unsigned int Noptions;
	vector<int> visits(Tatoms, 0);

	//Equilibration
	i=ini;	//i is the current position of the walker
	for(int k=0; k<Neqsteps; k++) {
		//find next position in the RW
		j=merw_step(S, neigbrs, i);
		//Update move
		i=j;
	}

	//Find new substitution
	while(done == false) {
		//find next position in the RW
		j=merw_step(S, neigbrs, i);
		//Update move
		i=j;
		//Update visit
		visits[i]++;
		//If maximum number of visits reached, return position
		if(visits[i] >= Nvisits) {
			done=true;
			next=i;
		}

	}

	//Error check
	if(next == -1) {
		cerr << "Next Al site has not been found" << endl;
		exit(1);
	}

	return next;
}

int merw_step(double **S, vector<vector<int> > neigbrs, int i) {

	int next, ctl=0;
	double acc=0;
	unsigned int Noptions;

	//Find prob to go to each neighbour
	Noptions=neigbrs[i].size();
	vector<double> probs(Noptions);

	for(int j=0; j<Noptions; j++) {
		probs[j]=S[i][neigbrs[i][j]];
	}

	//Generate uniform number p from (0,1)
	double p=((double) rand() / (RAND_MAX));

	//Choose randomly an element of neigbrs[i] with probability probs. Call it j
	for(int j=0; j<Noptions-1 && ctl==0; j++) {
		acc+=probs[j];
		if(p <= acc) {
			next=neigbrs[i][j];
			ctl=1;
		}
	}
	if (ctl==0) {
		next=neigbrs[i][Noptions-1];
	}

	return next;	
}



//Random
vector<int> generate_random (int Nsubst) {

	vector<int> Tpos(Tatoms), Als(Nsubst);
	int index;

	//Initialize vector containing all possible atoms to be substituted
	for(int i=0; i<Tatoms; i++) {
		Tpos[i]=i;
	}

	for(int k=0; k<Nsubst; k++) {
		//Find next Al
		index=rand()%Tpos.size();
		Als[k]=Tpos[index];

		//Remove from possible roots
		Tpos.erase(Tpos.begin()+index);
	}


	return Als;
}



void print_structure (atom *list, vector<int> Als, int struc, string name_zeo, string name_alg, string out_name) {
	
	int alcount, ctl;
	string fname, s;
	ofstream fout;

	fname = out_name + "/" + name_zeo + "_" + name_alg + "_" + to_string(struc) + ".cif";

    fout.open(fname.c_str());    
    if(fout.fail()) {
        cerr << "unable to open file " << fname.c_str() << " for reading" << endl;
        exit( 1 );
    }

	fout << setprecision(3) << fixed;

	fout << "data_" << name_zeo << endl << endl;

	fout << "_audit_creation_method RASPA-1.0" << endl;
	fout << "_audit_creation_date 2022-4-23" << endl;
	fout << "_audit_author_name 'Pablo_Romero'"<< endl << endl;

	fout << "_cell_length_a    " << a << endl;
	fout << "_cell_length_b    " << b << endl;
	fout << "_cell_length_c    " << c << endl;
	fout << "_cell_angle_alpha " << alpha << endl;
	fout << "_cell_angle_beta  " << beta << endl;
	fout << "_cell_angle_gamma " << gama << endl;
	fout << "_cell_volume      " << a*b*c << endl << endl;

	fout << "_symmetry_cell_setting          " << setting << endl;
	fout << "_symmetry_space_group_name_Hall 'P 1'" << endl;
	fout << "_symmetry_space_group_name_H-M  'P 1'" << endl;
	fout << "_symmetry_Int_Tables_number     1" << endl;

	fout << "_symmetry_equiv_pos_as_xyz 'x,y,z'" << endl << endl;

	fout << "loop_" << endl;
	fout << "_atom_site_label" << endl;
	fout << "_atom_site_type_symbol" << endl;
	fout << "_atom_site_fract_x" << endl;
	fout << "_atom_site_fract_y" << endl;
	fout << "_atom_site_fract_z" << endl;
	fout << "_atom_site_charge" << endl;

	alcount=1;

	int Tat=-1;
	for(int i=0; i<Natoms; i++) {
		if(strcmp(list[i].id,"Si") ==0 ) {
			Tat++;
			ctl=1;
			for(unsigned int j=0; j<Als.size() && ctl==1; j++) {
				if(Tat == Als[j]) {
					//Put Al
					fout << "Al" << alcount << "        Al     " << list[i].x << setw(10) << list[i].y << setw(10) << list[i].z << setw(10) << list[i].q << endl;  
					ctl=0;
					alcount++;
				}
			}
			//Put Si
			if(ctl == 1) {
				fout << list[i].at << setw(10) << list[i].id << setw(10) << list[i].x << setw(10) << list[i].y << setw(10) << list[i].z << setw(10) << list[i].q << endl;
			}
		} else {
			fout << list[i].at << setw(10) << list[i].id << setw(10) << list[i].x << setw(10) << list[i].y << setw(10) << list[i].z << setw(10) << list[i].q << endl;
		}
	}

	return;
}


