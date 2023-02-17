//Read data
atom* read_atom_sites (string fname);
void read_unit_cell (string fname);
//Generate adjacency matrix and neighbours
int connectivity(atom a1, atom a2);
int** get_amatrix_Tatoms(atom *list, int *Tids);
vector<vector<int> > get_neighbours(int **M_T);
//Generate Al substitutions: chains
int find_neigb_random(int n, int **M_T);
void delete_neigb(int n, int **M_T);
vector<int> generate_chains (int **M_T, vector<vector<int> > neigbrs, vector<int> chains);
//Generate Al substitutions: clusters
vector<int> clusters_substitutions (int **M_T, int Nsubst);
//Generate Al substitutions: merw
vector<int> merw_substitutions (int **M_T, int Nsubst);
void compute_S(int **M_T, double **S);
int merw(double **S, vector<vector<int> > neigbrs, int ini);
int merw_step(double **S, vector<vector<int> > neigbrs, int i);
bool in(int index, vector<int> Als);
//Generate Al substitutions: random
vector<int> generate_random (int Nsubst);
//Print generated structure
void print_structure (atom *list, vector<int> Als, int struc, string name_zeo, string name_alg, string out_name);