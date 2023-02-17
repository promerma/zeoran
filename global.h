//Unit cell parameters
int Natoms, Tatoms;
double a, b, c;

//Default merw parameters
int Neqsteps=100;	//Number of steps to find equilibrium position
int Nvisits=20;		//Number of visits needed to select a new Al

//atom structure
typedef struct{
  char id[5], at[5];
  double x, y, z, q;
} atom;