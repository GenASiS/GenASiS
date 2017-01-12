#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 19
//#define DEBUG 1

int nrho;
int ntemp;
int nye;

double *alltables;
double *logrho; 
double *logtemp;
double *yes;
double energy_shift;
double dtemp, dtempi;
double drho, drhoi;
double dye, dyei;

// min and max values

double eos_rhomax, eos_rhomin;
double eos_tempmin, eos_tempmax;
double eos_yemin, eos_yemax;

// table key
// 0 logpress 
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma

// some vectors for selecting variables for more
// efficient interpolation
int ivs_short[8];

// frontend function declarations
void nuc_eos_C_short(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs, double* xent,
		     double *xcs2, double* xdedt, double* xdpderho,
		     double *xdpdrhoe, double* xmunu, int keytemp,
		     int *keyerr,double rfeps); 

// core function declarations

void nuc_eos_C_ReadTable(char* nuceos_table_name);

void nuc_eos_C_linterp_many(double x, double y, double z,
			    double* f, double* ft, 
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt);

void nuc_eos_C_linterp_some(double x, double y, double z,
			    double* f, double* ft, 
			    int* ivs,
			    int nx, int ny, int nz, int nvars,
			    double* xt,double*yt, double* zt);

void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double* f, double* ft, 
				int nx, int ny, int nz, 
				double* xt, double*yt, double* zt,
				double* linterp_for_temp);

void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
			double leps, double prec, double *lt,
			int *keyerr);

void nuc_eos_C_testing();

