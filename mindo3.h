
#include <vector>
#include <windows.h>
#include "newmatap.h"          // newmat headers including advanced functions

using namespace std;

// basis function
/*struct bfunc {
  int index, 
      type;				// s,px,py,pz
  double u,ip;				// s and p AO 1e one-center integrals
  double exponent[6],			// basis
         coefficient[6];
  struct atom_mindo *atom;		// pointer to atom
};*/

struct bfnc {
  int    natom,				// atom number
  	 type;				// s,px,py,pz
  double u,ip;				// s and p AO 1e one-center integrals
  double expn[6],			// basis
         coef[6];
};

// atom structure
struct atom_mindo { 
  int atno,Z; 				// atomic number and charge
  double xyz[3];			// cartesian 
  double rho, Eref, Hf,
         gss, gsp, gpp, gppp,
         hsp, hppp;
  int nbf;				// number of basis functions
//  bfunc basis[4];			// basis functions, maximum 4
};

typedef vector <struct atom_mindo> atomvector_mindo;

typedef unsigned int uint;

class mindo3 {

  public:

    int    MaxSCF, SCFit, MO_count;
    double TolSCF, Ediff;
    BOOL   D_initialized;
    BOOL   F0_initialized;
    
    mindo3 (atomvector_mindo *atoms);
    double Energy ();
    void   forces ( Matrix *forces );
    void   num_forces_right   ( Matrix *forces, double dx );
    void   num_forces_left    ( Matrix *forces, double dx );
    void   num_forces_central ( Matrix *forces, double dx );
    void   num_hessian_right  ( SymmetricMatrix *hess, double dx );
    void   num_hessian_central( SymmetricMatrix *hess, double dx );
    void   print_orbs ();
    void   print_grads( Matrix G);
    void   print_hess ( SymmetricMatrix H);

    double gto(bfnc bf, double xyz[3]);
    double molorb(int numorb, double xyz[3]);
    double MOEnergy(int imo);

    double power(double x, double y);

  private:

    atomvector_mindo *atoms;

    DiagonalMatrix Orbe; 
    SquareMatrix Orbs, D;
//    SymmetricMatrix F0,F1,F2,F;
    int	   nclosed, nbf, nat;

    void   init ();
    void   init_parameters ();

    double overlap(bfnc bfi, bfnc bfj); 
    double doverlap(bfnc bfi, bfnc bfj, int dir);

    double dist2(double A[3],double B[3]);
    double gamma(struct atom_mindo ati, struct atom_mindo atj);
    double scale(int atnoi,int atnoj,double R);

    int    numel (int charge);
    double enuke ();
    int    numbf ();
    double refeng();

    double g ( bfnc bfi, bfnc bfj );
    double h ( bfnc bfi, bfnc bfj );

    void   guess_D ();
    void   calc_F0 ();
    void   calc_F1 ();
    void   calc_F2 ();
    void   mkdens  ();

    double SCF     ();

};

// orbitals names
const char or_nam[4][3] = {"S","PX","PY","PZ"};

// MINDO/3 Parameters: Thru Ar in eV

//s and p atomic orbital one-electron one-center integrals
const double Uss[19] = { 0.0, -12.505, 0.0,
        0.0, 0.0, -33.61, -51.79, -66.06,
        -91.73, -129.86, 0.0,
        0.0, 0.0, 0.0, -39.82, -56.23, -73.39, -98.99, 0.0}; 
const double Upp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, -25.11, -39.18, -56.40, -78.80, -105.93, 0.0,
        0.0, 0.0, 0.0, -29.15, -42.31, -57.25, -76.43, 0.0};

// s-s atomic orbital one center two electron repulsion integral
const double gss[19] = { 0.0, 12.848, 0.0,
        0.0, 0.0, 10.59, 12.23, 13.59, 15.42, 16.92, 0.0,
        0.0, 0.0, 0.0, 9.82, 11.56, 12.88, 15.03, 0.0};

// s-p atomic orbital one center two electron repulsion integral
const double gsp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 9.56, 11.47, 12.66, 14.48, 17.25, 0.0,
        0.0, 0.0, 0.0, 8.36, 10.08, 11.26, 13.16, 0.0};

// p-p atomic orbital one center two electron repulsion integral
const double gpp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 8.86, 11.08, 12.98, 14.52, 16.71, 0.0,
        0.0, 0.0, 0.0, 7.31, 8.64, 9.90, 11.30, 0.0};

// p-p' atomic orbital one center two electron repulsion integral
const double gppp[19] = { 0.0, 0.0, 0.0,
         0.0, 0.0, 7.86, 9.84, 11.59, 12.98, 14.91, 0.0,
         0.0, 0.0, 0.0, 6.54, 7.68, 8.83, 9.97, 0.0};

// s-p atomic orbital one-center two-electron exchange integral
const double hsp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, 1.81, 2.43, 3.14, 3.94, 4.83, 0.0,
        0.0, 0.0, 0.0, 1.32, 1.92, 2.26, 2.42, 0.0};

const double hppp[19] = { 0.0, 0.0, 0.0,
         0.0, 0.0, 0.50, 0.62, 0.70, 0.77, 0.90, 0.0,
         0.0, 0.0, 0.0, 0.38, 0.48, 0.54, 0.67, 0.0};

// averaged repulsion integral for use in gamma
const double f03[19] = { 0.0, 12.848, 10.0, 				
        10.0, 0.0, 8.958, 10.833, 12.377, 13.985, 16.250,
        10.000, 10.000, 0.000, 0.000,7.57 ,  9.00 ,10.20 , 11.73};

// s atomic orbital ionization potential for two center resonance integral term
const double IPs[19] = { 0.0, -13.605, 0.0,
        0.0, 0.0, -15.160, -21.340, -27.510, -35.300, -43.700, -17.820,
        0.0, 0.0, 0.0, 0.0, -21.100, -23.840, -25.260, 0.0};

// p atomic orbital ionization potential for two center resonance integral term
const double IPp[19] = { 0.0, 0.0, 0.0,
        0.0, 0.0, -8.520, -11.540, -14.340, -17.910, -20.890, -8.510,
        0.0, 0.0, 0.0, 0.0, -10.290, -12.410, -15.090, 0.0};

// s-type Slater atomic orbital exponent
const double zetas[19] = { 0.0, 1.30, 0.0,
          0.0, 0.0, 1.211156, 1.739391, 2.704546, 3.640575, 3.111270, 0.0,
          0.0, 0.0, 0.0, 1.629173, 1.926108, 1.719480, 3.430887, 0.0};

// p-type Slater atomic orbital exponent
const double zetap[19] = { 0.0, 0.0, 0.0,
          0.0, 0.0, 0.972826, 1.709645, 1.870839, 2.168448, 1.419860, 0.0,
          0.0, 0.0, 0.0, 1.381721, 1.590665, 1.403205, 1.627017, 0.0};

// Atomic heat of formations: Mopac got from CRC
const double Hfat[19] = { 0.0, 52.102, 0.0,
         0.0, 0.0, 135.7, 170.89, 113.0, 59.559, 18.86, 0.0,
         0.0, 0.0, 0.0, 106.0, 79.8, 65.65, 28.95, 0.0};

// Default isolated atomic energy values from Mopac:EISOL3
const double Eat[19] = {0.0, -12.505, 0.0,
       0.0 ,0.0,-61.70,-119.47,-187.51,-307.07,-475.00,0.0,
       0.0,0.0,0.0,-90.98,-150.81,-229.15,-345.93,0.0};

const double nbfat[19] = { 0.0, 1, 0.0,
          0.0, 0.0, 4, 4, 4, 4, 4, 0.0,
          0.0, 0.0, 0.0, 4, 4, 4, 4, 0.0};

const double CoreQ[19] = { 0.0, 1, 0.0,
          0.0, 0.0, 3, 4, 5, 6, 7, 0.0,
          0.0, 0.0, 0.0, 4, 5, 6, 7, 0.0};

const double atmas[19] = { 0.0,
  1.0, 0,0,0,0,0,0, 16.0, 0,0,0,0,0,0,0,0,0,0 };


/*
# Gaussian functions for fitting to Slaters. These functions are
# STO-6G fits to slater exponents with exponents of 1. To fit
# to exponents of \zeta, you need only multiply each
# exponent by \zeta^2
# The rest of these functions can be obtained from Stewart,
#  JCP 52, 431 (1970)
*/
//http://www.chem.arizona.edu/~lichtend/c518.dir/s2001.dir/basissets.dir/basisb.html

const int stong = 6;				// STO-6G

#define None      {0,		0,		0,		0,		0,		0}
#define gexps_1s  {2.310303149e01   ,4.235915534e00   ,1.185056519e00  ,4.070988982e-01,1.580884151e-01,6.510953954e-02}
#define gcoefs_1s {9.163596280e-03  ,4.936149294e-02  ,1.685383049e-01 ,3.705627997e-01,4.164915298e-01,1.303340841e-01}
#define gexps_2s  {2.768496241e01   ,5.077140627e00   ,1.426786050e00  ,2.040335729e-01,9.260298399e-02,4.416183978e-02}
#define gcoefs_2s {-4.151277819e-03 ,-2.067024148e-02 ,-5.150303337e-02,3.346271174e-01,5.621061301e-01,1.712994697e-01}
#define gexps_2p  {5.868285913e00   ,1.530329631e00   ,5.475665231e-01 ,2.288932733e-01,1.046655969e-01,4.948220127e-02}
#define gcoefs_2p {7.924233646e-03  ,5.144104825e-02  ,1.898400060e-01 ,4.049863191e-01,4.012362861e-01,1.051855189e-01}
#define gexps_3s  {3.273031938e00   ,9.200611311e-01  ,3.593349765e-01 ,8.636686991e-02,4.797373812e-02,2.724741144e-02}
#define gcoefs_3s {-6.775596947e-03 ,-5.639325779e-02 ,-1.587856086e-01,5.534527651e-01,5.015351020e-01,7.223633674e-02}
#define gexps_3p  {5.077973607e00   ,1.340786940e00   ,2.248434849e-01 ,1.131741848e-01,6.076408893e-02,3.315424265e-02}
#define gcoefs_3p {-3.329929840e-03 ,-1.419488340e-02 ,1.639395770e-01 ,4.485358256e-01,3.908813050e-01,7.411456232e-02}

// principle quantum number N
const int NQN[19] = { 0, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3 };

const double gexps[3][2][6] = {			// indexed by N,s_or_p:
  {gexps_1s,None    },  // N=1
  {gexps_2s,gexps_2p},  // N=2
  {gexps_3s,gexps_3p},  // N=3
}; //  s       p
    
    
const double gcoefs[3][2][6] = {		// indexed by N,s_or_p:
  {gcoefs_1s,None     },  // N=1
  {gcoefs_2s,gcoefs_2p},  // N=2
  {gcoefs_3s,gcoefs_3p},  // N=3
}; //  s       p

const int s_or_p[4] = {0,1,1,1};		// whether the func is s or p type, based on the L QN
