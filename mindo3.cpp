
// http://www.cachesoftware.com/mopac/Mopac2002manual/node439.html

#define WANT_STREAM            // include iostream and iomanipulators
#include "newmatap.h"          // newmat headers including advanced functions
#include "newmatio.h"          // newmat headers including output functions

#include <cmath>
#include "Constants.h"
#include "mindo3.h"

#include "ViewMol3D.h"

double mindo3::power(double x, double y) {
  if ((x==0) && (y==0)) return 1.0;
  if  (x==0)            return 0.0;
  if  (x==1)            return 1.0;
  if  (y==2)		return x*x;
  return pow(fabs(x),y);
}

double axy[19][19], Bxy[19][19];

mindo3::mindo3(atomvector_mindo *_atoms) {
  atoms = _atoms;
  init();
}

void mindo3::init_parameters() {
  int i, j;

// axy Core repulsion function terms
  axy[1][1 ] = 1.489450; axy[1][5 ] = 2.090352; axy[1][6 ] = 1.475836; axy[1][7 ] = 0.589380; axy[1][8 ] = 0.478901; 
  axy[1][9 ] = 3.771362; axy[1][14] = 0.940789; axy[1][15] = 0.923170; axy[1][16] = 1.700689; axy[1][17] = 2.089404; 
  axy[5][5 ] = 2.280544; axy[5][6 ] = 2.138291; axy[5][7 ] = 1.909763; axy[5][8 ] = 2.484827; axy[5][9 ] = 2.862183; 
  axy[6][6 ] = 1.371208; axy[6][7 ] = 1.635259; axy[6][8 ] = 1.820975; axy[6][9 ] = 2.725913; axy[6][14] = 1.101382;  
  axy[6][15] = 1.029693; axy[6][16] = 1.761370; axy[6][17] = 1.676222; axy[7][7 ] = 2.209618; axy[7][8 ] = 1.873859; 
  axy[7][9 ] = 2.861667; axy[8][8 ] = 1.537190; axy[8][9 ] = 2.266949; axy[9][9 ] = 3.864997; axy[14][14]= 0.918432;       
  axy[15][15]= 1.186652; axy[16][16]= 1.751617; axy[17][17]= 1.792125;      

  for(i=0;i<19;i++) for(j=i+1;j<19;j++) axy[j][i] = axy[i][j];

// Diatomic two center one-electron resonance integral multiplier
  Bxy[1 ][1 ] = 0.244770;  Bxy[1 ][5 ] = 0.185347;  Bxy[1 ][6 ] = 0.315011;  Bxy[1 ][7 ] = 0.360776;  Bxy[1 ][8 ] = 0.417759; 
  Bxy[1 ][9 ] = 0.195242;  Bxy[1 ][14] = 0.289647;  Bxy[1 ][15] = 0.320118;  Bxy[1 ][16] = 0.220654;  Bxy[1 ][17] = 0.231653;
  Bxy[5 ][5 ] = 0.151324;  Bxy[5 ][6 ] = 0.250031;  Bxy[5 ][7 ] = 0.310959;  Bxy[5 ][8 ] = 0.349745;  Bxy[5 ][9 ] = 0.219591;
  Bxy[6 ][6 ] = 0.419907;  Bxy[6 ][7 ] = 0.410886;  Bxy[6 ][8 ] = 0.464514;  Bxy[6 ][9 ] = 0.247494;  Bxy[6 ][14] = 0.411377; 
  Bxy[6 ][15] = 0.457816;  Bxy[6 ][16] = 0.284620;  Bxy[6 ][17] = 0.315480;  Bxy[7 ][7 ] = 0.377342;  Bxy[7 ][8 ] = 0.458110; 
  Bxy[7 ][9 ] = 0.205347;  Bxy[8 ][8 ] = 0.659407;  Bxy[8 ][9 ] = 0.334044;  Bxy[9 ][9 ] = 0.197464;  Bxy[14][14] = 0.291703; 
  Bxy[15][15] = 0.311790;  Bxy[16][16] = 0.202489;  Bxy[17][17] = 0.258969;

  for(i=0;i<19;i++) for(j=i+1;j<19;j++) Bxy[j][i] = Bxy[i][j];

}      

vector <bfnc> bfncs;

void mindo3::init() {

  init_parameters();

  int ibf = 0, i, j, na=0;

  bfncs.erase( bfncs.begin(), bfncs.end() );

  for ( vector <struct atom_mindo>::iterator at = atoms->begin(); at!=atoms->end(); at++, na++ ) {
     (*at).Z	= CoreQ [(*at).atno];
     (*at).rho  = e2/f03[(*at).atno];
     (*at).nbf  = nbfat [(*at).atno];
     (*at).Eref = Eat   [(*at).atno];
     (*at).Hf	= Hfat  [(*at).atno];
     (*at).gss  = gss   [(*at).atno];
     (*at).gsp  = gsp   [(*at).atno];
     (*at).gpp  = gpp   [(*at).atno];
     (*at).gppp = gppp  [(*at).atno];
     (*at).hsp  = hsp   [(*at).atno];
     (*at).hppp = hppp  [(*at).atno];
     for (i=0; i<(*at).nbf; i++) {
        struct bfnc bfnc;
        double zeta;
        bfnc.natom = na;
        bfnc.type  = i;
        if ( i==0 ) {
            zeta	= zetas[(*at).atno];
            bfnc.u	= Uss  [(*at).atno];
            bfnc.ip	= IPs  [(*at).atno];
        } else {
            zeta	= zetap[(*at).atno];
            bfnc.u	= Upp  [(*at).atno];
            bfnc.ip	= IPp  [(*at).atno];
        }
        for (j=0;j<stong;j++) {
            bfnc.expn[j] = gexps  [ NQN[(*at).atno]-1 ] [ s_or_p[i] ] [ j ] * zeta * zeta;
            bfnc.coef[j] = gcoefs [ NQN[(*at).atno]-1 ] [ s_or_p[i] ] [ j ];
        }
        bfncs.push_back(bfnc);
     }
  }

  MaxSCF = 128;
  TolSCF = 1e-9;//1.0e-7;
  D_initialized  = FALSE;
  F0_initialized = FALSE;
  MO_count = ibf-1;

  nclosed = numel(0)/2;
  nbf = numbf();
  nat = atoms->size();

}

double mindo3::dist2(double A[3],double B[3]) {
    return power(A[0]-B[0],2)+power(A[1]-B[1],2)+power(A[2]-B[2],2);
}

// Number of electrons in the system
int mindo3::numel(int charge) {
    int nel = 0;
    for ( vector <struct atom_mindo>::iterator at = atoms->begin(); at!=atoms->end(); at++ ) {
      nel += (*at).Z;
    }
    return nel-charge;
}

// Coulomb repulsion that goes to the proper limit at R=0"
double mindo3::gamma(struct atom_mindo ati, struct atom_mindo atj) {
    return e2/sqrt( dist2( ati.xyz, atj.xyz ) + 0.25*power(ati.rho + atj.rho,2) );
}

// Prefactor from the nuclear repulsion term
double mindo3::scale(int atnoi,int atnoj,double R) {
    double alpha = axy[atnoi][atnoj];					// Part of the scale factor for the nuclear repulsion
    if (atnoi == 1) {
        if ((atnoj == 7) || (atnoj == 8))
            return alpha*exp(-R);
    } else 
      if (atnoj == 1)
        if ((atnoi == 7) || (atnoi == 8))
            return alpha*exp(-R);
    return exp(-alpha*R);
}

// Compute the nuclear repulsion energy
double mindo3::enuke() {
    double enuke = 0.0;
    for ( int i=0; i<nat; i++ ) {
      struct atom_mindo ati = (*atoms)[i];
      for ( int j=i+1; j<nat; j++ ) {
            struct atom_mindo atj = (*atoms)[j];
            double R = sqrt( dist2( ati.xyz, atj.xyz ) );
            double sc = scale( ati.atno, atj.atno, R );
            double gammaij = gamma( ati, atj );
            enuke += ati.Z*atj.Z*(gammaij + (e2/R-gammaij)*sc);
       }
    }
    return enuke;
}

// Number of basis functions in an atom_mindo list
int mindo3::numbf() {
    int nbf = 0;
    for ( vector <struct atom_mindo>::iterator at = atoms->begin(); at!=atoms->end(); at++ )
    	nbf += (*at).nbf;
    return nbf;
}

// Ref = heat of formation - energy of atomization
double mindo3::refeng() {
    double eat = 0.0;
    double hfat = 0.0;
//    int nat = atoms->size();
    for ( int i=0; i<nat; i++ ) {
        eat  += (*atoms)[i].Eref;
        hfat += (*atoms)[i].Hf;
    }
    return hfat-eat*ev2kcal;
}

// from the routine gover.f
double mindo3::overlap(bfnc bfi, bfnc bfj) { 
    int i,j;
    double ri[3], rj[3];							// distance in bohr
    for ( i=0; i<3; i++ ) {
      ri[i]=(*atoms)[bfi.natom].xyz[i]/bohr2ang;
      rj[i]=(*atoms)[bfj.natom].xyz[i]/bohr2ang;
    }
    double RR = power(ri[0]-rj[0],2)+power(ri[1]-rj[1],2)+power(ri[2]-rj[2],2);
    int itype = bfi.type;
    int jtype = bfj.type;
    double Sij = 0.0;
    for ( i=0; i<stong; i++) {                                 
      for ( j=0; j<stong; j++) {                               
            double amb = bfi.expn[i] + bfj.expn[j];                               
            double apb = bfi.expn[i] * bfj.expn[j]; 
            double adb = apb/amb;                             
            double tomb, abn;
            if ((itype > 0) && (jtype > 0)) {			                //is = 4
                tomb = (ri[itype-1]-rj[itype-1])*(ri[jtype-1]-rj[jtype-1]);
                abn = -adb*tomb;
                if (itype == jtype) abn += 0.5;
                abn = 4*abn*sqrt(apb)/amb;
            } else if (itype > 0) {       			                //is = 3
                tomb = (ri[itype-1]-rj[itype-1]);
                abn = -2*tomb*bfj.expn[j] * sqrt(bfi.expn[i])/amb;
            } else if (jtype > 0) {                                             //is = 2
                tomb = (ri[jtype-1]-rj[jtype-1]);
                abn =  2*tomb*bfi.expn[i] * sqrt(bfj.expn[j])/amb;
            } else {                                                            //is = 1
                abn = 1.0;
            }    
            if (adb*RR < 90) {
                Sij += bfi.coef[i] * bfj.coef[j] * 
                       pow(2*sqrt(apb)/amb,1.5) * exp(-adb*RR)*abn;
            }
      }
    }
    return Sij;
}

// from the routine dcart.f
double mindo3::doverlap(struct bfnc bfi, struct bfnc bfj, int dir) {

    int i,j;
    double DS = 0.0, SS, amb, apb, adb, abn, adr, del1, del2 ,del3;
    double ri[3], rj[3];							// distance in bohr
    for ( i=0; i<3; i++ ) {
      ri[i]=(*atoms)[bfi.natom].xyz[i]/bohr2ang;
      rj[i]=(*atoms)[bfj.natom].xyz[i]/bohr2ang;
    }
    double RR = power(ri[0]-rj[0],2)+power(ri[1]-rj[1],2)+power(ri[2]-rj[2],2);
    del1 = ri[dir] - rj[dir];
    int itype = bfi.type;
    int jtype = bfj.type;

    for ( i=0; i<stong; i++) {                                 
      for ( j=0; j<stong; j++) {                               
            amb = bfi.expn[i] + bfj.expn[j];                               
            apb = bfi.expn[i] * bfj.expn[j]; 
            adb = apb/amb;
            adr = min(adb*RR,35.0);
            if ((itype == 0) && (jtype == 0)) {					// ss
                // is=1
                abn = -2.*adb*del1/bohr2ang;
            } else if ((itype == 0) && (jtype > 0)) {				// sp
                if ((jtype-1) == dir) 
                    //is = 3
                    abn = 2*adb/sqrt(bfj.expn[j])*(1-2*adb*del1*del1)/bohr2ang;
                else {
                    // is = 2
                    del2 = ri[jtype-1]-rj[jtype-1];
                    abn = -4*adb*adb*del1*del2/sqrt(bfj.expn[j])/bohr2ang;
                }
            } else if ((itype > 0) && (jtype == 0)) { 				//ps
                if ((itype-1) == dir) { 
                    // is = 5
                    abn = -2*adb/sqrt(bfj.expn[j])*(1-2*adb*del1*del1)/bohr2ang;
                } else {
                    // is = 4
                    del2 = ri[itype-1]-rj[itype-1];
                    abn = 4*adb*adb*del1*del2/sqrt(bfi.expn[i])/bohr2ang;
                }
            } else if (itype == jtype) { 
                if (dir == (itype-1)) {
                    // is = 9 (p|p)
                    abn=-8*adb*adb*del1/sqrt(apb)*(1.5-adb*del1*del1)/bohr2ang;
                } else {
                    // is = 8 (p'|p')
                    del2 = ri[jtype-1]-rj[jtype-1];
                    abn=-8*pow(adb,2)*del1/sqrt(apb)*(0.5-adb*del2*del2)/bohr2ang;
                }
            } else if ((dir != (itype-1)) && (dir != (jtype-1))) {
                // is = 7(p'|p")
                del2 = ri[itype-1] - rj[itype-1];
                del3 = ri[jtype-1] - rj[jtype-1];
                abn=8*pow(adb,3)*del1*del2*del3/sqrt(apb)/bohr2ang;
            } else {
                // is = 6 (p|p') or (p'|p)
                del2 = ri[itype+jtype-dir-2]-rj[itype+jtype-dir-2];
                abn=-4*adb*adb*del2/sqrt(apb)*(1-2*adb*del1*del1)/bohr2ang;
            }
            SS = power(2*sqrt(apb)/amb,1.5)*exp(-adr)*abn;
            DS += SS*bfi.coef[i]*bfj.coef[j];
      }
    }
    return DS;
}

SymmetricMatrix F0,F1,F2,F;

// Form the zero-iteration (density matrix independent) Fock matrix
void mindo3::calc_F0() {

  if ( !F0_initialized ) {

    F0.resize(bfncs.size()); F0 = 0;

    for ( uint i=1; i<=bfncs.size(); i++ ) {
       int iat = bfncs[i-1].natom;
       F0(i,i) = bfncs[i-1].u;
       int jato = -1;
       for ( uint j=1; j<=bfncs.size(); j++ ) {
          int jat = bfncs[j-1].natom;
          if ( iat != jat ) {
             if ( jat != jato ) { 
               F0(i,i) -= gamma( (*atoms)[iat], (*atoms)[jat] ) * (*atoms)[jat].Z; 
               jato = jat; 
             }
             double betaij = Bxy[ (*atoms)[iat].atno ][ (*atoms)[jat].atno ];	// Resonanace integral for coupling between different atoms
             double Sij    = overlap( bfncs[i-1], bfncs[j-1] );
             double IPij   = bfncs[i-1].ip+bfncs[j-1].ip;
             F0(i,j) = betaij*IPij*Sij;
          }
       }
    }

    F0_initialized = TRUE;
  }

}

// Average occupation density matrix
void mindo3::guess_D() {

  if ( !D_initialized ) {

    D.resize(nbf); D = 0;
    int ibf = 1;
    for ( vector <struct atom_mindo>::iterator at = atoms->begin(); at!=atoms->end(); at++ ) {
        int atno = (*at).atno;
        for (int i=0;i<(*at).nbf;i++) {
            if (atno == 1)
                D(ibf+i,ibf+i) = (*at).Z/1.0;
            else                 
                D(ibf+i,ibf+i) = (*at).Z/4.0;
        }
        ibf += (*at).nbf;
    }

  D_initialized = TRUE;

  }

}

// Coulomb-like term for orbitals on the same atom
double mindo3::g(bfnc bfi, bfnc bfj) {
    int i = bfi.type;
    int j = bfj.type;
//printf("%u - %u\n",i,j);
    if      ((i==0) && (j==0)) return (*atoms)[bfi.natom].gss;
    else if ((i==0) || (j==0)) return (*atoms)[bfi.natom].gsp;
    else if (i==j)             return (*atoms)[bfi.natom].gpp;
                               return (*atoms)[bfi.natom].gppp;
}

// Exchange-like term for orbitals on the same atom
double mindo3::h(bfnc bfi, bfnc bfj) {
    if ((bfi.type==0) || (bfj.type==0)) 
        return (*atoms)[bfi.natom].hsp;
    else       
        return (*atoms)[bfi.natom].hppp;
}              

// One-center corrections to the core fock matrix
void mindo3::calc_F1() {

    F1.resize(bfncs.size()); F1=0.0;

    for ( uint i=1; i<=bfncs.size(); i++ ) {
       F1(i,i) = 0.5 * g(bfncs[i-1],bfncs[i-1]) * D(i,i);
       for ( uint j=1; j<=bfncs.size(); j++ ) 
         if ( (bfncs[i-1].natom == bfncs[j-1].natom) && (i!=j) ) {
            double gij = g ( bfncs[i-1], bfncs[j-1] ),
                   hij = h ( bfncs[i-1], bfncs[j-1] );
            F1(i,i) += (gij - 0.5*hij) * D(j,j);
            if (i>j) F1(i,j) += 0.5 * (3*hij-gij) * D(i,j);
         }
    }

}

// Two-center corrections to the core fock matrix
void mindo3::calc_F2() {

    F2.resize(bfncs.size()); F2=0.0;

    for ( uint i=1; i<=bfncs.size(); i++ ) {
       for ( uint j=1; j<=bfncs.size(); j++ ) {
          if ( bfncs[i-1].natom != bfncs[j-1].natom ) {
             double gammaij = gamma( (*atoms)[ bfncs[i-1].natom ], (*atoms)[ bfncs[j-1].natom ] );
             F2(i,i) +=         gammaij * D(j,j);
             F2(i,j) += -0.25 * gammaij * D(i,j);
          }
       }
    }

}

// Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]
void mindo3::mkdens() {
   Matrix d;
   d << Orbs.columns(1,nclosed);
   D = d * d.t();
}

// SCF procedure for closed-shell molecules
double mindo3::SCF () { 
//puts("SCF:");
    calc_F0();
    guess_D();
    double Eel,Eold = 0;
    for ( SCFit=0; SCFit<MaxSCF; SCFit++ ) {
        calc_F1();
        calc_F2();
	F << F0 + F1 + F2;					// F = F0+F1+F2;
//cout << "F:" << F << endl;
        Eel = 0.5 * trace( D * ( F0 + F ) );		       	// Trace(D,*F0+F);
        Ediff = fabs(Eel-Eold);
//printf("%14.8lf %14.8lf %14.8le\n",Eold,Eel,Ediff);
        if ( Ediff < TolSCF ) break;
        Eold = Eel;
        EigenValues( F, Orbe, Orbs );				// or Jacobi(F, Dm, U), which is slower, but more reliable
        mkdens();
        D+=D;
    }

    return Eel;
}

double mindo3::Energy() {

  double Enuke	= enuke();
  double eref	= refeng();
  double Eel	= SCF();
  double Etot	= Eel+Enuke;
  double Hf	= Etot*ev2kcal+eref;

//printf("Enuke=%10.6lf, eref=%10.6lf, Eel=%10.6lf, Etot=%10.6lf, Hf=%10.6lf\n",Enuke,eref,Eel,Etot,Hf);

  return Hf;

}

void mindo3::print_orbs() {
  uint i;
  printf("%4s %4s "," "," ");
  for ( i=0; i<numbf();i++) printf("%11.6lf",Orbe(i+1)/*/ev2kcal*/); puts("");
  for ( i=0; i<bfncs.size(); i++ ) {
      printf("%4s %4s ",at_nam[(*atoms)[bfncs[i].natom].atno-1], or_nam[bfncs[i].type]);
      for(int j=0;j<nbf;j++)
        printf("%11.6lf",Orbs(i+1,j+1));
      puts("");
  }
}

void mindo3::print_grads(Matrix G) {
  printf("%16s%16s%16s\n","X","Y","Z");
  for(unsigned int i=0; i<atoms->size(); i++) {
     for(int j=0;j<3;j++)
        printf("%16.8le",G(i+1,j+1));
     puts("");
  }
}

void mindo3::print_hess(SymmetricMatrix H) {
//  printf("%16s%16s%16s\n","X","Y","Z");
  for(unsigned int i=0; i<atoms->size()*3; i++) {
    for(unsigned int j=0; j<atoms->size()*3; j++) {
        if (i>=j) 
           printf("%13.9lf",H(i+1,j+1)/2240.579577); 
        else 
           printf("%13s"," ");
     }
     puts("");
  }
}

// Gaussian type orbital, AO = GTO
double mindo3::gto(struct bfnc bf, double xyz[3]) {
  double R2 = dist2( (*atoms)[bf.natom].xyz, xyz );
  double sum = 0.0; 
  for ( int i=0; i<6; i++ ) {
    sum += bf.coef[i] * exp( -bf.expn[i] * R2 );
  }
  if (bf.type == 0)
     return sum;			// S
  else                
     return sum * xyz[ bf.type-1 ];	// Pxyz
}

// Molecular orbital, MO = LCAO = LCGTO
double mindo3::molorb(int numorb, double xyz[3]) {
  double sum = 0.0;
  for ( uint i=0; i<nbf; i++ ) {
      sum += Orbs(i+1,numorb+1) * gto(bfncs[i],xyz);
  }
  return sum;
}

double mindo3::MOEnergy(int imo) {
  return Orbe(imo+1);
}


// Compute analytic forces on list of atoms
void mindo3::forces(Matrix *forces) {
    int nat = atoms->size(), bfi,bfj;
    
    forces->resize(nat,3); *forces=0.0;

    // Loop over all pairs of atoms and compute the force between them
    for (int iat=0, iatb=0; iat<nat; iat++) {
        struct atom_mindo atomi = (*atoms)[iat];
        for (int jat=0, jatb=0; jat<iat; jat++ ) {
            struct atom_mindo atomj = (*atoms)[jat];
            double alpha = axy[atomi.atno][atomj.atno];
            double beta  = Bxy[atomi.atno][atomj.atno];
            double R2 = dist2(atomi.xyz,atomj.xyz);
            double R = sqrt(R2);
            double c2 = 0.25*power(atomi.rho+atomj.rho,2);

            for (int dir=0; dir<3; dir++ ){
                double Fij = 0.0;					// Force between atoms iat and jat in direction dir
                // initialize some constants
                double delta = atomi.xyz[dir]-atomj.xyz[dir];
                double c1  = delta*atomi.Z*atomj.Z*e2/R;
                double dr1 = e2*delta*power(R2+c2,-1.5);

                // Nuclear repulsion terms
                if ( ((atomi.atno == 1)
                      && ((atomj.atno == 7) || (atomj.atno == 8)))
                     || ((atomj.atno == 1)
                         && ((atomi.atno == 7) || (atomi.atno == 8))))
                    // Special case of NH or OH bonds
                    Fij += -c1*alpha*(1/R2 - R*power(R2+c2,-1.5)
                                      + 1/R - 1/sqrt(R2+c2))*exp(-R)
                                      - c1*R*pow(R2+c2,-1.5);
                else
                    Fij += -c1*(1/R2 - R*power(R2+c2,-1.5) + alpha/R 
                                - alpha/sqrt(R2+c2))*exp(-alpha*R)
                                - c1*R*power(R2+c2,-1.5);

                // Overlap terms
                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++) {
                    for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++) {
                        double Dij = D(bfi+1,bfj+1);
                        double dSij = doverlap( bfncs[bfi], bfncs[bfj], dir );
                        Fij += 2*beta*( bfncs[bfi].ip + bfncs[bfj].ip)*Dij*dSij;
                    }
                }

                // Core attraction terms
                for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++)
                    Fij += atomi.Z*D(bfj+1,bfj+1)*dr1;
                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++)
                    Fij += atomj.Z*D(bfi+1,bfi+1)*dr1;

                // Two-electron terms
                for ( bfi=iatb; bfi<iatb+atomi.nbf; bfi++) {
                   for ( bfj=jatb; bfj<jatb+atomj.nbf; bfj++) {
                        double Dii = D(bfi+1,bfi+1);
                        double Djj = D(bfj+1,bfj+1);
                        double Dij = D(bfi+1,bfj+1);
                        // exchange is the first term, coulomb is second:
                        Fij += 0.5*dr1*power(Dij,2)-dr1*Dii*Djj;
                   }
                }

                // Now sum total forces and convert to kcal/mol
                (*forces)(iat+1,dir+1) += ev2kcal*Fij;
                (*forces)(jat+1,dir+1) -= ev2kcal*Fij;
            }
            jatb += (*atoms)[jat].nbf;
        }
        iatb += (*atoms)[iat].nbf;
    }
}

// Compute numerical forces on list of atoms - Right side finite differencial
// dx = 1.0E-7 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_right(Matrix *forces, double dx) {
    int nat = atoms->size();
    forces->resize(nat,3); *forces=0.0;

    double fnc1 = Energy();

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] += dx;
		F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] -= dx;

                (*forces)(iat+1,dir+1) = (fnc2-fnc1)/dx;
            }
    }
}

// Compute numerical forces on list of atoms - Left side finite differencial
// dx = 1.0E-7 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_left(Matrix *forces, double dx) {
    int nat = atoms->size();
    forces->resize(nat,3); *forces=0.0;

    double fnc1 = Energy();

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] -= dx;
		F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] += dx;

                (*forces)(iat+1,dir+1) = (fnc1-fnc2)/dx;
            }
    }
}

// Compute numerical forces on list of atoms - Central finite differencial
// dx = 1.0E-4 ... 1.0E-8 is the best choice according to my experiments
void mindo3::num_forces_central(Matrix *forces, double dx) {

    forces->resize(nat,3); *forces=0.0;

    for (int iat=0; iat<nat; iat++) {
            for (int dir=0; dir<3; dir++ ){
                
                (*atoms)[iat].xyz[dir] += dx;
		F0_initialized = FALSE;
                double fnc2 = Energy();
                (*atoms)[iat].xyz[dir] -= dx;

                (*atoms)[iat].xyz[dir] -= dx;
		F0_initialized = FALSE;
                double fnc1 = Energy();
                (*atoms)[iat].xyz[dir] += dx;

                (*forces)(iat+1,dir+1) = (fnc2-fnc1)/(dx+dx);
            }
    }
}

// Compute numerical forces on list of atoms - Right side finite differencial
void mindo3::num_hessian_right(SymmetricMatrix *hess, double dx) {
    Matrix fr1, fr2;

    int nat = atoms->size(), i, j;
    hess->resize(nat*3); *hess=0.0;

    Energy();
    forces( &fr1 );

    init();

    for (i=0; i<nat*3; i++) {
        int iat=i/3, dir=i%3;
        (*atoms)[iat].xyz[dir] += dx;
        Energy();
        forces(&fr2);
        (*atoms)[iat].xyz[dir] -= dx;
        for (j=0; j<=i; j++) {
           int jat=j/3, djr=j%3;
           (*hess)(i+1,j+1) = (fr2(jat+1,djr+1)-fr1(jat+1,djr+1))/dx;
        }
    }
}

// Compute numerical forces on list of atoms - Both side finite differencial
void mindo3::num_hessian_central(SymmetricMatrix *hess, double dx) {
    Matrix fr1, fr2;

    int nat = atoms->size(), i, j;
    hess->resize(nat*3); *hess=0.0;

    init();

    for (i=0; i<nat*3; i++) {
        int iat=i/3, dir=i%3;
        
        (*atoms)[iat].xyz[dir] += dx;
        Energy();
        forces(&fr2);
        (*atoms)[iat].xyz[dir] -= dx;

        (*atoms)[iat].xyz[dir] -= dx;
        Energy();
        forces(&fr1);
        (*atoms)[iat].xyz[dir] += dx;

        for (j=0; j<=i; j++) {
           int jat=j/3, djr=j%3;
           (*hess)(i+1,j+1) = (fr2(jat+1,djr+1)-fr1(jat+1,djr+1))/(dx+dx);
        }
    }
}

