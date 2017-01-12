#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.h"

void nuc_eos_C_short(double xrho, double *xtemp, double xye,
		     double *xenr, double* xprs, double* xent,
		     double *xcs2, double* xdedt, double* xdpderho,
		     double *xdpdrhoe, double* xmunu, int keytemp,
		     int *keyerr,double rfeps) 

{
  // This routine expects rho in g/cm^3 and T in MeV.
  // It will strictly return values in cgs or MeV

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)

  *keyerr = 0;

  if(xrho > eos_rhomax) {
    *keyerr = 105;
    return;
  }
  if(xrho < eos_rhomin) {
    *keyerr = 106;
    return;
  }
  if(xye > eos_yemax) {
    *keyerr = 101;
    return;
  }
  if(xye < eos_yemin) {
    *keyerr = 102;
    return;
  }
  if(keytemp == 1) {
    if(*xtemp > eos_tempmax) {
      *keyerr = 103;
      return;
    }
    if(*xtemp < eos_tempmin) {
      *keyerr = 104;
      return;
    }
  }

  // set up local vars
  double lr = log10(xrho);
  double lt = log10(*xtemp);
  double xeps = *xenr + energy_shift;
  double leps = log10(MAX(xeps,1.0));



  // find temperature if need be
  if(keytemp == 0) {
    double nlt = 0.0;
    nuc_eos_C_findtemp(lr,lt,xye,leps,rfeps,&nlt,keyerr);
    if(*keyerr != 0) return;
    lt = nlt;
    *xtemp = pow(10.0,lt);
  } else if(keytemp == 2) {
    fprintf(stderr,"keytemp = 2 not yet implemented. Aborting!\n");
    abort();
  }

  double res[8];
  nuc_eos_C_linterp_some(lr,lt,xye,res,alltables,ivs_short,
			 nrho,ntemp,nye,8,logrho,logtemp,yes);
  
  // assign results

  if(keytemp != 0) {
    // set up internal energy correctly, correcting for shift
    *xenr = pow(10.0,res[1]) - energy_shift;
  }

  if(keytemp != 2) {
    // reset entropy
    *xent = res[2];
  }

  *xprs     = pow(10.0,res[0]);
  *xmunu    = res[3];
  *xcs2     = res[4];
  *xdedt    = res[5];
  *xdpdrhoe = res[6];
  *xdpderho = res[7];

  return;
}
