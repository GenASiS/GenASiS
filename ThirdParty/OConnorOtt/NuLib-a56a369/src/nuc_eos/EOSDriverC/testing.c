#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "nuc_eos.h"


void nuc_eos_C_testing() {


  int npoints = 5000000;

  // Randomly draw npoints from
  // range of logrho, logT, Y_e
  double lrhomin = log10(pow(10.0,logrho[0])*1.01);
  double lrhomax = log10(pow(10.0,logrho[nrho-1])*0.99);
  double ltempmin = log10(pow(10.0,logtemp[0])*1.01);
  double ltempmax = log10(pow(10.0,logtemp[ntemp-1])*0.99);
  double yemin = yes[0]*1.01;
  double yemax = yes[nye-1]*0.99;

  double lr_range = lrhomax-lrhomin;
  double lt_range = ltempmax-ltempmin;
  double ye_range = yemax-yemin;


  int failcount1 = 0;
  int failcount2 = 0;
  for(int i=0;i<npoints;i++) {

    double rand1 = rand() / (double)RAND_MAX;
    double rand2 = rand() / (double)RAND_MAX;
    double rand3 = rand() / (double)RAND_MAX;
    double rand4 = rand() / (double)RAND_MAX;
    
    double xlr = lrhomin + lr_range*rand1;
    double xlt = ltempmin + lt_range*rand2;
    double xye = yemin + ye_range*rand3;

    // get eps
    double xeps;
    int ivs[1] = {1};
    int nvars = 1;

    nuc_eos_C_linterp_some(xlr, xlt, xye, &xeps, alltables, 
			   ivs, nrho, ntemp, nye, nvars,
			   logrho, logtemp, yes);

    
    
    // now set temperature to something crazy and try to
    // recover the correct temperature
    double savedtemp = xlt;
    double nlt;
    int keyerrt = 0;
    xlt = log10(pow(10.0,xlt)+ pow(10.0,xlt)*(-1+2.0*rand4));

    nuc_eos_C_findtemp(xlr,xlt,xye,xeps,1.0e-10,&nlt,&keyerrt);
    if(keyerrt != 0) {
      failcount1++;
      //      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n", i, xlr, xlt, xye, xeps);
    } else {
      if(abs(nlt-savedtemp) > 1.0e-10) {
	failcount2++;
      }
      //      assert( abs(nlt-savedtemp) <= 1.0e-10);
    }


  }
  fprintf(stdout,"failcount1: %d fraction: %15.6E\n",failcount1,failcount1/(1.0*npoints));
  fprintf(stdout,"failcount2: %d fraction: %15.6E\n",failcount2,failcount2/(1.0*npoints));
  
  
  return;
}

