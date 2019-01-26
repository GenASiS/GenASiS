#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "nuc_eos.h"

int main(void) {


  nuc_eos_C_ReadTable("eostable.h5");

  printf("nrho: %d\n",nrho);
  printf("ntemp: %d\n",ntemp);
  printf("nye: %d\n",nye);

  // test EOS call

#if 0
  double xrho = 1.0e12;
  double lr = log10(xrho);
  double xtemp = 2.0e0;
  double lt = log10(xtemp);
  double xye = 0.4e0;

  double res[2];
  int ivs[2];

  ivs[0] = 0;
  ivs[1] = 1;

  nuc_eos_C_linterp_some(lr,lt,xye,res,
			 alltables, ivs, nrho, ntemp, nye, 2,
			 logrho, logtemp, yes);

  double lpress = res[0];
  double leps = res[1];
  
  printf("r,t,ye: %15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  printf("press eps: %15.6E %15.6E\n",lpress,leps);
  

  double nlt;
  int keyerrt = 0;
  leps = log10(pow(10.0,leps) * 1.01);
  lt = lt;
  nuc_eos_C_findtemp(lr,lt,xye,leps,1.0e-10,&nlt,&keyerrt);
  printf("t tnew: %15.6E %15.6E\n",lt,nlt);
  printf("keyerr: %d\n",keyerrt);


  ivs[0] = 0;
  ivs[1] = 1;

  nuc_eos_C_linterp_some(lr,nlt,xye,res,
			 alltables, ivs, nrho, ntemp, nye, 2,
			 logrho, logtemp, yes);

  lpress = res[0];
  leps = res[1];
  printf("lr,lt,ye,leps: %15.6E %15.6E %15.6E %15.6E\n",lr,nlt,xye,leps);
  printf("press eps: %15.6E %15.6E\n",lpress,leps);

#endif
  nuc_eos_C_testing();



  double xrho = 1.0e12;
  double xtemp = 2.0;
  double xye = 0.35;

  // outvars
  double xeps,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu;
  const double prec = 1.0e-10;
  int keytemp,keyerr;

  keytemp = 1;
  nuc_eos_C_short(xrho,&xtemp,xye,&xeps,&xprs,&xent,&xcs2,
		  &xdedt,&xdpderho,&xdpdrhoe,&xmunu,keytemp,
		  &keyerr,prec);

  printf("keytemp: %d\n",keytemp);
  printf("keyerr: %d\n",keyerr);

  printf("%15.6E %15.6E %15.6E %15.6E\n",xrho,xtemp,xye,xeps);
  printf("press: %15.6E     ent: %15.6E\n",xprs,xent);
  printf("  cs2: %15.6E    dedt: %15.6E\n",xcs2,xdedt);

  printf("*********\n");
  double xtemp0 = xtemp;
  xtemp = 10.0;
  keytemp = 0;
  nuc_eos_C_short(xrho,&xtemp,xye,&xeps,&xprs,&xent,&xcs2,
		  &xdedt,&xdpderho,&xdpdrhoe,&xmunu,keytemp,
		  &keyerr,prec);

  printf("keytemp: %d\n",keytemp);
  printf("keyerr: %d\n",keyerr);

  printf("%15.6E %15.6E %15.6E %15.6E %15.6E\n",xrho,xtemp,xtemp0,xye,xeps);
  printf("press: %15.6E     ent: %15.6E\n",xprs,xent);
  printf("  cs2: %15.6E    dedt: %15.6E\n",xcs2,xdedt);



  return 0;
}

