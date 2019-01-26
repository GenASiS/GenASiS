#include<stdio.h>
#include<assert.h>
#include "nuc_eos.h"


void nuc_eos_C_linterp_for_temp(double x, double y, double z,
				double* f, double* ft, 
				int nx, int ny, int nz, 
				double* xt,double*yt, double* zt, 
				double* dlepsdlt)

{

//!
//!     purpose: interpolation of a function of three variables in an
//!              equidistant(!!!) table.
//!
//!     method:  8-point Lagrange linear interpolation formula
//!
//!     x        input vector of first  variable
//!     y        input vector of second variable
//!     z        input vector of third  variable
//!
//!     f        output vector of interpolated function values
//!
//!     ft       3d array of tabulated function values
//!     nx       x-dimension of table
//!     ny       y-dimension of table
//!     nz       z-dimension of table
//!     xt       vector of x-coordinates of table
//!     yt       vector of y-coordinates of table
//!     zt       vector of z-coordinates of table
//!
//!---------------------------------------------------------------------


  // helper variables
  double fh[8], delx, dely, delz, a[8];
  double dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi;
  int n,ix,iy,iz;
  
  // determine spacing parameters of equidistant (!!!) table
  dx = (xt[nx-1] - xt[0]) / (1.0*(nx-1));
  dy = (yt[ny-1] - yt[0]) / (1.0*(ny-1));
  dz = (zt[nz-1] - zt[0]) / (1.0*(nz-1));
  
  dxi = 1.0 / dx;
  dyi = 1.0 / dy;
  dzi = 1.0 / dz;
  
  dxyi = dxi * dyi;
  dxzi = dxi * dzi;
  dyzi = dyi * dzi;

  dxyzi = dxi * dyi * dzi;

  // determine location in table

  ix = 1 + (int)( (x - xt[0] - 1.0e-10) * dxi );
  iy = 1 + (int)( (y - yt[0] - 1.0e-10) * dyi );
  iz = 1 + (int)( (z - zt[0] - 1.0e-10) * dzi );

  ix = MAX( 1, MIN( ix, nx-1 ) );
  iy = MAX( 1, MIN( iy, ny-1 ) );
  iz = MAX( 1, MIN( iz, nz-1 ) );

  // set up aux vars for interpolation

  delx = xt[ix] - x;
  dely = yt[iy] - y;
  delz = zt[iz] - z;

  int idx[8];

  idx[0] = NTABLES*(ix + nx*(iy + ny*iz));
  idx[1] = NTABLES*((ix-1) + nx*(iy + ny*iz));
  idx[2] = NTABLES*(ix + nx*((iy-1) + ny*iz));
  idx[3] = NTABLES*(ix + nx*(iy + ny*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nx*((iy-1) + ny*iz));
  idx[5] = NTABLES*((ix-1) + nx*(iy + ny*(iz-1)));
  idx[6] = NTABLES*(ix + nx*((iy-1) + ny*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nx*((iy-1) + ny*(iz-1)));

  
  int iv = 1; // log(eps)
    
  // set up aux vars for interpolation
  // assuming array ordering (iv, ix, iy, iz)

  fh[0] = ft[iv+idx[0]];
  fh[1] = ft[iv+idx[1]];
  fh[2] = ft[iv+idx[2]];
  fh[3] = ft[iv+idx[3]];
  fh[4] = ft[iv+idx[4]];
  fh[5] = ft[iv+idx[5]];
  fh[6] = ft[iv+idx[6]];
  fh[7] = ft[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values

  a[0] = fh[0];
  a[1] = dxi *   ( fh[1] - fh[0] );
  a[2] = dyi *   ( fh[2] - fh[0] );
  a[3] = dzi *   ( fh[3] - fh[0] );
  a[4] = dxyi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = dxzi *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dyzi *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = dxyzi * ( fh[7] - fh[0] + fh[1] + fh[2] + 
		   fh[3] - fh[4] - fh[5] - fh[6] );

  *dlepsdlt = -a[2];

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;


  return;
}
  
