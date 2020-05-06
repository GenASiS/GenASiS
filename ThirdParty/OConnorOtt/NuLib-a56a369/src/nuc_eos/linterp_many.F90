module linterp

  implicit none

contains

      SUBROUTINE intp3d_many ( x, y, z, f, kt, ft, nx, ny, nz, nvars, xt, yt, zt)
!
!                                                          
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula          
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------


      integer, intent ( in ) :: &
         kt,nx,ny,nz,nvars
      real*8, intent ( in ), dimension ( :, :, :, : ) :: &
        ft

      real*8, intent ( in ), dimension ( : ) :: &
         x, y, z, xt, yt, zt
      real*8, intent ( out ), dimension ( :, : ) :: &
        f

#ifdef ENABLE_OMP_OFFLOAD
      !$OMP declare target
#endif
      
      integer :: iv      
      real*8 d1,d2,d3
!
!
!      integer,parameter :: ktx = 1

      real*8  fh(kt,8,nvars), delx(kt), dely(kt), delz(kt), &
           a1(kt,nvars), a2(kt,nvars), a3(kt,nvars), a4(kt,nvars), &
           a5(kt,nvars), a6(kt,nvars), a7(kt,nvars), a8(kt,nvars)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      !-- FIXME: commented this out since STOP is not supported on device 
      !IF (kt .GT. ktx)  STOP '***KTX**'
!
!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / FLOAT(nx-1)
      dy    = (yt(ny) - yt(1)) / FLOAT(ny-1)
      dz    = (zt(nz) - zt(1)) / FLOAT(nz-1)
!
      dxi   = 1. / dx
      dyi   = 1. / dy
      dzi   = 1. / dz
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!
      dxyzi = dxi * dyi * dzi
!
!
!------- loop over all points to be interpolated
!
      do  n = 1, kt                                            
!
!------- determine location in (equidistant!!!) table 
!                                                                  
         ix = 2 + INT( (x(n) - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y(n) - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z(n) - zt(1) - 1.e-10) * dzi )
!                                                     
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
!
!         write(*,*) iy-1,iy,iy+1
!
!------- set-up auxiliary arrays for Lagrange interpolation
!                                                                 
         delx(n) = xt(ix) - x(n)
         dely(n) = yt(iy) - y(n)
         delz(n) = zt(iz) - z(n)
!      
         do iv = 1, nvars
            fh(n,1,iv) = ft(ix  , iy  , iz, iv  )                             
            fh(n,2,iv) = ft(ix-1, iy  , iz, iv  )                             
            fh(n,3,iv) = ft(ix  , iy-1, iz, iv  )                             
            fh(n,4,iv) = ft(ix  , iy  , iz-1, iv)                             
            fh(n,5,iv) = ft(ix-1, iy-1, iz, iv  )                             
            fh(n,6,iv) = ft(ix-1, iy  , iz-1, iv)                             
            fh(n,7,iv) = ft(ix  , iy-1, iz-1, iv)                             
            fh(n,8,iv) = ft(ix-1, iy-1, iz-1, iv)                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
            !                                                    
            a1(n,iv) = fh(n,1,iv)                             
            a2(n,iv) = dxi   * ( fh(n,2,iv) - fh(n,1,iv) )       
            a3(n,iv) = dyi   * ( fh(n,3,iv) - fh(n,1,iv) )       
            a4(n,iv) = dzi   * ( fh(n,4,iv) - fh(n,1,iv) )       
            a5(n,iv) = dxyi  * ( fh(n,5,iv) - fh(n,2,iv) - fh(n,3,iv) + fh(n,1,iv) )
            a6(n,iv) = dxzi  * ( fh(n,6,iv) - fh(n,2,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a7(n,iv) = dyzi  * ( fh(n,7,iv) - fh(n,3,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a8(n,iv) = dxyzi * ( fh(n,8,iv) - fh(n,1,iv) + fh(n,2,iv) + fh(n,3,iv) + &
                 fh(n,4,iv) - fh(n,5,iv) - fh(n,6,iv) - fh(n,7,iv) )
!
            f(n,iv)  = a1(n,iv) +  a2(n,iv) * delx(n)                         &
                 +  a3(n,iv) * dely(n)                         &
                 +  a4(n,iv) * delz(n)                         &
                 +  a5(n,iv) * delx(n) * dely(n)               &
                 +  a6(n,iv) * delx(n) * delz(n)               &
                 +  a7(n,iv) * dely(n) * delz(n)               &
                 +  a8(n,iv) * delx(n) * dely(n) * delz(n)     
!
         enddo

      enddo
!                                                                    
      
    end SUBROUTINE intp3d_many
 
end module linterp