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
    
    
    SUBROUTINE intp3d ( x, y, z, f, kt, ft, nx, ny, nz, i_ft, xt, yt, zt, &
                          d1, d2, d3 )     
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
!     d1       centered derivative of ft with respect to x
!     d2       centered derivative of ft with respect to y
!     d3       centered derivative of ft with respect to z
!     Note that d? only make sense when intp3d is called with kt=1
!---------------------------------------------------------------------
!
! 

!
      integer, intent ( in ) :: &
        kt, nx, ny, nz, i_ft
      real*8, dimension ( : ), intent ( in ) :: &
        x, y, z, xt, yt, zt
      real*8, dimension ( : ), intent ( out ) :: &
        f
      real*8, dimension ( :, :, :, : ), intent ( in ) :: &
        ft
      real*8, intent ( out ) :: &
        d1,d2,d3
!
!
      integer, PARAMETER :: ktx = 1
      real*8 :: &
        fh(ktx,8), delx(ktx), dely(ktx), delz(ktx), &
        a1(ktx), a2(ktx), a3(ktx), a4(ktx), &
        a5(ktx), a6(ktx), a7(ktx), a8(ktx)

      double precision dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz
#ifdef ENABLE_OMP_OFFLOAD
      !$OMP declare target
#endif

      !IF (kt .GT. ktx)  STOP'***KTX**'
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
      DO  n = 1, kt                                            
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
         fh(n,1) = ft(ix  , iy  , iz  , i_ft)                             
         fh(n,2) = ft(ix-1, iy  , iz  , i_ft)                             
         fh(n,3) = ft(ix  , iy-1, iz  , i_ft)                             
         fh(n,4) = ft(ix  , iy  , iz-1, i_ft)                             
         fh(n,5) = ft(ix-1, iy-1, iz  , i_ft)                             
         fh(n,6) = ft(ix-1, iy  , iz-1, i_ft)                             
         fh(n,7) = ft(ix  , iy-1, iz-1, i_ft)                             
         fh(n,8) = ft(ix-1, iy-1, iz-1, i_ft)                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
!                                                    
         a1(n) = fh(n,1)                             
         a2(n) = dxi   * ( fh(n,2) - fh(n,1) )       
         a3(n) = dyi   * ( fh(n,3) - fh(n,1) )       
         a4(n) = dzi   * ( fh(n,4) - fh(n,1) )       
         a5(n) = dxyi  * ( fh(n,5) - fh(n,2) - fh(n,3) + fh(n,1) )
         a6(n) = dxzi  * ( fh(n,6) - fh(n,2) - fh(n,4) + fh(n,1) )
         a7(n) = dyzi  * ( fh(n,7) - fh(n,3) - fh(n,4) + fh(n,1) )
         a8(n) = dxyzi * ( fh(n,8) - fh(n,1) + fh(n,2) + fh(n,3) + &
                           fh(n,4) - fh(n,5) - fh(n,6) - fh(n,7) )
!
         d1 = -a2(n)
         d2 = -a3(n)
         d3 = -a4(n)
         f(n)  = a1(n) +  a2(n) * delx(n) &
                       +  a3(n) * dely(n) &
                       +  a4(n) * delz(n) &
                       +  a5(n) * delx(n) * dely(n) &
                       +  a6(n) * delx(n) * delz(n) &
                       +  a7(n) * dely(n) * delz(n) &
                       +  a8(n) * delx(n) * dely(n) * delz(n)     
!
      ENDDO
!                                                                    

      END SUBROUTINE intp3d

end module linterp
