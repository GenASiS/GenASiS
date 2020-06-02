#include "Preprocessor"

submodule ( EOS_P_HN_OConnorOtt__Form ) EOS_P_HN_OConnorOtt__Kernel

  use Basics
  
  implicit none

contains

  
  module procedure Interpolate_3D_Kernel 
    
    integer ( KDI ) :: &
      iValue, &
      iV, &  !-- iVariable
      iS, &  !-- iSelected
      iF, &  !-- iFluid
      nx,ny,nz,nvars, &
      nValues, &
      ix,iy,iz
    real ( KDR ) :: &
      L_x, &
      L_y, &
      delx, dely, delz, &
      dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi,&
      a1, a2, a3, a4, a5, a6, a7, a8 
    real ( KDR ), dimension ( 8 ) :: &
      fh
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nx = size ( XT )
    ny = size ( YT )
    nz = size ( ZT )
    nvars = size ( T, dim = 4 )
    nValues = size ( F, dim = 1 )
    
    !--  determine spacing parameters of (equidistant!!!) table
    dx    = (xt(nx) - xt(1)) / real(nx-1, kind=KDR)
    dy    = (yt(ny) - yt(1)) / real(ny-1, kind=KDR)
    dz    = (zt(nz) - zt(1)) / real(nz-1, kind=KDR)

    dxi   = 1. / dx
    dyi   = 1. / dy
    dzi   = 1. / dz

    dxyi  = dxi * dyi
    dxzi  = dxi * dzi
    dyzi  = dyi * dzi

    dxyzi = dxi * dyi * dzi
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( L_x, L_y, ix, iy, iz, &
      !$OMP&           delx, dely, delz, iV, iF, iS, fh, &
      !$OMP&           a1, a2, a3, a4, a5, a6, a7, a8 ) &
      !$OMP& shared ( nx, ny, nz, nvars, nvalues, dx, dy, dz, &
      !$OMP&          dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi, E_Shift )
      do  iValue = 1, nValues
        do iS = 1, size ( ia_E )

          !-- Convert to log space for the x and y
          L_x = log10 ( F ( iValue, ia_F_I ( 1 ) ) )
          L_y = log10 ( F ( iValue, ia_F_I ( 2 ) ) )
          
          !-- determine location in (equidistant!!!) table 
          ix = 2 + INT( (L_x - xt(1) - 1.e-10_KDR) * dxi )
          iy = 2 + INT( (L_y - yt(1) - 1.e-10_KDR) * dyi )
          iz = 2 + INT( ( F ( iValue, ia_F_I ( 3 ) ) - zt(1) - 1.e-10_KDR) * dzi )
                                                     
          ix = MAX( 2, MIN( ix, nx ) )
          iy = MAX( 2, MIN( iy, ny ) )   
          iz = MAX( 2, MIN( iz, nz ) )

          !-- set-up auxiliary arrays for Lagrange interpolation
                                                                 
          delx = xt(ix) - L_x
          dely = yt(iy) - L_y
          delz = zt(iz) - F ( iValue, ia_F_I ( 3 ) )
          
          iV = ia_E ( iS )
          iF = ia_F_O ( iS )
          
          fh(1) = T(ix  , iy  , iz,   iv)
          fh(2) = T(ix-1, iy  , iz,   iv)   
          fh(3) = T(ix  , iy-1, iz,   iv)   
          fh(4) = T(ix  , iy  , iz-1, iv)
          fh(5) = T(ix-1, iy-1, iz,   iv)
          fh(6) = T(ix-1, iy  , iz-1, iv)
          fh(7) = T(ix  , iy-1, iz-1, iv)
          fh(8) = T(ix-1, iy-1, iz-1, iv)
          
          !-- set up coefficients of the interpolation polynomial and 
          !   evaluate function values 
          a1 = fh(1)
          a2 = dxi   * ( fh(2) - fh(1) )
          a3 = dyi   * ( fh(3) - fh(1) )
          a4 = dzi   * ( fh(4) - fh(1) )
          a5 = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
          a6 = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
          a7 = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
          a8 = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
               fh(4) - fh(5) - fh(6) - fh(7) )

          f(iValue, iF)  &
            = a1 +  a2 * delx               &
               +  a3 * dely                      &  
               +  a4 * delz                      &  
               +  a5 * delx * dely               &  
               +  a6 * delx * delz               &
               +  a7 * dely * delz               &
               +  a8 * delx * dely * delz
          
          !-- T ( :, :, :, 1 ) and T ( :, :, :, 2 ) is in log space
          if ( iV == 1 .or. iV == 2 ) then
            f ( iValue, iF ) = 10.0e0_KDR ** f ( iValue, iF )
          end if
          if ( iV == 2 ) then
            f ( iValue, iF ) = f ( iValue, iF ) - E_Shift
          end if
                    
        enddo
        
      enddo
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP  parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( L_x, L_y, ix, iy, iz, &
      !$OMP&           delx, dely, delz, iV, iF, iS, fh, &
      !$OMP&           a1, a2, a3, a4, a5, a6, a7, a8 ) &
      !$OMP& shared ( nx, ny, nz, nvars, nvalues, dx, dy, dz, &
      !$OMP&          dxi, dyi, dzi, dxyi, dxzi, dyzi, dxyzi, E_Shift )
      do  iValue = 1, nValues
        do iS = 1, size ( ia_E )

          !-- Convert to log space for the x and y
          L_x = log10 ( F ( iValue, ia_F_I ( 1 ) ) )
          L_y = log10 ( F ( iValue, ia_F_I ( 2 ) ) )
          
          !-- determine location in (equidistant!!!) table 
          ix = 2 + INT( (L_x - xt(1) - 1.e-10_KDR) * dxi )
          iy = 2 + INT( (L_y - yt(1) - 1.e-10_KDR) * dyi )
          iz = 2 + INT( ( F ( iValue, ia_F_I ( 3 ) ) - zt(1) - 1.e-10_KDR) * dzi )
                                                     
          ix = MAX( 2, MIN( ix, nx ) )
          iy = MAX( 2, MIN( iy, ny ) )   
          iz = MAX( 2, MIN( iz, nz ) )

          !-- set-up auxiliary arrays for Lagrange interpolation
                                                                 
          delx = xt(ix) - L_x
          dely = yt(iy) - L_y
          delz = zt(iz) - F ( iValue, ia_F_I ( 3 ) )
          
          iV = ia_E ( iS )
          iF = ia_F_O ( iS )
          
          fh(1) = T(ix  , iy  , iz,   iv)
          fh(2) = T(ix-1, iy  , iz,   iv)   
          fh(3) = T(ix  , iy-1, iz,   iv)   
          fh(4) = T(ix  , iy  , iz-1, iv)
          fh(5) = T(ix-1, iy-1, iz,   iv)
          fh(6) = T(ix-1, iy  , iz-1, iv)
          fh(7) = T(ix  , iy-1, iz-1, iv)
          fh(8) = T(ix-1, iy-1, iz-1, iv)
          
          !-- set up coefficients of the interpolation polynomial and 
          !   evaluate function values 
          a1 = fh(1)
          a2 = dxi   * ( fh(2) - fh(1) )
          a3 = dyi   * ( fh(3) - fh(1) )
          a4 = dzi   * ( fh(4) - fh(1) )
          a5 = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
          a6 = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
          a7 = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
          a8 = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
               fh(4) - fh(5) - fh(6) - fh(7) )

          f(iValue, iF)  &
            = a1 +  a2 * delx               &
               +  a3 * dely                      &  
               +  a4 * delz                      &  
               +  a5 * delx * dely               &  
               +  a6 * delx * delz               &
               +  a7 * dely * delz               &
               +  a8 * delx * dely * delz
          
          !-- T ( :, :, :, 1 ) and T ( :, :, :, 2 ) is in log space
          if ( iV == 1 .or. iV == 2 ) then
            f ( iValue, iF ) = 10.0e0_KDR ** f ( iValue, iF )
          end if
          if ( iV == 2 ) then
            f ( iValue, iF ) = f ( iValue, iF ) - E_Shift
          end if
                    
        enddo
        
      enddo
      !$OMP  end parallel do
      
    end if !-- Use device
                      
  end procedure Interpolate_3D_Kernel
  
  
  module procedure FindTemperatureKernel 
                 
    integer ( KDI ) :: &
      iV, &
      iI, &   !-- iIteration
      nValues, &
      nIterations
    real ( KDR ) :: &
      Shift, &
      Tolerance, &
      T_L_T_Max, T_L_T_Min, &
      L_N, Ye, &
      SV_R, D2, &
      L_T, L_T_1, &   !-- LogTemperature (from Fluid input)
      SV_0, SV_1, &
      L_DT
    logical ( KDL ) :: &
      LogScale, &
      UseDevice
    
    Shift = 0.0_KDR
    if ( present ( ShiftOption ) ) &
      Shift = ShiftOption
    
    nIterations = 20
    if ( present ( nIterationsOption ) ) &
      nIterations = nIterationsOption
    
    Tolerance = 1e-10_KDR
    if ( present ( ToleranceOption ) ) &
      Tolerance = ToleranceOption
      
    LogScale = .false.
    if ( present ( LogScaleOption ) ) &
      LogScale = .true.
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = .true.
    
    nValues = size ( F, dim = 1 )
    
    T_L_T_Max = T_L_T ( size ( T_L_T ) )
    T_L_T_Min = T_L_T ( 1 )
    
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& shared ( Shift, Tolerance, T_L_T_Max, T_L_T_Min, LogScale ), &
      !$OMP& private ( L_N, Ye, SV_R, D2, L_T, L_T_1, SV_0, SV_1, L_DT )
      Value_Device: do iV = 1, nValues
        
        L_N   = log10 ( F ( iV, ia_F_I ( 1 ) ) )
        L_T   = log10 ( F ( iV, ia_F_I ( 2 ) ) )
        Ye    = F ( iV, ia_F_I ( 3 )  )
        
        L_T_1 = L_T
        
        if ( LogScale ) then
          SV_0 = log10 ( max ( ( F ( iV, i_SF ) + Shift ), 1.0_KDR ) )
        else
          SV_0 = F ( iV, i_SF ) + Shift
        end if 
        SV_1 = SV_0
        
        call InterpolateTableKernel &
               ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
        
        if ( abs ( SV_R - SV_0 ) < Tolerance * abs ( SV_0 ) ) then
          cycle Value_Device
        end if
        
        do iI = 1, nIterations
          
          L_DT  = - ( SV_R - SV_0 ) / D2
          L_T_1 = L_T
          L_T   = max ( min ( ( L_T + L_DT ), T_L_T_Max ), T_L_T_Min )
          SV_1  = SV_R
          
          call InterpolateTableKernel &
                 ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
        
          if ( abs ( SV_R - SV_0 )  <  Tolerance * abs ( SV_0 ) ) then
            F ( iV, ia_F_I ( 2 ) ) = 10.0_KDR ** L_T
            cycle Value_Device
          endif 
          
          ! if we are closer than 10^-2  to the 
          ! root (eps-eps0)=0, we are switching to 
          ! the secant method, since the table is rather coarse and the
          ! derivatives may be garbage.

          if ( abs ( SV_R - SV_0 )  <  1e-3_KDR * abs ( SV_0 ) ) then
            D2 = ( SV_R - SV_1 ) / ( L_T - L_T_1 )
          end if
          
          !-- FIXME: Need error handling
          !if ( iI == nIterations ) &
          !  call Show ( 'Error, Max iteration reached' )
      
        end do 

      end do Value_Device
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& shared ( Shift, Tolerance, T_L_T_Max, T_L_T_Min, LogScale ), &
      !$OMP& private ( L_N, Ye, SV_R, D2, L_T, L_T_1, SV_0, SV_1, L_DT )
      Value_Host: do iV = 1, nValues
        
        L_N   = log10 ( F ( iV, ia_F_I ( 1 ) ) )
        L_T   = log10 ( F ( iV, ia_F_I ( 2 ) ) )
        Ye    = F ( iV, ia_F_I ( 3 )  )
        
        L_T_1 = L_T
        
        if ( LogScale ) then
          SV_0 = log10 ( max ( ( F ( iV, i_SF ) + Shift ), 1.0_KDR ) )
        else
          SV_0 = F ( iV, i_SF ) + Shift
        end if 
        SV_1 = SV_0
        
        call InterpolateTableKernel &
               ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
        
        if ( abs ( SV_R - SV_0 ) < Tolerance * abs ( SV_0 ) ) then
          cycle Value_Host
        end if
        
        do iI = 1, nIterations
          
          L_DT  = - ( SV_R - SV_0 ) / D2
          L_T_1 = L_T
          L_T   = max ( min ( ( L_T + L_DT ), T_L_T_Max ), T_L_T_Min )
          SV_1  = SV_R
          
          call InterpolateTableKernel &
                 ( L_N, L_T, Ye, T, T_L_N, T_L_T, T_Ye, i_ST, SV_R, D2 )
        
          if ( abs ( SV_R - SV_0 )  <  Tolerance * abs ( SV_0 ) ) then
            F ( iV, ia_F_I ( 2 ) ) = 10.0_KDR ** L_T
            cycle Value_Host
          endif 
          
          ! if we are closer than 10^-2  to the 
          ! root (eps-eps0)=0, we are switching to 
          ! the secant method, since the table is rather coarse and the
          ! derivatives may be garbage.

          if ( abs ( SV_R - SV_0 )  <  1e-3_KDR * abs ( SV_0 ) ) then
            D2 = ( SV_R - SV_1 ) / ( L_T - L_T_1 )
          end if
          
          !-- FIXME: Need error handling
          !if ( iI == nIterations ) &
          !  call Show ( 'Error, Max iteration reached' )
      
        end do 

      end do Value_Host
      !$OMP  end parallel do
    
    end if
      
  end procedure FindTemperatureKernel
  
  
  module procedure InterpolateTableKernel 
    
#ifdef ENABLE_OMP_OFFLOAD
    !$OMP declare target
#endif
    
    integer ( KDI ) :: &
      iValue, &
      nx,ny,nz,nvars, &
      ix,iy,iz
    real ( KDR ) :: &
      delx, dely, delz, &
      dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi, &
      a1, a2, a3, a4, a5, a6, a7, a8  
    real ( KDR ), dimension ( 8 ) :: &
      fh

    nx = size ( XT )
    ny = size ( YT )
    nz = size ( ZT )
    nvars = size ( T, dim = 4 )  

    !--  determine spacing parameters of (equidistant!!!) table
    dx    = (xt(nx) - xt(1)) / real(nx-1, kind=KDR)
    dy    = (yt(ny) - yt(1)) / real(ny-1, kind=KDR)
    dz    = (zt(nz) - zt(1)) / real(nz-1, kind=KDR)

    dxi   = 1. / dx
    dyi   = 1. / dy
    dzi   = 1. / dz

    dxyi  = dxi * dyi
    dxzi  = dxi * dzi
    dyzi  = dyi * dzi

    dxyzi = dxi * dyi * dzi
      
    !-- determine location in (equidistant!!!) table 
    ix = 2 + INT( (X - xt(1) - 1.e-10_KDR) * dxi ) 
    iy = 2 + INT( (Y - yt(1) - 1.e-10_KDR) * dyi ) 
    iz = 2 + INT( (Z - zt(1) - 1.e-10_KDR) * dzi )
                                               
    ix = MAX( 2, MIN( ix, nx ) )               
    iy = MAX( 2, MIN( iy, ny ) )               
    iz = MAX( 2, MIN( iz, nz ) )               

    !-- set-up auxiliary arrays for Lagrange interpolation
                                                          
    delx = xt(ix) - X                                   
    dely = yt(iy) - Y
    delz = zt(iz) - Z

    fh(1) = T(ix  , iy  , iz,   i_ST)
    fh(2) = T(ix-1, iy  , iz,   i_ST)
    fh(3) = T(ix  , iy-1, iz,   i_ST)
    fh(4) = T(ix  , iy  , iz-1, i_ST)
    fh(5) = T(ix-1, iy-1, iz,   i_ST)
    fh(6) = T(ix-1, iy  , iz-1, i_ST)
    fh(7) = T(ix  , iy-1, iz-1, i_ST)
    fh(8) = T(ix-1, iy-1, iz-1, i_ST)

    !-- set up coefficients of the interpolation polynomial and 
    !   evaluate function values 
    a1 = fh(1)
    a2 = dxi   * ( fh(2) - fh(1) )
    a3 = dyi   * ( fh(3) - fh(1) )
    a4 = dzi   * ( fh(4) - fh(1) )
    a5 = dxyi  * ( fh(5) - fh(2) - fh(3) + fh(1) )
    a6 = dxzi  * ( fh(6) - fh(2) - fh(4) + fh(1) )
    a7 = dyzi  * ( fh(7) - fh(3) - fh(4) + fh(1) )
    a8 = dxyzi * ( fh(8) - fh(1) + fh(2) + fh(3) + &
         fh(4) - fh(5) - fh(6) - fh(7) )

    D2 = -a3
    SV_R  &
      = a1 +  a2 * delx               &
         +  a3 * dely                      &
         +  a4 * delz                      &
         +  a5 * delx * dely               &
         +  a6 * delx * delz               &
         +  a7 * dely * delz               &
         +  a8 * delx * dely * delz

  end procedure InterpolateTableKernel
  
  
end submodule EOS_P_HN_OConnorOtt__Kernel
