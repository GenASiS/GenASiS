#include "Preprocessor"

submodule ( Slope_DFV_C_F__Form ) Slope_DFV_C_F__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure Compute_C_Kernel

    !-- Compute_Cylindrical_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( S_M_1 )  >  1 )
      lV  =  oV + 1
    end where
    
    uV  =  1
    where ( shape ( S_M_1 )  >  1 )
      uV  =  shape ( S_M_1 )  -  oV
    end where
      
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

            S_M_1 ( iV, jV, kV )  &
              =  S_UD_33 ( iV, jV, kV )  &
                 *  ( A_I_1 ( iV + 1, jV, kV )  -  A_I_1 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else !-- use host
              
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

            S_M_1 ( iV, jV, kV )  &
              =  S_UD_33 ( iV, jV, kV )  &
                 *  ( A_I_1 ( iV + 1, jV, kV )  -  A_I_1 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure Compute_C_Kernel


  module procedure Compute_S_Kernel

    !-- Compute_Spherical_Kernel

    integer ( KDI ) :: &
      iV, jV, kV, &
      dJ
    integer ( KDI ), dimension ( 3 ) :: &
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption

    lV  =  1
    where ( shape ( S_M_1 )  >  1 )
      lV  =  oV + 1
    end where
    
    uV  =  1
    where ( shape ( S_M_1 )  >  1 )
      uV  =  shape ( S_M_1 )  -  oV
    end where

    if ( count ( lV > 1 )  ==  1 ) then
      dJ  =  0  !-- 1D
    else
      dJ  =  1  !-- 2D or 3D
    end if
  
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) firstprivate ( dJ )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

            S_M_1 ( iV, jV, kV )  &
              =  ( S_UD_22 ( iV, jV, kV )  +  S_UD_33 ( iV, jV, kV ) )  &
                 *  0.5_KDR  &
                 *  ( A_I_1 ( iV + 1, jV, kV )  -  A_I_1 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
            S_M_2 ( iV, jV, kV )  &
              =  S_UD_33 ( iV, jV, kV )  &
                 *  ( A_I_2 ( iV, jV + dJ, kV )  -  A_I_2 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else !-- use host
              
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) firstprivate ( dJ )
      do kV  =  lV ( 3 ),  uV ( 3 ) 
        do jV  =  lV ( 2 ),  uV ( 2 )
          do iV  =  lV ( 1 ),  uV ( 1 )

            S_M_1 ( iV, jV, kV )  &
              =  ( S_UD_22 ( iV, jV, kV )  +  S_UD_33 ( iV, jV, kV ) )  &
                 *  0.5_KDR  &
                 *  ( A_I_1 ( iV + 1, jV, kV )  -  A_I_1 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
            S_M_2 ( iV, jV, kV )  &
              =  S_UD_33 ( iV, jV, kV )  &
                 *  ( A_I_2 ( iV, jV + dJ, kV )  -  A_I_2 ( iV, jV, kV ) )  &
                 /  V ( iV, jV, kV )
 
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if !-- UseDevice
        
  end procedure Compute_S_Kernel


end submodule Slope_DFV_C_F__Kernel
