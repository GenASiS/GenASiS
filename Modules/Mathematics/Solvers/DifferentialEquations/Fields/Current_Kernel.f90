#include "Preprocessor"

submodule ( Current_Template ) Current_Kernel

  use Basics
  
  implicit none

contains


  module procedure SetDiffusionFactorUnity

    integer ( KDI ) :: &
      iV, jV, &
      nValues, &
      nVariables
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues    = size ( DFV_I, dim = 1 )
    nVariables = size ( DFV_I, dim = 2 )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV ) 
      do jV = 1, nVariables
        do iV = 1, nValues
          DFV_I ( iV, jV )  =  1.0_KDR
        end do !-- iV
      end do !-- jV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
      
      !$OMP  parallel do collapse ( 2 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV ) 
      do jV = 1, nVariables
        do iV = 1, nValues
          DFV_I ( iV, jV )  =  1.0_KDR
        end do !-- iV
      end do !-- jV
      !$OMP end parallel do
      
    end if

  end procedure SetDiffusionFactorUnity


  module procedure ComputeSolverSpeeds_HLL_Kernel
  
    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( AP_I )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        AP_I ( iV ) = max ( 0.0_KDR, + LP_IL ( iV ), + LP_IR ( iV ) )
        AM_I ( iV ) = max ( 0.0_KDR, - LM_IL ( iV ), - LM_IR ( iV ) )
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        AP_I ( iV ) = max ( 0.0_KDR, + LP_IL ( iV ), + LP_IR ( iV ) )
        AM_I ( iV ) = max ( 0.0_KDR, - LM_IL ( iV ), - LM_IR ( iV ) )
      end do !-- iV
      !$OMP end parallel do
      
    end if

  end procedure ComputeSolverSpeeds_HLL_Kernel


  module procedure ComputeFluxes_HLL_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
     
    nV = size ( F_I )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nV
        F_I ( iV ) &
          =  (    AP_I ( iV ) * F_IL ( iV ) &
               +  AM_I ( iV ) * F_IR ( iV ) &
               -  DF_I ( iV ) * AP_I ( iV ) * AM_I ( iV ) &
                  * ( U_IR ( iV ) - U_IL ( iV ) ) ) &
             /  max ( AP_I ( iV ) + AM_I ( iV ), tiny ( 0.0_KDR ) )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do private ( iV )
      do iV = 1, nV
        F_I ( iV ) &
          =  (    AP_I ( iV ) * F_IL ( iV ) &
               +  AM_I ( iV ) * F_IR ( iV ) &
               -  DF_I ( iV ) * AP_I ( iV ) * AM_I ( iV ) &
                  * ( U_IR ( iV ) - U_IL ( iV ) ) ) &
             /  max ( AP_I ( iV ) + AM_I ( iV ), tiny ( 0.0_KDR ) )
      end do
      !$OMP end parallel do
    
    end if

  end procedure ComputeFluxes_HLL_Kernel


end submodule Current_Kernel
