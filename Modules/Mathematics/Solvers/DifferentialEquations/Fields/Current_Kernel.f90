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

    nValues = size ( DFV_I )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV ) 
      do iV = 1, nValues
        DFV_I ( iV )  =  1.0_KDR
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
      
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV ) 
      do iV = 1, nValues
        DFV_I ( iV )  =  1.0_KDR
      end do !-- iV
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
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nValues
        AP_I ( iV ) = max ( 0.0_KDR, + LP_IL ( iV ), + LP_IR ( iV ) )
        AM_I ( iV ) = max ( 0.0_KDR, - LM_IL ( iV ), - LM_IR ( iV ) )
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST ) private ( iV ) 
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
    real ( KDR ) :: &
      RealTiny
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
     
    nV = size ( F_I )
    RealTiny = tiny ( 0.0_KDR )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
      !$OMP& firstprivate ( RealTiny )
      do iV = 1, nV
        F_I ( iV ) &
          =  (    AP_I ( iV ) * F_IL ( iV ) &
               +  AM_I ( iV ) * F_IR ( iV ) &
               -  DF_I ( iV ) * AP_I ( iV ) * AM_I ( iV ) &
                  * ( U_IR ( iV ) - U_IL ( iV ) ) ) &
             /  max ( AP_I ( iV ) + AM_I ( iV ), RealTiny )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP& firstprivate ( RealTiny )
      do iV = 1, nV
        F_I ( iV ) &
          =  (    AP_I ( iV ) * F_IL ( iV ) &
               +  AM_I ( iV ) * F_IR ( iV ) &
               -  DF_I ( iV ) * AP_I ( iV ) * AM_I ( iV ) &
                  * ( U_IR ( iV ) - U_IL ( iV ) ) ) &
             /  max ( AP_I ( iV ) + AM_I ( iV ), RealTiny )
      end do
      !$OMP  end parallel do
    
    end if

  end procedure ComputeFluxes_HLL_Kernel


end submodule Current_Kernel
