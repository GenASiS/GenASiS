#include "Preprocessor"

submodule ( Fluid_P_P__Form ) Fluid_P_P__Kernel

  use Basics
  
  implicit none
  
contains

  module procedure Apply_EOS_P_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( P )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        Gamma ( iV )  =  Gamma_0
        P     ( iV )  =  E ( iV )  *  ( Gamma_0 - 1.0_KDR ) 
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV ) > 0.0_KDR ) then
          K ( iV )  =  P ( iV ) / ( N ( iV ) ** Gamma_0 )
        else
          K ( iV )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( K ( iV ) > 0.0_KDR ) then
          SB ( iV )  =    log ( K ( iV ) / K0 ) / ( Gamma_0 - 1.0_KDR )
        else
          SB ( iV )  =  - 0.1 * huge ( 1.0_KDR )
        end if
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        Gamma ( iV )  =  Gamma_0
        P     ( iV )  =  E ( iV )  *  ( Gamma_0 - 1.0_KDR ) 
      end do !-- iV
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV ) > 0.0_KDR ) then
          K ( iV )  =  P ( iV ) / ( N ( iV ) ** Gamma_0 )
        else
          K ( iV )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( K ( iV ) > 0.0_KDR ) then
          SB ( iV )  =    log ( K ( iV ) / K0 ) / ( Gamma_0 - 1.0_KDR )
        else
          SB ( iV )  =  - 0.1 * huge ( 1.0_KDR )
        end if
      end do !-- iV
      !$OMP end parallel do
      
    end if

  end procedure Apply_EOS_P_Kernel


end submodule Fluid_P_P__Kernel
