#include "Preprocessor"

submodule ( Universe_R_CC__Form ) Universe_R_CC__Kernel

  use Basics
  
  implicit none

contains


  module procedure Compute_dT_RT_CGS_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( Q )
    
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP reduction ( min : dT_E, dT_N )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then
          dT_E  =  min ( dT_E,  &
                         E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          dT_N  =  min ( dT_N,  &
                         N ( iV )  /  max ( SqrtTiny, abs ( R ( iV ) ) ) )
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP reduction ( min : dT_E, dT_N )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then
          dT_E  =  min ( dT_E,  &
                         E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          dT_N  =  min ( dT_N,  &
                         N ( iV )  /  max ( SqrtTiny, abs ( R ( iV ) ) ) )
        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_RT_CGS_Kernel


end submodule Universe_R_CC__Kernel
