#include "Preprocessor"

submodule ( Universe_R_B__Form ) Universe_R_B__Kernel

  use Basics
  
  implicit none

contains


  module procedure Compute_dT_R_E_CGS_Kernel

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
      !$OMP reduction ( min : dT )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then
          if ( abs ( E ( iV )  -  E_Eq ( iV ) ) &
                   /  max ( SqrtTiny, E_Eq ( iV ) )  <  0.1_KDR & 
               .and. E ( iV )  >  SqrtTiny ) & 
          then
            dT  =  min ( dT,  &
                         E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          end if
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP reduction ( min : dT )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then
          if ( abs ( E ( iV )  -  E_Eq ( iV ) ) &
                   /  max ( SqrtTiny, E_Eq ( iV ) )  <  0.1_KDR & 
               .and. E ( iV )  >  SqrtTiny ) & 
          then
            dT  =  min ( dT,  &
                         E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          end if
        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_R_E_CGS_Kernel


  module procedure Compute_dT_ET_CGS_Kernel

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
      !$OMP reduction ( min : dT )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) &
          dT  =  min ( dT,  &
                       E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP reduction ( min : dT )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) &
          dT  =  min ( dT,  &
                       E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_ET_CGS_Kernel


end submodule Universe_R_B__Kernel
