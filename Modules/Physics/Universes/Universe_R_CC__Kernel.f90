#include "Preprocessor"

submodule ( Universe_R_CC__Form ) Universe_R_CC__Kernel

  use Basics
  
  implicit none

contains


  module procedure Compute_dT_RI_CGS_Kernel

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
          if ( abs ( E ( iV )  -  E_Eq ( iV ) ) &
                   /  max ( SqrtTiny, E_Eq ( iV ) )  <  0.1_KDR & 
               .and. E ( iV )  >  SqrtTiny ) & 
          then
            dT_E  =  min ( dT_E,  &
                           E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          end if
          if ( abs ( N ( iV )  -  N_Eq ( iV ) ) &
                   /  max ( SqrtTiny, N_Eq ( iV ) )  <  0.1_KDR &
               .and. N ( iV )  >  SqrtTiny ) & 
          then
            dT_N  =  min ( dT_N,  &
                           N ( iV )  /  max ( SqrtTiny, abs ( R ( iV ) ) ) )
          end if
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
          if ( abs ( E ( iV )  -  E_Eq ( iV ) ) &
                   /  max ( SqrtTiny, E_Eq ( iV ) )  <  0.1_KDR & 
               .and. E ( iV )  >  SqrtTiny ) & 
          then
            dT_E  =  min ( dT_E,  &
                           E ( iV )  /  max ( SqrtTiny, abs ( Q ( iV ) ) ) )
          end if
          if ( abs ( N ( iV )  -  N_Eq ( iV ) ) &
                   /  max ( SqrtTiny, N_Eq ( iV ) )  <  0.1_KDR &
               .and. N ( iV )  >  SqrtTiny ) & 
          then
            dT_N  =  min ( dT_N,  &
                           N ( iV )  /  max ( SqrtTiny, abs ( R ( iV ) ) ) )
          end if
        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_RI_CGS_Kernel


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


  module procedure Compute_dT_RK_R_CGS_Kernel

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
      
    nV  =  size ( E )
    
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP reduction ( min : dT_E, dT_N, dT_S )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then
          dT_E  =  min ( dT_E,  & 
                         sqrt ( 1.0e-4  *  abs ( E ( iV ) ) &
                                /  max ( SqrtTiny, abs ( E_E ( iV ) ) ) ) &
                         *  dT )
          dT_N  =  min ( dT_N,  &
                         sqrt ( 1.0e-4  *  abs ( N ( iV ) ) &
                                /  max ( SqrtTiny, abs ( E_N ( iV ) ) ) ) &
                         *  dT )
          dT_S  =  min ( dT_S,  &
                         sqrt ( 1.0e-4  *  ( 1.0e-2 * abs ( E ( iV ) ) &
                                             +  abs ( S ( iV ) ) ) &
                                /  max ( SqrtTiny, abs ( E_S ( iV ) ) ) ) &
                         *  dT )
        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_RK_R_CGS_Kernel


end submodule Universe_R_CC__Kernel
