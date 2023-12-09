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
      SqrtTiny, &
      Tolerance, &
      FE_E, FE_N, FE_S, &  !-- Fractional errors
       F_E,  F_N,  F_S     !-- Factors
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( E )
    
    SqrtTiny   =  sqrt ( tiny ( 0.0_KDR ) )
    Tolerance  =  1.0e-4_KDR

    if ( UseDevice ) then
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( FE_E, FE_N, FE_S, F_E, F_N, F_S ) &
      !$OMP reduction ( min : dT_E, dT_N, dT_S )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then

          FE_E  =  max ( Tolerance * SqrtTiny, abs ( E_E ( iV ) ) )  &
                   /  max ( SqrtTiny, abs ( E ( iV ) ) )
          FE_N  =  max ( Tolerance * SqrtTiny, abs ( E_N ( iV ) ) )  &
                   /  max ( SqrtTiny, abs ( N ( iV ) ) )
          FE_S  =  max ( Tolerance * SqrtTiny, abs ( E_S ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S ( iV ) ) )

          F_E  =  sqrt ( Tolerance / FE_E )
          F_N  =  sqrt ( Tolerance / FE_N )
          F_S  =  sqrt ( Tolerance / FE_S )

! call Show ( iV, '>>> iV' )
! call Show ( F_E, '>>>>>> F_E' )
! call Show ( F_N, '>>>>>> F_N' )
! call Show ( F_S, '>>>>>> F_S' )
! call Show ( E_S ( iV ), '>>>>>> E_S ( iV )' )
! call Show ( S ( iV ), '>>>>>> S ( iV )' )

          ! if ( F_E  >  1.0_KDR ) &
          !   F_E  =  min ( F_E, 10.0_KDR ) 
          ! if ( F_N  >  1.0_KDR ) &
          !   F_N  =  min ( F_N, 10.0_KDR ) 
          ! if ( F_S  >  1.0_KDR ) &
          !   F_S  =  min ( F_S, 10.0_KDR ) 
          
          ! if ( F_E  <  1.0_KDR ) &
          !   F_E  =  max ( F_E, 0.2_KDR ) 
          ! if ( F_N  <  1.0_KDR ) &
          !   F_N  =  max ( F_N, 0.2_KDR ) 
          ! if ( F_S  <  1.0_KDR ) &
          !   F_S  =  max ( F_S, 0.2_KDR ) 
          
          dT_E  =  min ( dT_E, F_E * dT )
          dT_N  =  min ( dT_N, F_N * dT )
          dT_S  =  min ( dT_S, F_S * dT )

        end if
      end do
      !$OMP  end parallel do
! call Show ( maxval ( max ( Tolerance * SqrtTiny, abs ( E_S ( 3: ) ) )  &
!                    /  max ( SqrtTiny, abs ( S ( 3: ) ) ) ), '>>> maxval FE_S' )
! call Show ( maxloc ( max ( Tolerance * SqrtTiny, abs ( E_S ) )  &
!                    /  max ( SqrtTiny, abs ( S ( 3: ) ) ) ), '>>> FE_S' )
    end if
      
  end procedure Compute_dT_RK_R_CGS_Kernel


end submodule Universe_R_CC__Kernel
