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


  module procedure Compute_dT_RK_F_CGS_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      Tolerance, &
      FE_E, FE_S_1, FE_S_2, FE_S_3, &  !-- Fractional errors
       F_E,  F_S_1,  F_S_2,  F_S_3     !-- Factors
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( E )
    
    SqrtTiny   =  sqrt ( tiny ( 0.0_KDR ) )
    Tolerance  =  1.0e-4_KDR

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( FE_E, FE_S_1, FE_S_2, FE_S_3 ) &
      !$OMP private ( F_E, F_S_1, F_S_2, F_S_3 ) &
      !$OMP reduction ( min : dT_E, dT_S_1, dT_S_2, dT_S_3 )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then

          FE_E  =  max ( Tolerance * SqrtTiny, abs ( E_E ( iV ) ) )  &
                 /  max ( SqrtTiny, abs ( E ( iV ) ) )
           F_E  =  sqrt ( Tolerance / FE_E )
          if ( F_E  <  1.0_KDR ) &
            F_E  =  max ( F_E, 0.2_KDR ) 
          dT_E  =  min ( dT_E, F_E * dT )

          FE_S_1  =  max ( Tolerance * SqrtTiny, abs ( E_S_1 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_1 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_1 ( iV ) ) )
           F_S_1  =  sqrt ( Tolerance / FE_S_1 )
          if ( F_S_1  <  1.0_KDR ) &
            F_S_1  =  max ( F_S_1, 0.2_KDR ) 
          dT_S_1  =  min ( dT_S_1, F_S_1 * dT )

          FE_S_2  =  max ( Tolerance * SqrtTiny, abs ( E_S_2 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_2 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_2 ( iV ) ) )
           F_S_2  =  sqrt ( Tolerance / FE_S_2 )
          if ( F_S_2  <  1.0_KDR ) &
            F_S_2  =  max ( F_S_2, 0.2_KDR ) 
          dT_S_2  =  min ( dT_S_2, F_S_2 * dT )

          FE_S_3  =  max ( Tolerance * SqrtTiny, abs ( E_S_3 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_3 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_3 ( iV ) ) )
           F_S_3  =  sqrt ( Tolerance / FE_S_3 )
          if ( F_S_3  <  1.0_KDR ) &
            F_S_3  =  max ( F_S_3, 0.2_KDR ) 
          dT_S_3  =  min ( dT_S_3, F_S_3 * dT )

        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( FE_E, FE_S_1, FE_S_2, FE_S_3 ) &
      !$OMP private ( F_E, F_S_1, F_S_2, F_S_3 ) &
      !$OMP reduction ( min : dT_E, dT_S_1, dT_S_2, dT_S_3 )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then

          FE_E  =  max ( Tolerance * SqrtTiny, abs ( E_E ( iV ) ) )  &
                 /  max ( SqrtTiny, abs ( E ( iV ) ) )
           F_E  =  sqrt ( Tolerance / FE_E )
          if ( F_E  <  1.0_KDR ) &
            F_E  =  max ( F_E, 0.2_KDR ) 
          dT_E  =  min ( dT_E, F_E * dT )

          FE_S_1  =  max ( Tolerance * SqrtTiny, abs ( E_S_1 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_1 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_1 ( iV ) ) )
           F_S_1  =  sqrt ( Tolerance / FE_S_1 )
          if ( F_S_1  <  1.0_KDR ) &
            F_S_1  =  max ( F_S_1, 0.2_KDR ) 
          dT_S_1  =  min ( dT_S_1, F_S_1 * dT )

          FE_S_2  =  max ( Tolerance * SqrtTiny, abs ( E_S_2 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_2 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_2 ( iV ) ) )
           F_S_2  =  sqrt ( Tolerance / FE_S_2 )
          if ( F_S_2  <  1.0_KDR ) &
            F_S_2  =  max ( F_S_2, 0.2_KDR ) 
          dT_S_2  =  min ( dT_S_2, F_S_2 * dT )

          FE_S_3  =  max ( Tolerance * SqrtTiny, abs ( E_S_3 ( iV ) ) )  &
                   /  max ( SqrtTiny, &
                            1.0e-4 * abs ( E ( iV ) )  +  abs ( S_3 ( iV ) ) )
!                   /  max ( SqrtTiny, abs ( S_3 ( iV ) ) )
           F_S_3  =  sqrt ( Tolerance / FE_S_3 )
          if ( F_S_3  <  1.0_KDR ) &
            F_S_3  =  max ( F_S_3, 0.2_KDR ) 
          dT_S_3  =  min ( dT_S_3, F_S_3 * dT )

        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_RK_F_CGS_Kernel


  module procedure Compute_dT_RK_R_CGS_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      Tolerance, &
      FE_E, FE_S_1, FE_S_2, FE_S_3, &  !-- Fractional errors
       F_E,  F_S_1,  F_S_2,  F_S_3     !-- Factors
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( E )
    
    SqrtTiny   =  sqrt ( tiny ( 0.0_KDR ) )
    Tolerance  =  1.0e-4_KDR

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( FE_E, FE_S_1, FE_S_2, FE_S_3 ) &
      !$OMP private ( F_E, F_S_1, F_S_2, F_S_3 ) &
      !$OMP reduction ( min : dT_E, dT_S_1, dT_S_2, dT_S_3 )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then

            FE_E  =  max ( Tolerance * SqrtTiny, abs ( E_E ( iV ) ) )  &
                   /  max ( SqrtTiny, abs ( E ( iV ) ) )
             F_E  =  sqrt ( Tolerance / FE_E )
            if ( F_E  <  1.0_KDR ) &
              F_E  =  max ( F_E, 0.2_KDR ) 
            dT_E  =  min ( dT_E, F_E * dT )

            FE_S_1  =  max ( Tolerance * SqrtTiny, abs ( E_S_1 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_1 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_1 ( iV ) ) )
             F_S_1  =  sqrt ( Tolerance / FE_S_1 )
            if ( F_S_1  <  1.0_KDR ) &
              F_S_1  =  max ( F_S_1, 0.2_KDR ) 
            dT_S_1  =  min ( dT_S_1, F_S_1 * dT )

            FE_S_2  =  max ( Tolerance * SqrtTiny, abs ( E_S_2 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_2 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_2 ( iV ) ) )
             F_S_2  =  sqrt ( Tolerance / FE_S_2 )
            if ( F_S_2  <  1.0_KDR ) &
              F_S_2  =  max ( F_S_2, 0.2_KDR ) 
            dT_S_2  =  min ( dT_S_2, F_S_2 * dT )

            FE_S_3  =  max ( Tolerance * SqrtTiny, abs ( E_S_3 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_3 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_3 ( iV ) ) )
             F_S_3  =  sqrt ( Tolerance / FE_S_3 )
            if ( F_S_3  <  1.0_KDR ) &
              F_S_3  =  max ( F_S_3, 0.2_KDR ) 
            dT_S_3  =  min ( dT_S_3, F_S_3 * dT )

        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( FE_E, FE_S_1, FE_S_2, FE_S_3 ) &
      !$OMP private ( F_E, F_S_1, F_S_2, F_S_3 ) &
      !$OMP reduction ( min : dT_E, dT_S_1, dT_S_2, dT_S_3 )
      do iV = 1, nV
        if ( ProperCell ( iV ) ) then

            FE_E  =  max ( Tolerance * SqrtTiny, abs ( E_E ( iV ) ) )  &
                   /  max ( SqrtTiny, abs ( E ( iV ) ) )
             F_E  =  sqrt ( Tolerance / FE_E )
            if ( F_E  <  1.0_KDR ) &
              F_E  =  max ( F_E, 0.2_KDR ) 
            dT_E  =  min ( dT_E, F_E * dT )

            FE_S_1  =  max ( Tolerance * SqrtTiny, abs ( E_S_1 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_1 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_1 ( iV ) ) )
             F_S_1  =  sqrt ( Tolerance / FE_S_1 )
            if ( F_S_1  <  1.0_KDR ) &
              F_S_1  =  max ( F_S_1, 0.2_KDR ) 
            dT_S_1  =  min ( dT_S_1, F_S_1 * dT )

            FE_S_2  =  max ( Tolerance * SqrtTiny, abs ( E_S_2 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_2 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_2 ( iV ) ) )
             F_S_2  =  sqrt ( Tolerance / FE_S_2 )
            if ( F_S_2  <  1.0_KDR ) &
              F_S_2  =  max ( F_S_2, 0.2_KDR ) 
            dT_S_2  =  min ( dT_S_2, F_S_2 * dT )

            FE_S_3  =  max ( Tolerance * SqrtTiny, abs ( E_S_3 ( iV ) ) )  &
                     /  max ( SqrtTiny, &
                              1.0e-4 * abs ( E ( iV ) )  +  abs ( S_3 ( iV ) ) )
  !                   /  max ( SqrtTiny, abs ( S_3 ( iV ) ) )
             F_S_3  =  sqrt ( Tolerance / FE_S_3 )
            if ( F_S_3  <  1.0_KDR ) &
              F_S_3  =  max ( F_S_3, 0.2_KDR ) 
            dT_S_3  =  min ( dT_S_3, F_S_3 * dT )

        end if
      end do
      !$OMP  end parallel do
    end if
      
  end procedure Compute_dT_RK_R_CGS_Kernel


end submodule Universe_R_B__Kernel
