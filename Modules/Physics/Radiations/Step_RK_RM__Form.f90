module Step_RK_RM__Form

  !-- Step_RungeKutta_RadiationMoments_Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use RadiationMoments_BM__Form
  use Interactions_BM__Form

  implicit none
  private

  type, public, extends ( Step_RK_CS_CS_Form ) :: Step_RK_RM_Form
    real ( KDR ), dimension ( : ), allocatable :: &
      Residual_R_E, &
      Residual_F_E
  contains
    procedure, private, pass :: &
      Initialize_CS_CS
    final :: &
      Finalize
    procedure, public, pass :: &
      SolveUpdateImplicit
  end type Step_RK_RM_Form


contains


  subroutine Initialize_CS_CS &
               ( S, CS_1, CS_2, NameOption, ImplicitExplicitOption, &
                 nStagesOption )

    class ( Step_RK_RM_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS_1, &
      CS_2
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_RM'

    call S % Step_RK_CS_CS_Form % Initialize &
           ( CS_1, CS_2, NameOption, ImplicitExplicitOption, nStagesOption )

    allocate ( S % Residual_R_E ( S % MaxImplicitIterations ) )
    allocate ( S % Residual_F_E ( S % MaxImplicitIterations ) )
    S % Residual_R_E  =  0.0_KDR
    S % Residual_F_E  =  0.0_KDR

  end subroutine Initialize_CS_CS


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_RM_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Residual_F_E ) ) &
      deallocate ( S % Residual_F_E )
    if ( allocated ( S % Residual_R_E ) ) &
      deallocate ( S % Residual_R_E )

  end subroutine Finalize


  subroutine SolveUpdateImplicit  ( S, T, dT, iS )

    class ( Step_RK_RM_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC, &
      iV, &
      iI, &
      nV
    integer ( KDI ) :: &
      iEnergy_R, &
      iEnergy_F
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, &
      iMomentum_F
    real ( KDR ) :: &
       E_R_P, &
       E_F_P, &
      dE_R, &
      dE_F, &
      SqrtTiny

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    associate &
      ( S_R  =>  S % Step_CS_1, &
        S_F  =>  S % Step_CS_2 )
    select type ( R  =>  S_R % CurrentSet )
      class is ( RadiationMoments_BM_Form )
    select type ( F  =>  S_F % CurrentSet )
      class is ( Fluid_P_Form )
    select type ( I  =>  R % Interactions )
      class is ( Interactions_BM_Form )
    associate &
      ( Y_I_R     =>  S_R % Intermediate, &
        Y_I_F     =>  S_F % Intermediate, &
         KK_R     =>  S_R % SlopeStageImplicit ( iS ) % Element, &
         KK_F     =>  S_F % SlopeStageImplicit ( iS ) % Element, &
         AA       =>  S_R % AA ( iS ) % Value ( iS ), &
         Tol      =>  S   % ImplicitTolerance, &
         Max_I    =>  S   % MaxImplicitIterations, &
         Res_R_E  =>  S   % Residual_R_E, &
         Res_F_E  =>  S   % Residual_F_E, &
         ID       =>  S   % ImplicitDiagnostics ( iS ) )

    call Search &
           ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )

    call Search &
           ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )
      associate &
        (     I_V  =>      I % Storage ( iC ) % Value, &
              R_V  =>      R % Storage ( iC ) % Value, &
          Y_I_R_V  =>  Y_I_R % Storage ( iC ) % Value, &
          Y_I_F_V  =>  Y_I_F % Storage ( iC ) % Value, &
           KK_R_V  =>   KK_R % Storage ( iC ) % Value, &
           KK_F_V  =>   KK_F % Storage ( iC ) % Value, &
             ID_V  =>     ID % Storage ( iC ) % Value )
      associate &
        (   Xi_J      =>  I_V ( :, I % EMISSIVITY_J ), &
            Xi_H      =>  I_V ( :, I % EMISSIVITY_H ), &
           Chi_J      =>  I_V ( :, I % OPACITY_J ), &
           Chi_H      =>  I_V ( :, I % OPACITY_H ), &
               J      =>  R_V ( :, R % ENERGY_DENSITY_C ), &
               H_1    =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_1 ), &
               H_2    =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_2 ), &
               H_3    =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_3 ), &
             E_R      =>  R_V ( :, R % ENERGY_DENSITY_B ), &
             S_R_1    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
             S_R_2    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
             S_R_3    =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
             E_F      =>  R_V ( :, F % ENERGY_DENSITY_B ), &
             S_F_1    =>  R_V ( :, F % MOMENTUM_DENSITY_D_1 ), &
             S_F_2    =>  R_V ( :, F % MOMENTUM_DENSITY_D_2 ), &
             S_F_3    =>  R_V ( :, F % MOMENTUM_DENSITY_D_3 ), &
             E_R_0    =>  Y_I_R_V ( :, iEnergy_R ), &
             S_R_1_0  =>  Y_I_R_V ( :, iMomentum_R ( 1 ) ), &
             S_R_2_0  =>  Y_I_R_V ( :, iMomentum_R ( 2 ) ), &
             S_R_3_0  =>  Y_I_R_V ( :, iMomentum_R ( 3 ) ), &
             E_F_0    =>  Y_I_F_V ( :, iEnergy_F ), &
             S_F_1_0  =>  Y_I_F_V ( :, iMomentum_F ( 1 ) ), &
             S_F_2_0  =>  Y_I_F_V ( :, iMomentum_F ( 2 ) ), &
             S_F_3_0  =>  Y_I_F_V ( :, iMomentum_F ( 3 ) ), &
          KK_R_E      =>  KK_R_V ( :, iEnergy_R ), &
          KK_R_S_1    =>  KK_R_V ( :, iMomentum_R ( 1 ) ), &
          KK_R_S_2    =>  KK_R_V ( :, iMomentum_R ( 2 ) ), &
          KK_R_S_3    =>  KK_R_V ( :, iMomentum_R ( 3 ) ), &
          KK_F_E      =>  KK_F_V ( :, iEnergy_F ), &
          KK_F_S_1    =>  KK_F_V ( :, iMomentum_F ( 1 ) ), &
          KK_F_S_2    =>  KK_F_V ( :, iMomentum_F ( 2 ) ), &
          KK_F_S_3    =>  KK_F_V ( :, iMomentum_F ( 3 ) ), &
             N_I      =>  ID_V ( :, ID % N_ITERATIONS ), &
             R_Max    =>  ID_V ( :, ID % RESIDUAL_MAX ), &
          ProperCell  =>  C % ProperCell )

      select type ( G  =>  R % Geometry )
      class is ( Gravitation_G_Form )

        associate &
          ( GSV  =>  G % Storage ( iC ) % Value )
        associate &
          ( M_DD_11  =>  GSV ( :, G % METRIC_F_DD_11 ), &
            M_DD_22  =>  GSV ( :, G % METRIC_F_DD_22 ), &
            M_DD_33  =>  GSV ( :, G % METRIC_F_DD_33 ) )

        nV  =  size ( ProperCell )

!call Show ( '>>> Stage' )
!call Show ( iS, '>>> iS' )

        do iV = 1, nV
          if ( ProperCell ( iV ) ) then      

!if ( iV == 3 ) &
!  call Show ( iV, '>>> iV' )

            !-- Iterate radiation and fluid energy to convergence

            iI  =  0
            Implicit: do 

              iI  =  iI + 1

!if ( iV == 3 ) then
!  call Show ( '>>> Iteration' )
!  call Show ( iI, '>>> iI' )
!end if

              !-- Compute interactions

              call I % Compute ( iC, iV )

              !-- Compute energy updates

              KK_R_E ( iV )  &
                =  ( Xi_J ( iV )  -  Chi_J ( iV )  *  J ( iV ) ) &
                   /  ( 1.0_KDR  +  Chi_J ( iV ) * dT )

              KK_F_E ( iV )  =  - KK_R_E ( iV )

              !-- Apply energy updates

              E_R ( iV )  =  E_R_0 ( iV )  +  dT * AA * KK_R_E ( iV )
              E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )

              !-- Exit test

              if ( iI  >  1 ) then

                dE_R  =  abs ( E_R ( iV )  -  E_R_P )  &
                         /  max ( abs ( E_R_0 ( iV ) ), SqrtTiny )
                dE_F  =  abs ( E_F ( iV )  -  E_F_P )  &
                         /  max ( abs ( E_F_0 ( iV ) ), SqrtTiny )

                Res_R_E ( iI )  =  dE_R
                Res_F_E ( iI )  =  dE_F

!if ( iV == 3 ) then
!  call Show ( dE_R, '>>> dE_R' )
!  call Show ( dE_F, '>>> dE_F' )
!end if

                if ( dE_R  <  Tol .and. dE_F  <  Tol ) then
!                if ( dE_F  <  Tol ) then
!if ( iV == 3 ) then
!  call Show ( '>>> Solution reached' )
!  call Show ( iI, '>>> iI' )
!end if
                  N_I   ( iV )  =  iI
                  R_Max ( iV )  =  max ( dE_R, dE_F )

                  exit Implicit
                end if

                if ( iI  ==  Max_I ) then
                  call Show ( 'Max iterations reached', CONSOLE % ERROR )
                  call Show ( iV, 'iV', CONSOLE % ERROR )
                  call Show ( Res_R_E ( : iI ), 'Residual_R_E' )
                  call Show ( Res_F_E ( : iI ), 'Residual_F_E' )
                  call Show ( 'Step_RK_RM__Form', 'subroutine', &
                              CONSOLE % ERROR )
                  call Show ( 'SolveUpdateImplicit', 'subroutine', &
                              CONSOLE % ERROR )
                  call PROGRAM_HEADER % Abort ( )
                end if

              end if

              !-- Prepare for next iteration

              E_R_P  =  E_R ( iV )
              E_F_P  =  E_F ( iV )

              call R % ComputeFromBalanced ( iC, iV )
              call F % ComputeFromBalanced ( iC, iV )

            end do Implicit

            !--- Momentum update

            KK_R_S_1 ( iV )  &
              =  ( Xi_H ( iV )  &
                   -  Chi_H ( iV )  *  M_DD_11 ( iV ) * H_1 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
            KK_R_S_2 ( iV )  &
              =  ( Xi_H ( iV )  &
                   -  Chi_H ( iV )  *  M_DD_22 ( iV ) * H_2 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )
            KK_R_S_3 ( iV )  &
              =  ( Xi_H ( iV )  &
                   -  Chi_H ( iV )  *  M_DD_33 ( iV ) * H_3 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_H ( iV ) * dT )

            KK_F_S_1 ( iV )  =  - KK_R_S_1 ( iV )
            KK_F_S_2 ( iV )  =  - KK_R_S_2 ( iV )
            KK_F_S_3 ( iV )  =  - KK_R_S_3 ( iV )

          else
            KK_R_E   ( iV )  =  0.0_KDR
            KK_R_S_1 ( iV )  =  0.0_KDR
            KK_R_S_2 ( iV )  =  0.0_KDR
            KK_R_S_3 ( iV )  =  0.0_KDR
            KK_F_E   ( iV )  =  0.0_KDR
            KK_F_S_1 ( iV )  =  0.0_KDR
            KK_F_S_2 ( iV )  =  0.0_KDR
            KK_F_S_3 ( iV )  =  0.0_KDR
          end if !-- ProperCell
        end do !-- iV

        end associate !-- M_DD_11, etc.
        end associate !-- GSV

      class default
        call Show ( 'Gravitation type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_RM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- G

      end associate !-- Xi_J, etc.
      end associate !-- IntV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_RM__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- KK_R, etc.
    end select !-- I
    end select !-- F
    end select !-- R
    end associate !-- S_1, etc.

  end subroutine SolveUpdateImplicit


end module Step_RK_RM__Form
