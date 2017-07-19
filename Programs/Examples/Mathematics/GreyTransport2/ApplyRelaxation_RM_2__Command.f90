module ApplyRelaxation_RM_2__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    ApplyRelaxation_RM_2

    private :: &
      ComputeCoefficients, &
      ComputeComovingIncrements, &
      ComputeConservedIncrements

contains

  
  subroutine ApplyRelaxation_RM_2 &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 RadiationMoments, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      RadiationMoments
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iV, &
      nV, &
      iEnergy, &
      iMomentum_1
    real ( KDR ) :: &
      dJ, dH_1
    real ( KDR ), dimension ( 3 ) :: &
      k_D  !-- Eddington tensor components ( K / J, Stress / EnergyDensity )
    real ( KDR ), dimension ( 2 ) :: &
      B
    real ( KDR ), dimension ( 2, 2 ) :: &
      A
    class ( GeometryFlatForm ), pointer :: &
      G

    call Show ( 'ApplyRelaxation_RM_2', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    nV = size ( IncrementExplicit % Value, dim = 1 )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    G => Chart % Geometry ( )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, &
                  iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1 )

    !$OMP parallel do private ( iV, A, B, dJ, dH_1, k_D )
    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

call Show ( iV, '>>> iV' )
call Show ( RadiationMoments % Name, 'RadiationMoments' )

      call ComputeCoefficients &
             ( RM, &
               I  % Value ( iV, I % EMISSIVITY_J ), &
               I  % Value ( iV, I % OPACITY_J ), &
               I  % Value ( iV, I % OPACITY_H ), &
               RM % Value ( iV, RM % COMOVING_ENERGY ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               RM % Value ( iV, RM % STRESS_FACTOR ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
               G  % Value ( iV, G % METRIC_DD_22 ), &
               G  % Value ( iV, G % METRIC_DD_33 ), &
               TimeStep, A, B, k_D )

call Show ( A, '>>> A' )
call Show ( B, '>>> B' )
call Show ( k_D ( 1 ), '>>> k_D ( 1 )' )

      call ComputeComovingIncrements &
             ( A, B, dJ, dH_1 )
      !-- FIXME: total hack, use Sources_RM to accessibly store comoving 
      !          increments 
      call ComputeConservedIncrements &
             ( IncrementExplicit % Value ( iV, iEnergy ), &
               IncrementExplicit % Value ( iV, iMomentum_1 ), &
               SRM % Value ( iV, SRM % NET_EMISSION_E ), &
               SRM % Value ( iV, SRM % NET_EMISSION_S_D ( 1 ) ), &
               k_D, dJ, dH_1, &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ) )

call Show ( RM % Value ( iV, RM % COMOVING_ENERGY ), '>>> J' )
call Show ( dJ, '>>> dJ' )
call Show ( RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), '>>> H_1' )
call Show ( dH_1, '>>> dH_1' )
call Show ( RM % Value ( iV, RM % CONSERVED_ENERGY ), '>>> E' )
call Show ( IncrementExplicit % Value ( iV, iEnergy ), '>>> dE' )
call Show ( RM % Value ( iV, RM % CONSERVED_MOMENTUM_D ( 1 ) ), '>>> S_1' )
call Show ( IncrementExplicit % Value ( iV, iMomentum_1 ), '>>> dS_1' )

    end do
    !$OMP end parallel do

    end associate !-- I
    end select !-- Chart
    end select !-- SRM
    end select !-- RM
    nullify ( G )
    
  end subroutine ApplyRelaxation_RM_2


  subroutine ComputeCoefficients &
               ( RM, Xi_J, Chi_J, Chi_H, J, H_1, H_2, H_3, SF, V_1, &
                 M_DD_22, M_DD_33, dt, A, B, k_D )

    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_1, H_2, H_3, SF, V_1, &
      M_DD_22, M_DD_33, &
      dt
    real ( KDR ), dimension ( 2, 2 ), intent ( out ) :: &
      A
    real ( KDR ), dimension ( 2 ), intent ( out ) :: &
      B
    real ( KDR ), dimension ( 3 ), intent ( out ) :: &
      k_D

    call RM % ComputeComovingStress_D &
           ( k_D, [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 1 )

    if ( J > 0.0_KDR ) then
      k_D  =  k_D / J
    else
      k_D  =  0.0_KDR
    end if

    A ( 1, 1 )  =  1.0_KDR  +  Chi_J * dt

    A ( 2, 1 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_1

    A ( 1, 2 )  =  V_1 * ( 1.0_KDR  +  Chi_J * dt )  +  k_D ( 1 ) * V_1

    A ( 2, 2 )  =  1.0_KDR  +  Chi_H * dt

    B ( 1 )  =  ( Xi_J  - Chi_J * J  -  Chi_H * V_1 * H_1 ) * dt

    B ( 2 )  =  ( - Chi_H * H_1  +  V_1 * ( Xi_J  -  Chi_J * J ) ) * dt

  end subroutine ComputeCoefficients


  subroutine ComputeComovingIncrements ( A, B, dJ, dH_1 )

    real ( KDR ), dimension ( 2, 2 ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      B
    real ( KDR ), intent ( out ) :: &
      dJ, dH_1

    real ( KDR ) :: &
      Det

    Det   =  A ( 1, 1 )  *  A ( 2, 2 )   -   A ( 1, 2 )  *  A ( 2, 1 )

    dJ    =  ( B ( 1 )  *  A ( 2, 2 )   -   B ( 2 )  *  A ( 1, 2 ) )  /  Det
    dH_1  =  ( A ( 1, 1 )  *  B ( 2 )   -   A ( 2, 1 )  *  B ( 1 ) )  /  Det

  end subroutine ComputeComovingIncrements


  subroutine ComputeConservedIncrements &
               ( dE, dS_1, SVNE_E, SVNE_S_1, k_D, dJ, dH_1, V_1 )

    real ( KDR ), intent ( inout ) :: &
      dE, dS_1, &
      SVNE_E, SVNE_S_1
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      k_D
    real ( KDR ), intent ( in ) :: &
      dJ, dH_1, V_1

    dE    =  dE    +  dJ    +  2.0_KDR  *  V_1  *  dH_1
    dS_1  =  dS_1  +  dH_1  +  ( V_1  +  k_D ( 1 ) * V_1 )  *  dJ

    !-- FIXME: total hack, use this to accessibly store comoving increments 
    SVNE_E    =  dJ
    SVNE_S_1  =  dH_1

  end subroutine ComputeConservedIncrements


end module ApplyRelaxation_RM_2__Command
