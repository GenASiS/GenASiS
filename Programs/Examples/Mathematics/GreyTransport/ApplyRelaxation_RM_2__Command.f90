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
      ComputeConservedIncrements, &
      ComputeSources, &
      ComputeDecoupled, &
      ComputeNumber

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
    class ( CurrentTemplate ), intent ( in ), target :: &
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
      iMomentum_1, &
      iNumber
    real ( KDR ) :: &
      dJ, dH_1, &
      Q_Old, A_1_Old, V_1_Old, Div_E
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

    associate ( I => RM % Interactions )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )

    select type ( Chart )
    class is ( Chart_SL_Template )

    G => Chart % Geometry ( )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, &
                  iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1 )

    !$OMP parallel do private &
    !$OMP   ( iV, A, B, dJ, dH_1, k_D, Q_Old, A_1_Old, V_1_Old, Div_E )
    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

      call ComputeCoefficients &
             ( RM, &
               IncrementExplicit % Value ( iV, iEnergy ), &
               IncrementExplicit % Value ( iV, iMomentum_1 ), &
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
               TimeStep, A, B, k_D, Q_Old, A_1_Old, V_1_Old, Div_E )

if ( iStage == 1 ) then
  SRM % Value ( iV, SRM % A_11 )  =  A ( 1, 1 )
  SRM % Value ( iV, SRM % A_21 )  =  A ( 2, 1 )
  SRM % Value ( iV, SRM % A_12 )  =  A ( 1, 2 )
  SRM % Value ( iV, SRM % A_22 )  =  A ( 2, 2 )
end if

!      if ( abs ( V_1_Old * A_1_Old )  <  abs ( Q_Old ) +  abs ( Div_E ) ) then
!      if ( A ( 1, 2 )  <  A ( 1, 1 ) ) then
        call ComputeComovingIncrements &
               ( A, B, dJ, dH_1 )
        call ComputeConservedIncrements &
               ( k_D, dJ, dH_1, &
                 RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
                 IncrementExplicit % Value ( iV, iEnergy ), &
                 IncrementExplicit % Value ( iV, iMomentum_1 ) )
!       else
! call Show ( RM % Name, '>>> Decoupling', CONSOLE % ERROR )
! call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
! call Show ( iV, '>>> iV', CONSOLE % ERROR )
!         call ComputeDecoupled &
!                ( IncrementExplicit % Value ( iV, iEnergy ), &
!                  IncrementExplicit % Value ( iV, iMomentum_1 ), &
!                  I  % Value ( iV, I % EMISSIVITY_J ), &
!                  I  % Value ( iV, I % OPACITY_J ), &
!                  I  % Value ( iV, I % OPACITY_H ), &
!                  RM % Value ( iV, RM % COMOVING_ENERGY ), &
!                  RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
!                  TimeStep, dJ, dH_1 )
!       end if

      call ComputeSources &
             ( SRM % Value ( iV, SRM % INTERACTIONS_J ), &
               SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 1 ) ), &
               I  % Value ( iV, I % EMISSIVITY_J ), &
               I  % Value ( iV, I % OPACITY_J ), &
               I  % Value ( iV, I % OPACITY_H ), &
               RM % Value ( iV, RM % COMOVING_ENERGY ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               dJ, dH_1, S % B ( iStage ) )

    end do
    !$OMP end parallel do

    end select !-- Chart
    end select !-- SRM
    end associate !-- I
    end select !-- RM

    nullify ( G )
    
  end subroutine ApplyRelaxation_RM_2


  subroutine ComputeCoefficients &
               ( RM, dE, dS_1, Xi_J, Chi_J, Chi_H, J, H_1, H_2, H_3, SF, V_1, &
                 M_DD_22, M_DD_33, dt, A, B, k_D, &
                 Q_Old, A_1_Old, V_1_Old, Div_E )

    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      dE, dS_1, &
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
    real ( KDR ), intent ( out ) :: &
      Q_Old, A_1_Old, V_1_Old, Div_E

    call RM % ComputeComovingStress_D &
           ( k_D, [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 1 )

    if ( J > 0.0_KDR ) then
      k_D  =  k_D / J
    else
      k_D  =  0.0_KDR
    end if

    A ( 1, 1 )  =  1.0_KDR  +  Chi_J * dt

    A ( 2, 1 )  =  V_1 * ( 1.0_KDR  +  Chi_J * dt )  +  k_D ( 1 ) * V_1

    A ( 1, 2 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_1

    A ( 2, 2 )  =  1.0_KDR  +  Chi_H * dt

    B ( 1 )  =  dE    +  (    Xi_J  -  Chi_J * J  &
                           -  Chi_H * V_1 * H_1 ) * dt

    B ( 2 )  =  dS_1  +  ( -  Chi_H * H_1  &
                           +  V_1 * ( Xi_J  -  Chi_J * J ) ) * dt

    Q_Old    =  Xi_J  -  Chi_J * J
    A_1_Old  =        -  Chi_H * H_1
    V_1_Old  =  V_1
    Div_E    =  dE

  end subroutine ComputeCoefficients


  ! subroutine ComputeComovingIncrements ( A, B, dJ, dH_1 )

  !   real ( KDR ), dimension ( 2, 2 ), intent ( in ) :: &
  !     A
  !   real ( KDR ), dimension ( 2 ), intent ( in ) :: &
  !     B
  !   real ( KDR ), intent ( out ) :: &
  !     dJ, dH_1

  !   real ( KDR ) :: &
  !     Det

  !   Det   =  A ( 1, 1 )  *  A ( 2, 2 )   -   A ( 1, 2 )  *  A ( 2, 1 )

  !   dJ    =  ( B ( 1 )  *  A ( 2, 2 )   -   B ( 2 )  *  A ( 1, 2 ) )  /  Det
  !   dH_1  =  ( A ( 1, 1 )  *  B ( 2 )   -   A ( 2, 1 )  *  B ( 1 ) )  /  Det

  ! end subroutine ComputeComovingIncrements


  subroutine ComputeComovingIncrements ( A, B, dJ, dH_1 )

    real ( KDR ), dimension ( 2, 2 ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( 2 ), intent ( in ) :: &
      B
    real ( KDR ), intent ( out ) :: &
      dJ, dH_1

    integer ( KDI ), dimension ( 2 ) :: &
      Pivot

    Pivot  =  maxloc ( abs ( A ) )

    if ( all ( Pivot == [ 1, 1 ] ) ) then

      dH_1  =  ( B ( 2 )  -  A ( 2, 1 ) * B ( 1 ) / A ( 1, 1 ) )  &
               /  ( A ( 2, 2 )  -  A ( 2, 1 ) * A ( 1, 2 ) / A ( 1, 1 ) )

      dJ    =  B ( 1 ) / A ( 1, 1 )  -  A ( 1, 2 ) / A ( 1, 1 ) * dH_1

    else if ( all ( Pivot == [ 1, 2 ] ) ) then

      dJ    =  ( B ( 2 )  -  A ( 2, 2 ) * B ( 1 ) / A ( 1, 2 ) )  &
               /  ( A ( 2, 1 )  -  A ( 2, 2 ) * A ( 1, 1 ) / A ( 1, 2 ) )

      dH_1  =  B ( 1 ) / A ( 1, 2 )  -  A ( 1, 1 ) / A ( 1, 2 ) * dJ

    else if ( all ( Pivot == [ 2, 1 ] ) ) then

      dH_1  =  ( B ( 1 )  -  A ( 1, 1 ) * B ( 2 ) / A ( 2, 1 ) )  &
               /  ( A ( 1, 2 )  -  A ( 1, 1 ) * A ( 2, 2 ) / A ( 2, 1 ) )

      dJ    =  B ( 2 ) / A ( 2, 1 )  -  A ( 2, 2 ) / A ( 2, 1 ) * dH_1

    else if ( all ( Pivot == [ 2, 2 ] ) ) then

      dJ    =  ( B ( 1 )  -  A ( 1, 2 ) * B ( 2 ) / A ( 2, 2 ) )  &
               /  ( A ( 1, 1 )  -  A ( 1, 2 ) * A ( 2, 1 ) / A ( 2, 2 ) )

      dH_1  =  B ( 2 ) / A ( 2, 2 )  -  A ( 2, 1 ) / A ( 2, 2 ) * dJ

    end if

  end subroutine ComputeComovingIncrements


  subroutine ComputeConservedIncrements ( k_D, dJ, dH_1, V_1, dE, dS_1 )

    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      k_D
    real ( KDR ), intent ( in ) :: &
      dJ, dH_1, V_1
    real ( KDR ), intent ( out ) :: &
      dE, dS_1

    dE    =  dJ    +  2.0_KDR  *  V_1  *  dH_1
    dS_1  =  dH_1  +  ( V_1  +  k_D ( 1 ) * V_1 )  *  dJ

  end subroutine ComputeConservedIncrements


  subroutine ComputeSources &
               ( Q, A_1, Xi_J, Chi_J, Chi_H, J, H_1, dJ, dH_1, Weight_RK )

    real ( KDR ), intent ( inout ) :: &
      Q, A_1
    real ( KDR ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_1, &
      dJ, dH_1, &
      Weight_RK

    Q    =  Q    +  Weight_RK * ( Xi_J  -  Chi_J  *  ( J    +  dJ ) )
    A_1  =  A_1  +  Weight_RK * (       -  Chi_H  *  ( H_1  +  dH_1 ) )

  end subroutine ComputeSources


  subroutine ComputeDecoupled &
               ( dE, dS_1, Xi_J, Chi_J, Chi_H, J, H_1, dt, dJ, dH_1 )

    real ( KDR ), intent ( inout ) :: &
      dE, dS_1
    real ( KDR ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_1, &
      dt
    real ( KDR ), intent ( out ) :: &
      dJ, dH_1

    dJ    =  ( dE  +  ( Xi_J  -  Chi_J * J ) * dt )  &
             /  ( 1.0_KDR  +  Chi_J * dt )

    dH_1  =  ( dS_1           -  Chi_H * H_1 * dt )  &
             /  ( 1.0_KDR  +  Chi_H * dt )

    dE    =  dJ
    dS_1  =  dH_1

  end subroutine ComputeDecoupled


  subroutine ComputeNumber &
               ( dD, R, Xi_N, Chi_N, J, H_1, N, V_1, dt, Weight_RK )

    real ( KDR ), intent ( inout ) :: &
      dD, &
      R
    real ( KDR ), intent ( in ) :: &
      Xi_N, Chi_N, &
      J, H_1, N, V_1, &
      dt, Weight_RK

    real ( KDR ) :: &
      dN
    real ( KDR ), dimension ( 3 ) :: &
      f_D

    !-- Increments

    f_D  =  0
    if ( J > 0.0_KDR ) then
      f_D ( 1 )  =  H_1 / J
    end if

    dN  =  ( dD  +  ( Xi_N  -  Chi_N * N ) * dt )  &
           /  ( 1.0_KDR  +  V_1 * f_D ( 1 )  +  Chi_N * dt )

    dD  =  ( 1.0_KDR  +  V_1 * f_D ( 1 ) )  *  dN

    !-- Source

    R  =  R  +  Weight_RK * ( Xi_N  -  Chi_N  *  ( N  +  dN ) )

  end subroutine ComputeNumber


end module ApplyRelaxation_RM_2__Command
