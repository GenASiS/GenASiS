module ApplyRelaxation_RM_G__Command

  !-- ApplyRelaxation_RadiationMoments_Grey__Command

  use OMP_LIB
  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    ApplyRelaxation_RM_G

    private :: &
      ComputeCoefficients, &
      ComputeComovingIncrements, &
      ComputeConservedIncrements, &
      ComputeSources

contains


  subroutine ApplyRelaxation_RM_G &
               ( S, Sources_RM, Increment, RadiationMoments, Chart, &
                 TimeStep, iStage, GeometryOption, iStrgeometryValueOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( StorageForm ), intent ( inout ) :: &
      Increment
    class ( CurrentTemplate ), intent ( in ), target :: &
      RadiationMoments
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage
    class ( GeometryFlatForm ), intent ( in ), optional :: &
      GeometryOption
    integer ( KDI ), intent ( in ), optional :: &
      iStrgeometryValueOption

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV, &  !-- nValues
      iThread, &
      iEnergy, &
      iMomentum_1, iMomentum_2, iMomentum_3 !&
!      iNumber
    real ( KDR ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3
    real ( KDR ), dimension ( 3, 3 ) :: &
      k_DD  !-- Eddington tensor components ( K / J, Stress / EnergyDensity )
    class ( GeometryFlatForm ), pointer :: &
      G
!    class ( NeutrinoMoments_G_Form ), pointer :: &
!      NM
    type ( LinearEquations_LAPACK_Form ), dimension ( : ), allocatable :: &
      LE

    call Show ( 'ApplyRelaxation_RM_G', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    nV = size ( Increment % Value, dim = 1 )

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
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 2 ), &
                  iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 3 ), &
                  iMomentum_3 )
!    if ( associated ( NM ) ) &
!      call Search ( NM % iaConserved, NM % CONSERVED_NUMBER, iNumber )

    allocate ( LE ( 0 : PROGRAM_HEADER % MaxThreads - 1 ) )

    !$OMP parallel do private &
    !$OMP ( iV, dJ, dH_D_1, dH_D_2, dH_D_3, k_DD )
    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

      iThread = OMP_GET_THREAD_NUM ( )
      if ( .not. allocated ( LE ( iThread ) % Matrix ) ) &
        call LE ( iThread ) % Initialize ( nEquations = 4, nSolutions = 1 )

      call ComputeCoefficients &
             ( RM, &
               Increment % Value ( iV, iEnergy ), &
               Increment % Value ( iV, iMomentum_1 ), &
               Increment % Value ( iV, iMomentum_2 ), &
               Increment % Value ( iV, iMomentum_3 ), &
               I  % Value ( iV, I % EMISSIVITY_J ), &
               I  % Value ( iV, I % OPACITY_J ), &
               I  % Value ( iV, I % OPACITY_H ), &
               RM % Value ( iV, RM % COMOVING_ENERGY ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               RM % Value ( iV, RM % STRESS_FACTOR ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 2 ) ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 3 ) ), &
               G  % Value ( iV, G % METRIC_DD_22 ), &
               G  % Value ( iV, G % METRIC_DD_33 ), &
               TimeStep, LE ( iThread ) % Matrix, &
               LE ( iThread ) % RightHandSide, k_DD )

      call ComputeComovingIncrements &
             ( LE ( iThread ), dJ, dH_D_1, dH_D_2, dH_D_3 )

      call ComputeConservedIncrements &
             ( k_DD, dJ, dH_D_1, dH_D_2, dH_D_3, &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 2 ) ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 3 ) ), &
               G  % Value ( iV, G % METRIC_DD_22 ), &
               G  % Value ( iV, G % METRIC_DD_33 ), &
               Increment % Value ( iV, iEnergy ), &
               Increment % Value ( iV, iMomentum_1 ), &
               Increment % Value ( iV, iMomentum_2 ), &
               Increment % Value ( iV, iMomentum_3 ) )

      call ComputeSources &
             ( SRM % Value ( iV, SRM % INTERACTIONS_J ), &
               SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 1 ) ), &
               SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 2 ) ), &
               SRM % Value ( iV, SRM % INTERACTIONS_H_D ( 3 ) ), &
                 I % Value ( iV, I % EMISSIVITY_J ), &
                 I % Value ( iV, I % OPACITY_J ), &
                 I % Value ( iV, I % OPACITY_H ), &
                RM % Value ( iV, RM % COMOVING_ENERGY ), &
                RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
                RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
                RM % Value ( iV, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
                dJ, dH_D_1, dH_D_2, dH_D_3, & 
                 G % Value ( iV, G % METRIC_DD_22 ), &
                 G % Value ( iV, G % METRIC_DD_33 ), &
                S % B ( iStage ) )

    end do
    !$OMP end parallel do

    end select !-- Chart
    end select !-- SRM
    end associate !-- I
    end select !-- RM

    nullify ( G )!, NM )
    
  end subroutine ApplyRelaxation_RM_G


  subroutine ComputeCoefficients &
               ( RM, dE, dS_D_1, dS_D_2, dS_D_3, Xi_J, Chi_J, Chi_H, &
                 J, H_1, H_2, H_3, SF, V_1, V_2, V_3, M_DD_22, M_DD_33, dt, &
                 A, B, k_DD )

    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      dE, dS_D_1, dS_D_2, dS_D_3, &
      Xi_J, Chi_J, Chi_H, &
      J, H_1, H_2, H_3, SF, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33, &
      dt
    real ( KDR ), dimension ( 4, 4 ), intent ( out ) :: &
      A
    real ( KDR ), dimension ( 4, 1 ), intent ( out ) :: &
      B
    real ( KDR ), dimension ( 3, 3 ), intent ( out ) :: &
      k_DD
    
    call RM % ComputeComovingStress_D &
           ( k_DD ( 1, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 1 )
    call RM % ComputeComovingStress_D &
           ( k_DD ( 2, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 2 )
    call RM % ComputeComovingStress_D &
           ( k_DD ( 3, : ), [ H_1, H_2, H_3 ], J, SF, M_DD_22, M_DD_33, iD = 3 )

    if ( J > 0.0_KDR ) then
      k_DD  =  k_DD / J
    else
      k_DD  =  0.0_KDR
    end if

    A ( 1, 1 )  =  1.0_KDR  +  Chi_J * dt
    A ( 2, 1 )  =            V_1 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 1, : ), [ V_1, V_2, V_3 ] )
    A ( 3, 1 )  =  M_DD_22 * V_2 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 2, : ), [ V_1, V_2, V_3 ] )
    A ( 4, 1 )  =  M_DD_33 * V_3 * ( 1.0_KDR  +  Chi_J * dt )  &
                   +  dot_product ( k_DD ( 3, : ), [ V_1, V_2, V_3 ] )

    A ( 1, 2 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_1
    A ( 2, 2 )  =  1.0_KDR  +  Chi_H * dt
    A ( 3, 2 )  =  0.0_KDR
    A ( 4, 2 )  =  0.0_KDR

    A ( 1, 3 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_2
    A ( 2, 3 )  =  0.0_KDR
    A ( 3, 3 )  =  1.0_KDR  +  Chi_H * dt
    A ( 4, 3 )  =  0.0_KDR

    A ( 1, 4 )  =  ( 2.0_KDR  +  Chi_H * dt )  *  V_3
    A ( 2, 4 )  =  0.0_KDR
    A ( 3, 4 )  =  0.0_KDR
    A ( 4, 4 )  =  1.0_KDR  +  Chi_H * dt

    B ( 1, 1 )  =  dE    +  (    Xi_J  -  Chi_J * J  &
                              -  Chi_H * (             V_1 * H_1  &
                                           * M_DD_22 * V_2 * H_2  &
                                           * M_DD_33 * V_3 * H_3 ) ) * dt
    B ( 2, 1 )  =  dS_D_1  +            ( -  Chi_H * H_1  &
                                        +  V_1 * ( Xi_J  -  Chi_J * J ) ) * dt
    B ( 3, 1 )  =  dS_D_2  +  M_DD_22 * ( -  Chi_H * H_2  &
                                        +  V_2 * ( Xi_J  -  Chi_J * J ) ) * dt
    B ( 4, 1 )  =  dS_D_3  +  M_DD_33 * ( -  Chi_H * H_3  &
                                        +  V_3 * ( Xi_J  -  Chi_J * J ) ) * dt

  end subroutine ComputeCoefficients


  subroutine ComputeComovingIncrements ( LE, dJ, dH_D_1, dH_D_2, dH_D_3 )

    type ( LinearEquations_LAPACK_Form ), intent ( inout ) :: &
      LE
    real ( KDR ), intent ( out ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3

    call LE % Solve ( )

    dJ      =  LE % Solution ( 1, 1 )
    dH_D_1  =  LE % Solution ( 2, 1 )
    dH_D_2  =  LE % Solution ( 3, 1 )
    dH_D_3  =  LE % Solution ( 4, 1 )

  end subroutine ComputeComovingIncrements


  subroutine ComputeConservedIncrements &
               ( k_DD, dJ, dH_D_1, dH_D_2, dH_D_3, V_1, V_2, V_3, &
                 M_DD_22, M_DD_33, dE, dS_D_1, dS_D_2, dS_D_3 )

    real ( KDR ), dimension ( 3, 3 ), intent ( in ) :: &
      k_DD
    real ( KDR ), intent ( in ) :: &
      dJ, dH_D_1, dH_D_2, dH_D_3, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33
    real ( KDR ), intent ( out ) :: &
      dE, dS_D_1, dS_D_2, dS_D_3

    dE    =  dJ    +  2.0_KDR  *  dot_product ( [ V_1, V_2, V_3 ], &
                                                [ dH_D_1, dH_D_2, dH_D_3 ] )
    dS_D_1  =  dH_D_1  &
             +  ( V_1  &
                  +  dot_product ( k_DD ( 1, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

    dS_D_2  =  dH_D_2  &
             +  ( M_DD_22 * V_2  &
                  +  dot_product ( k_DD ( 2, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

    dS_D_3  =  dH_D_3  &
             +  ( M_DD_33 * V_3  &
                  +  dot_product ( k_DD ( 3, : ), [ V_1, V_2, V_3 ] ) )  *  dJ  

  end subroutine ComputeConservedIncrements


  subroutine ComputeSources &
               ( Q, A_D_1, A_D_2, A_D_3, Xi_J, Chi_J, Chi_H, &
                 J, H_U_1, H_U_2, H_U_3, dJ, dH_D_1, dH_D_2, dH_D_3, &
                 M_DD_22, M_DD_33, Weight_RK )

    real ( KDR ), intent ( inout ) :: &
      Q, A_D_1, A_D_2, A_D_3
    real ( KDR ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_U_1, H_U_2, H_U_3, &
      dJ, dH_D_1, dH_D_2, dH_D_3, &
      M_DD_22, M_DD_33, &
      Weight_RK

    Q      =  Q      &
              +  Weight_RK * ( Xi_J  &
                               -  Chi_J  *  (               J  +      dJ ) )
    A_D_1  =  A_D_1  &
              +  Weight_RK * (  &
                               -  Chi_H  *  (           H_U_1  +  dH_D_1 ) )
    A_D_2  =  A_D_2  &
              +  Weight_RK * (  &
                               -  Chi_H  *  ( M_DD_22 * H_U_2  +  dH_D_2 ) )
    A_D_3  =  A_D_3  &
              +  Weight_RK * (  &
                               -  Chi_H  *  ( M_DD_33 * H_U_3  +  dH_D_3 ) )

  end subroutine ComputeSources


end module ApplyRelaxation_RM_G__Command
