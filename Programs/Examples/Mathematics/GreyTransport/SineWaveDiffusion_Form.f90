module SineWaveDiffusion_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use Interactions_F__Form
  use Interactions_ASC__Form
  use ApplyRelaxation_RM_G__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: SineWaveDiffusionForm
    type ( Interactions_ASC_Form ), allocatable :: &
      Interactions_ASC
    type ( RadiationMoments_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type SineWaveDiffusionForm

    private :: &
      SetRadiation, &
      SetInteractions, &
      SetReference

      private :: &
        SetRadiationKernel, &
        SetInteractionsKernel

    real ( KDR ), private :: &
      EquilibriumDensity, &
      EffectiveOpacity, &
      TransportOpacity, &
      TimeScale

contains


  subroutine Initialize ( SWD, Name )

    class ( SineWaveDiffusionForm ), intent ( inout ) :: &
      SWD
    character ( * ), intent ( in )  :: &
      Name

    class ( RadiationMomentsForm ), pointer :: &
      RM
    class ( Interactions_F_Form ), pointer :: &
      I

    EquilibriumDensity = 1.0_KDR
    EffectiveOpacity   = 0.0_KDR
    TransportOpacity   = 1.0e3_KDR
    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    TimeScale  =  3.0_KDR * ( TransportOpacity - EffectiveOpacity ) &
                  /  CONSTANT % SPEED_OF_LIGHT

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: SWD % PositionSpace )
    select type ( PS => SWD % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    call PS % CreateChart &
           ( MinCoordinateOption = [ -3.0_KDR, 0.0_KDR, 0.0_KDR ], &
             MaxCoordinateOption = [ +3.0_KDR, 0.0_KDR, 0.0_KDR ] )
    call PS % SetGeometry ( )

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: SWD % Current_ASC )
    select type ( RMA => SWD % Current_ASC )  !-- RadiationMomentsAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )

    !-- Interactions

    allocate ( SWD % Interactions_ASC )
    associate ( IA => SWD % Interactions_ASC )
    call IA % Initialize ( PS, 'FIXED' )
    call RMA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: SWD % Step )
    select type ( S => SWD % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( SWD, RMA, Name )
    S % ApplyRelaxation % Pointer => ApplyRelaxation_RM_G
    end select !-- S

    !-- Diagnostics

    allocate ( SWD % Reference )
    allocate ( SWD % Difference )
    call SWD % Reference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call SWD % Difference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    SWD % SetReference => SetReference

    !-- Initial conditions

    RM => RMA % RadiationMoments ( )
    call SetRadiation ( SWD, RM, Time = 0.0_KDR )

    I => IA % Interactions_F ( )
    call SetInteractions ( I )

    !-- Initialize template

    call SWD % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption  =  TimeScale )

    !-- Cleanup

    end associate !-- IA
    end select !-- RMA
    end select !-- PS
    nullify ( RM, I )

  end subroutine Initialize


  subroutine ComputeError ( SWD )

    class ( SineWaveDiffusionForm ), intent ( in ) :: &
      SWD

    real ( KDR ) :: &
      L1
    class ( RadiationMomentsForm ), pointer :: &
      RM
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => SWD % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    RM => SWD % Difference % RadiationMoments ( )

    associate &
      ( Difference => RM % Value ( :, RM % COMOVING_ENERGY ) )
    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )
    end associate !-- Difference

    associate &
      ( DifferenceSum => CO % Incoming % Value ( 1 ), &
        nValues => CO % Incoming % Value ( 2 ) )
    L1 = DifferenceSum / nValues
    end associate

    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end select !-- PS

  end subroutine ComputeError


  impure elemental subroutine Finalize ( SWD )
    
    type ( SineWaveDiffusionForm ), intent ( inout ) :: &
      SWD

    if ( allocated ( SWD % Difference ) ) &
      deallocate ( SWD % Difference )
    if ( allocated ( SWD % Reference ) ) &
      deallocate ( SWD % Reference )
    if ( allocated ( SWD % Interactions_ASC) ) &
      deallocate ( SWD % Interactions_ASC )

    call SWD % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetRadiation ( SWD, RM, Time )

    class ( SineWaveDiffusionForm ), intent ( in ) :: &
      SWD
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => SWD % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetRadiationKernel &
           ( X   = G % Value ( :, G % CENTER_U ( 1 ) ), &
             T   = Time, &
             Tau = TimeScale, &
             Pi  = CONSTANT % PI, &
             J   = RM % Value ( :, RM % COMOVING_ENERGY ), &
             HX  = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
             HY  = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
             HZ  = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ) )

    call RM % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetRadiation


  subroutine SetInteractions ( I )

    class ( Interactions_F_Form ), intent ( inout ) :: &
      I

    call SetInteractionsKernel &
           ( I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ) )

  end subroutine SetInteractions


  subroutine SetReference ( SWD )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      SWD

    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( SWD )
    class is ( SineWaveDiffusionForm )

    select type ( RMA => SWD % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- FA

    RM_R => SWD % Reference % RadiationMoments ( )
    call SetRadiation ( SWD, RM_R, SWD % Time )

    RM_D => SWD % Difference % RadiationMoments ( )
!    RM_D % Value  =  RM % Value  -  RM_R % Value
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- SWD
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


  subroutine SetRadiationKernel ( X, T, Tau, Pi, J, HX, HY, HZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X
    real ( KDR ), intent ( in ) :: &
      T, &
      Tau, &
      Pi
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues

    nValues = size ( X )
    
    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      J ( iV )  =  3.0_KDR  *  sqrt ( 4.0_KDR * Pi )  &
                   * ( 1.0_KDR  &
                       +  exp ( - Pi ** 2 / 9.0_KDR * T / Tau ) &
                        *  sin ( Pi * X ( iV ) / 3.0_KDR ) )

      HX ( iV )  =  - 1.0_KDR / Tau &
                      *  sqrt ( 4.0_KDR  *  Pi ** 3 )  &
                      *  exp ( - Pi ** 2 / 9.0_KDR  *  T / Tau ) &
                        *  cos ( Pi * X ( iV ) / 3.0_KDR )

      HY ( iV )  =  0.0_KDR
      HZ ( iV )  =  0.0_KDR

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


  subroutine SetInteractionsKernel ( Xi_J, Chi_J, Chi_H )

    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues  =  size ( Xi_J )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      Xi_J ( iV )   =  EffectiveOpacity * EquilibriumDensity
      Chi_J ( iV )  =  EffectiveOpacity
      Chi_H ( iV )  =  TransportOpacity
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetInteractionsKernel


end module SineWaveDiffusion_Form
