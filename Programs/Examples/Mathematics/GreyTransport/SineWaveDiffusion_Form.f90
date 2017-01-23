module SineWaveDiffusion_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use Interactions_F__Form
  use Interactions_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_Template ) :: SineWaveDiffusionForm
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
      TransportOpacity

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

    if ( SWD % Type == '' ) &
      SWD % Type = 'a SineWaveDiffusion' 

    EquilibriumDensity = 1.0_KDR
    EffectiveOpacity   = 0.0_KDR
    TransportOpacity   = 1.0e3_KDR

    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: SWD % PositionSpace )
    select type ( PS => SWD % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart &
           ( MinCoordinateOption = [ -3.0_KDR, 0.0_KDR, 0.0_KDR ], &
             MaxCoordinateOption = [ +3.0_KDR, 0.0_KDR, 0.0_KDR ] )

    !-- Geometry of PositionSpace

    allocate ( SWD % Geometry_ASC )
    associate ( GA => SWD % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: SWD % Current_ASC )
    select type ( RMA => SWD % Current_ASC )  !-- FluidAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )

    !-- Interactions

    allocate ( SWD % Interactions_ASC )
    associate ( IA => SWD % Interactions_ASC )
    call IA % Initialize ( PS, 'FIXED' )
    call RMA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_Form :: SWD % Step )
    select type ( S => SWD % Step )
    class is ( Step_RK2_C_Form )
    call S % Initialize ( Name )
    S % ApplyRelaxation  =>  ApplyRelaxation_Interactions
    end select !-- S

    !-- Diagnostics

    allocate ( SWD % Reference )
    allocate ( SWD % Difference )
    call SWD % Reference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Reference' )
    call SWD % Difference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Difference' )
    SWD % SetReference => SetReference

    !-- Initial conditions

    RM => RMA % RadiationMoments ( )
    call SetRadiation ( SWD, RM, Time = 0.0_KDR )

    I => IA % Interactions_F ( )
    call SetInteractions ( I )

    !-- Initialize template

    call SWD % InitializeTemplate_C &
           ( Name, FinishTimeOption  =  3.0_KDR * TransportOpacity )

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
      ( Difference => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ) )
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

    call SWD % FinalizeTemplate_C ( )

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
           ( X  = G % Value ( :, G % CENTER ( 1 ) ), &
             T  = Time, &
             TO = TransportOpacity, &
             Pi = CONSTANT % PI, &
             J  = RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
             HX = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
             HY = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
             HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ) )

    call RM % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetRadiation


  subroutine SetInteractions ( I )

    class ( Interactions_F_Form ), intent ( inout ) :: &
      I

    call SetInteractionsKernel &
           ( I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
             I % Value ( :, I % EFFECTIVE_OPACITY ), &
             I % Value ( :, I % TRANSPORT_OPACITY ) )

  end subroutine SetInteractions


  subroutine SetReference ( SWD )

    class ( IntegratorTemplate ), intent ( in ) :: &
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


  subroutine SetRadiationKernel ( X, T, TO, Pi, J, HX, HY, HZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X
    real ( KDR ), intent ( in ) :: &
      T, &
      TO, &
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
                       +  exp ( - Pi ** 2 / 9.0_KDR * T / ( 3.0_KDR * TO ) ) &
                        *  sin ( Pi * X ( iV ) / 3.0_KDR ) )

      HX ( iV )  =  - 1.0_KDR / ( 3.0_KDR * TO ) &
                      *  sqrt ( 4.0_KDR  *  Pi ** 3 )  &
                      *  exp ( - Pi ** 2 / 9.0_KDR  *  T / ( 3.0_KDR * TO ) ) &
                        *  cos ( Pi * X ( iV ) / 3.0_KDR )

      HY ( iV )  =  0.0_KDR
      HZ ( iV )  =  0.0_KDR

    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


  subroutine SetInteractionsKernel ( EDV, EOV, TOV )

    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EDV, &
      EOV, &
      TOV

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues  =  size ( EDV )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      EDV ( iV )  =  EquilibriumDensity
      EOV ( iV )  =  EffectiveOpacity
      TOV ( iV )  =  TransportOpacity
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetInteractionsKernel


end module SineWaveDiffusion_Form
