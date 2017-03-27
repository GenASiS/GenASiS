module PlaneWaveStreaming_Template

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    PlaneWaveStreamingTemplate
      real ( KDR ) :: &
        Speed
      real ( KDR ), dimension ( 3 ) :: &
        Wavenumber
      type ( RadiationMoments_ASC_Form ), allocatable :: &
        Reference, &
        Difference
  contains
    procedure, public, pass :: &
      InitializeTemplate_PWS
    procedure, public, pass :: &
      ComputeError
    procedure, public, pass :: &
      FinalizeTemplate_PWS
    procedure ( Waveform ), private, pass, deferred :: &
      Waveform
  end type PlaneWaveStreamingTemplate

  abstract interface
    elemental function Waveform ( PWS, X ) result ( W )
      !-- Waveform with a full period in the range 0 < X < 1
      use Basics
      import PlaneWaveStreamingTemplate
      class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
        PWS
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ) :: &
        W
    end function Waveform
  end interface

    private :: &
      SetRadiation, &
      SetReference

      private :: &
        SetRadiationKernel

contains


  subroutine InitializeTemplate_PWS ( PWS, Name )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ), target :: &
      PWS
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    class ( RadiationMomentsForm ), pointer :: &
      RM

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: PWS % PositionSpace )
    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( PWS % Geometry_ASC )
    associate ( GA => PWS % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: PWS % Current_ASC )
    select type ( RMA => PWS % Current_ASC )  !-- FluidAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: PWS % Step )
    select type ( S => PWS % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( RMA, Name )
    end select !-- S

    !-- Diagnostics

    allocate ( PWS % Reference )
    allocate ( PWS % Difference )
    call PWS % Reference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Reference', &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call PWS % Difference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Difference', &
             IgnorabilityOption = CONSOLE % INFO_2 )
    PWS % SetReference => SetReference

    !-- Initial conditions

    nWavelengths = 0
    nWavelengths ( 1 : PS % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( C => PS % Chart )
    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      PWS % Wavenumber = nWavelengths / BoxSize
    elsewhere
      PWS % Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize
    end associate !-- C

    PWS % Speed = CONSTANT % SPEED_OF_LIGHT

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )

    associate &
      ( K     => PWS % Wavenumber, &
        Abs_K => sqrt ( dot_product &
                          ( PWS % Wavenumber, PWS % Wavenumber ) ), &
        V     => PWS % Speed )
    Period = 1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )
    end associate !-- K, etc.

    RM => RMA % RadiationMoments ( )
    call SetRadiation ( PWS, RM, Time = 0.0_KDR )

    !-- Initialize template

    call PWS % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = nPeriods * Period )

    !-- Cleanup

    end select !-- RMA
    end select !-- PS
    nullify ( RM )

  end subroutine InitializeTemplate_PWS


  subroutine ComputeError ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS

    real ( KDR ) :: &
      L1
    class ( RadiationMomentsForm ), pointer :: &
      RM
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    RM => PWS % Difference % RadiationMoments ( )

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


  impure elemental subroutine FinalizeTemplate_PWS ( PWS )

    class ( PlaneWaveStreamingTemplate ), intent ( inout ) :: &
      PWS

    if ( allocated ( PWS % Difference ) ) &
      deallocate ( PWS % Difference )
    if ( allocated ( PWS % Reference ) ) &
      deallocate ( PWS % Reference )

    call PWS % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_PWS


  subroutine SetRadiation ( PWS, RM, Time )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => PWS % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    call SetRadiationKernel &
           ( PWS, &
             X  = G % Value ( :, G % CENTER ( 1 ) ), &
             Y  = G % Value ( :, G % CENTER ( 2 ) ), &
             Z  = G % Value ( :, G % CENTER ( 3 ) ), &
             K  = PWS % Wavenumber, &
             V  = PWS % Speed, &
             T  = Time, &
             J  = RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
             HX = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
             HY = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
             HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ) )

    call RM % ComputeFromPrimitive ( G )

    end select    !-- PS
    nullify ( G )

  end subroutine SetRadiation


  subroutine SetReference ( PWS )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      PWS

    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( PWS )
    class is ( PlaneWaveStreamingTemplate )

    select type ( RMA => PWS % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- FA

    RM_R => PWS % Reference % RadiationMoments ( )
    call SetRadiation ( PWS, RM_R, PWS % Time )

    RM_D => PWS % Difference % RadiationMoments ( )
!    RM_D % Value  =  RM % Value  -  RM_R % Value
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- PWS
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


  subroutine SetRadiationKernel ( PWS, X, Y, Z, K, V, T, J, HX, HY, HZ )

    class ( PlaneWaveStreamingTemplate ), intent ( in ) :: &
      PWS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      K
    real ( KDR ), intent ( in ) :: &
      V, &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      Abs_K, &
      VX, VY, VZ

    nValues = size ( X )
    
    Abs_K  =  sqrt ( dot_product ( K, K ) )
       VX  =  V * K ( 1 ) / Abs_K
       VY  =  V * K ( 2 ) / Abs_K
       VZ  =  V * K ( 3 ) / Abs_K

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      J ( iV )  =  PWS % Waveform &
                     (    K ( 1 ) * ( X ( iV )  -  VX * T ) &
                       +  K ( 2 ) * ( Y ( iV )  -  VY * T ) &
                       +  K ( 3 ) * ( Z ( iV )  -  VZ * T ) )
    end do !-- iV
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      HX ( iV )  =  VX  *  J ( iV )
      HY ( iV )  =  VY  *  J ( iV )
      HZ ( iV )  =  VZ  *  J ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


end module PlaneWaveStreaming_Template

