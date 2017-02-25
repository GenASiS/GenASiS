module DensityWaveIntegrator_Form

  use Basics
  use Manifolds
  use Steps
  use Integrator_Template
  use Integrator_C_PS__Template
  use ProtoFields

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: &
    DensityWaveIntegratorForm
      real ( KDR ) :: &
        Offset, &
        Amplitude, &
        Speed
      real ( KDR ), dimension ( 3 ) :: &
        Wavenumber
      type ( ProtoCurrent_ASC_Form ), allocatable :: &
        Reference, &
        Difference
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
    procedure, private, pass :: &
      Waveform
  end type DensityWaveIntegratorForm

    private :: &
      SetCurrent, &
      SetReference

contains


  subroutine Initialize ( DW, Name )

    class ( DensityWaveIntegratorForm ), intent ( inout ) :: &
      DW
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      nPeriods
    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Period
    class ( ProtoCurrentForm ), pointer :: &
      PC

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: DW % PositionSpace )
    select type ( PS => DW % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( DW % Geometry_ASC )
    associate ( GA => DW % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Conserved current

    allocate ( ProtoCurrent_ASC_Form :: DW % Current_ASC )
    select type ( PCA => DW % Current_ASC )
    class is ( ProtoCurrent_ASC_Form )
    call PCA % Initialize ( PS )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: DW % Step )
    select type ( S => DW % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( Name )
    end select !-- Step

    !-- Diagnostics

    allocate ( DW % Reference )
    allocate ( DW % Difference )
    call DW % Reference % Initialize ( PS, NameOutputOption = 'Reference' )
    call DW % Difference % Initialize ( PS, NameOutputOption = 'Difference' )
    DW % SetReference => SetReference

    !-- Initial conditions

    nWavelengths = 0
    nWavelengths ( 1 : PS % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    associate ( C => PS % Chart )
    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      DW % Wavenumber = nWavelengths / BoxSize
    elsewhere
      DW % Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize
    end associate !-- C

    DW % Offset    = 2.0_KDR
    DW % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( DW % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( DW % Amplitude, 'Amplitude' )

    DW % Speed = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( DW % Speed, 'Speed' )

    nPeriods = 1
    call PROGRAM_HEADER % GetParameter ( nPeriods, 'nPeriods' )

    associate &
      ( K     => DW % Wavenumber, &
        Abs_K => sqrt ( dot_product ( DW % Wavenumber, DW % Wavenumber ) ), &
        V     => DW % Speed )
    Period = 1.0_KDR / ( Abs_K * V )
    call Show ( Period, 'Period' )
    end associate !-- K, etc.

    PC => PCA % ProtoCurrent ( )
    call SetCurrent ( DW, PC, Time = 0.0_KDR )

    !-- Initialize template

    call DW % InitializeTemplate_C_PS &
           ( Name, FinishTimeOption = nPeriods * Period )

    !-- Cleanup

    end select !-- PCA
    end select !-- PS

    nullify ( PC )

  end subroutine Initialize


  subroutine ComputeError ( DW )

    class ( DensityWaveIntegratorForm ), intent ( in ) :: &
      DW

    real ( KDR ) :: &
      L1
    class ( ProtoCurrentForm ), pointer :: &
      PC_D
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => DW % PositionSpace )
    class is ( Atlas_SC_Form )
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    PC_D => DW % Difference % ProtoCurrent ( )

    associate ( D => PC_D % Value ( :, PC_D % COMOVING_DENSITY ) )
    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( D ) )
    CO % Outgoing % Value ( 2 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )
    end associate !-- D

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


  impure elemental subroutine Finalize ( DW )

    type ( DensityWaveIntegratorForm ), intent ( inout ) :: &
      DW

    if ( allocated ( DW % Difference ) ) &
      deallocate ( DW % Difference )
    if ( allocated ( DW % Reference ) ) &
      deallocate ( DW % Reference )

    call DW % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  elemental function Waveform ( DW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( DensityWaveIntegratorForm ), intent ( in ) :: &
      DW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    associate ( TwoPi => 2.0_KDR * CONSTANT % PI )

    W = DW % Offset + DW % Amplitude * sin ( TwoPi * X )

    end associate !-- TwoPi

  end function Waveform


  subroutine SetCurrent ( DW, PC, Time )

    class ( DensityWaveIntegratorForm ), intent ( in ) :: &
      DW
    class ( ProtoCurrentForm ), intent ( inout ) :: &
      PC
    real ( KDR ), intent ( in ) :: &
      Time
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => DW % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    PC % Speed      = DW % Speed
    PC % Wavenumber = DW % Wavenumber
    
    associate &
      ( K     => PC % Wavenumber, &
        Abs_K => sqrt ( dot_product ( PC % Wavenumber, PC % Wavenumber ) ), &
        V     => PC % Speed, &
        T     => Time, &
        X     => G % Value ( :, G % CENTER ( 1 ) ), &
        Y     => G % Value ( :, G % CENTER ( 2 ) ), &
        Z     => G % Value ( :, G % CENTER ( 3 ) ), &
        VX    => PC % Value ( :, PC % VELOCITY ( 1 ) ), &
        VY    => PC % Value ( :, PC % VELOCITY ( 2 ) ), &
        VZ    => PC % Value ( :, PC % VELOCITY ( 3 ) ), &
        N     => PC % Value ( :, PC % COMOVING_DENSITY ) )

    VX = V * K ( 1 ) / Abs_K
    VY = V * K ( 2 ) / Abs_K
    VZ = V * K ( 3 ) / Abs_K

    N = DW % Waveform &
          (    K ( 1 ) * ( X  -  VX * T ) &
            +  K ( 2 ) * ( Y  -  VY * T ) &
            +  K ( 3 ) * ( Z  -  VZ * T ) )

    call PC % ComputeFromPrimitive ( G )

    end associate !-- K, etc.
    end select    !-- PS

    nullify ( G )

  end subroutine SetCurrent


  subroutine SetReference ( DW )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      DW

    class ( ProtoCurrentForm ), pointer :: &
      PC, &
      PC_R, &  !-- PC_Reference
      PC_D     !-- PC_Difference

    select type ( DW )
    class is ( DensityWaveIntegratorForm )

    select type ( PCA => DW % Current_ASC )
    class is ( ProtoCurrent_ASC_Form )
    PC => PCA % ProtoCurrent ( )
    end select !-- PCA

    PC_R => DW % Reference % ProtoCurrent ( )
    call SetCurrent ( DW, PC_R, DW % Time )

    PC_D => DW % Difference % ProtoCurrent ( )

    PC_D % Value  =  PC % Value  -  PC_R % Value

    end select !-- DW
    nullify ( PC, PC_R, PC_D )

  end subroutine SetReference


end module DensityWaveIntegrator_Form
