module DensityWave_Form

  use Basics
  use Manifolds
  use ProtoCurrent_Form
  use ProtoCurrent_ASC__Form

  implicit none
  private

  type, public :: DensityWaveForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iCycle = 0
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Time = 0.0_KDR
    character ( LDF ) :: &
      Name = ''
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
    type ( ProtoCurrent_ASC_Form ), allocatable :: &
      ProtoCurrent
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Write
    final :: &
      Finalize
    procedure, private, pass :: &
      Waveform
  end type DensityWaveForm

contains


  subroutine Initialize ( DW, Name )

    class ( DensityWaveForm ), intent ( inout ) :: &
      DW
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    logical ( KDL ) :: &
      VerboseStream
    character ( LDF ) :: &
      OutputDirectory
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( ProtoCurrentForm ), pointer :: &
      PC

    DW % IGNORABILITY = CONSOLE % INFO_1
    DW % Name = Name

    call Show ( 'Initializing a DensityWave', DW % IGNORABILITY )
    call Show ( DW % Name, 'Name', DW % IGNORABILITY )

    allocate ( DW % Atlas )
    associate ( A => DW % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart ( )

    allocate ( DW % Geometry )
    call DW % Geometry % Initialize ( A )
    call A % SetGeometry ( DW % Geometry )

    allocate ( DW % ProtoCurrent )
    associate ( PCA => DW % ProtoCurrent )
    call PCA % Initialize ( A )

    nWavelengths = 0
    nWavelengths ( 1 : A % nDimensions ) = 1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    DW % Offset    = 2.0_KDR
    DW % Amplitude = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( DW % Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( DW % Amplitude, 'Amplitude' )

    G  => A % Geometry ( )
    PC => PCA % ProtoCurrent ( )

    PC % Speed = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( PC % Speed, 'Speed' )

    associate ( C => A % Chart )
    associate ( BoxSize => C % MaxCoordinate - C % MinCoordinate )
    where ( BoxSize > 0.0_KDR )
      PC % Wavenumber = nWavelengths / BoxSize
    elsewhere
      PC % Wavenumber = 0.0_KDR
    end where
    end associate !-- BoxSize
    end associate !-- C

    associate &
      ( X  => G % Value ( :, G % CENTER ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER ( 3 ) ), &
        N  => PC % Value ( :, PC % COMOVING_DENSITY ), &
        K  => PC % Wavenumber )

    N = DW % Waveform ( K ( 1 ) * X  +  K ( 2 ) * Y  +  K ( 3 ) * Z )

    call PC % ComputeFromPrimitive ( G )

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )

    allocate ( DW % GridImageStream )
    associate ( GIS => DW % GridImageStream )
    call GIS % Initialize &
           ( DW % Name, CommunicatorOption = A % Communicator, &
             WorkingDirectoryOption = OutputDirectory )

    VerboseStream = .false.
    call PROGRAM_HEADER % GetParameter ( VerboseStream, 'VerboseStream' )

    call A % OpenStream &
           ( GIS, 'Time', iStream = 1, VerboseOption = VerboseStream )

    end associate !-- GIS
    end associate !-- X, etc.
    end associate !-- PCA
    end associate !-- A

  end subroutine Initialize


  subroutine Write ( DW )

    class ( DensityWaveForm ), intent ( inout ) :: &
      DW

    associate &
      ( A => DW % Atlas, &
        GIS => DW % GridImageStream )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write &
           ( iStream = 1, TimeOption = DW % Time * UNIT % IDENTITY, &
             CycleNumberOption = DW % iCycle )
    call GIS % Close ( )

    end associate !-- A, etc.

  end subroutine Write


  subroutine Finalize ( DW )

    type ( DensityWaveForm ), intent ( inout ) :: &
      DW

    if ( allocated ( DW % ProtoCurrent ) ) &
      deallocate ( DW % ProtoCurrent )
    if ( allocated ( DW % Atlas ) ) &
      deallocate ( DW % Atlas )
    if ( allocated ( DW % Geometry ) ) &
      deallocate ( DW % Geometry )
    if ( allocated ( DW % GridImageStream ) ) &
      deallocate ( DW % GridImageStream )

    call Show ( 'Finalizing a DensityWave', DW % IGNORABILITY )
    call Show ( DW % Name, 'Name', DW % IGNORABILITY )

  end subroutine Finalize


  elemental function Waveform ( DW, X ) result ( W )

    !-- Waveform with a full period in the range 0 < X < 1

    class ( DensityWaveForm ), intent ( in ) :: &
      DW
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ) :: &
      W
    
    associate ( TwoPi => 2.0_KDR * CONSTANT % PI )

    W = DW % Offset + DW % Amplitude * sin ( TwoPi * X )

    end associate !-- TwoPi

  end function Waveform


end module DensityWave_Form
