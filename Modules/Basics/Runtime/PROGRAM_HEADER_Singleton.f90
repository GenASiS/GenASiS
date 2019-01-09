!-- PROGRAM_HEADER provides functionalities commonly required by all programs
!   (drivers), including initialization of parallel environment, obtaining
!   program parameters, and displaying basic runtime statistics.

module PROGRAM_HEADER_Singleton

  use ISO_FORTRAN_ENV
  use OMP_LIB
  use Specifiers
  use Devices
  use Display
  use MessagePassing
  use FileSystem
  use Timer_Form
  use GetMemoryUsage_Command
  use CommandLineOptions_Form
  !  use petsc

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_TIMERS = 128

  type, public :: ProgramHeaderSingleton
    integer ( KDI ) :: &
      nTimers = 0, &
      MaxThreads = 0, &
      TimerLevel = 0, &
      ExecutionTimeHandle
    real ( KDR ) :: &
      TimerDisplayFraction
    character ( LDL ) :: &
      Dimensionality = ''
    character ( LDF ) :: &
      Name
    type ( CommunicatorForm ), allocatable :: &
      Communicator 
    type ( ParameterStreamForm ), allocatable :: &
      ParameterStream
    type ( TimerForm ), dimension ( : ), allocatable :: &
      Timer
    type ( CommandLineOptionsForm ), allocatable :: &
      CommandLineOptions
  contains
    procedure, public, nopass :: &
      Initialize
    procedure, private, nopass :: &
      GetParameter_0D_Integer
    procedure, private, nopass :: &
      GetParameter_0D_Real
    procedure, private, nopass :: &
      GetParameter_0D_MeasuredValue
    procedure, private, nopass :: &
      GetParameter_0D_Logical
    procedure, private, nopass :: &
      GetParameter_0D_Character
    procedure, private, nopass :: &
      GetParameter_1D_Integer
    procedure, private, nopass :: &
      GetParameter_1D_Real
    procedure, private, nopass :: &
      GetParameter_1D_MeasuredValue
    procedure, private, nopass :: &
      GetParameter_1D_Logical
    procedure, private, nopass :: &
      GetParameter_1D_Character
    generic :: &
      GetParameter &
        => GetParameter_0D_Integer, GetParameter_0D_Real, &
           GetParameter_0D_MeasuredValue, GetParameter_0D_Logical, &
           GetParameter_0D_Character, &
           GetParameter_1D_Integer, GetParameter_1D_Real, &
           GetParameter_1D_MeasuredValue, GetParameter_1D_Logical, &
           GetParameter_1D_Character
    procedure, public, nopass :: &
      AddTimer
    procedure, public, nopass :: &
      TimerPointer
    procedure, public, nopass :: &
      ShowStatistics
    procedure, public, nopass :: &
      Abort => Abort_PH  !-- avoids conflict with intrinsic "abort"
    final :: &
      Finalize
  end type ProgramHeaderSingleton
  
  
  type ( ProgramHeaderSingleton ), public, target, allocatable :: &
    PROGRAM_HEADER

    private :: &
      PrepareAndShow_OMP_Environment, &
      ReadTimers
      
contains

 
  subroutine Initialize &
               ( Name, DimensionalityOption, AppendDimensionalityOption ) 

    character ( * ), intent ( in )  :: &
      Name
    character ( * ), intent ( in ), optional :: &
      DimensionalityOption
    logical ( KDL ), intent ( in ), optional :: &
      AppendDimensionalityOption
      
    integer ( KDI )  :: &
      DisplayRank, &
      OMP_ScheduleChunkSize
    integer ( OMP_SCHED_KIND ) :: &
      OMP_ScheduleKind
    character ( 5 ) :: &
      Encoding
    character ( LDL )  :: &
      Verbosity, &
      OMP_ScheduleLabel, &
      OMP_ScheduleLabelPrefix
    character ( LDF ) :: &
      Filename
    logical ( KDL ) :: &
      AppendDimensionality, &
      DimensionalityFound
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    procedure ( ), pointer :: &
      Abort
      
!-- Runtime error with CCE
!    if ( KBCH == selected_char_kind ( 'ASCII' ) ) then
!      open ( OUTPUT_UNIT, encoding = 'DEFAULT' )
!    else if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
    if ( KBCH == selected_char_kind ( 'ISO_10646' ) ) then
      Encoding = 'UTF-8'
      open ( OUTPUT_UNIT, encoding = Encoding )
    end if
    
    AppendDimensionality = .true.
    if ( present ( AppendDimensionalityOption ) ) &
      AppendDimensionality = AppendDimensionalityOption
      
    PH => PROGRAM_HEADER 
      
    allocate ( PH % Communicator )
    call PH % Communicator % Initialize ( )
    
    call UNIT % Initialize ( )
    
    Abort => Abort_PH
    call CONSOLE % Initialize ( PH % Communicator % Rank, AbortOption = Abort )

    allocate ( PH % ParameterStream )
    allocate ( PH % CommandLineOptions )
    associate &
      ( PS   => PH % ParameterStream, &
        CLO => PH % CommandLineOptions )

    call CLO % Initialize ( )

    Verbosity = CONSOLE % LABEL ( CONSOLE % Verbosity ) 
    call CLO % Read ( Verbosity, 'Verbosity', CONSOLE % INFO_1 )
    call CONSOLE % SetVerbosity ( Verbosity )

    DisplayRank  = CONSOLE % DisplayRank
    call CLO % Read ( DisplayRank, 'DisplayRank', CONSOLE % INFO_1 )
    call PH % Communicator % Synchronize ( )
    call CONSOLE % SetDisplayRank ( DisplayRank )
    
    call PrepareAndShow_OMP_Environment ( )
    
    if ( AppendDimensionality ) then
      if ( present ( DimensionalityOption ) ) &
        PH % Dimensionality = DimensionalityOption
      call CLO % Read &
             ( PH % Dimensionality, 'Dimensionality', &
               IgnorabilityOption = CONSOLE % INFO_1, &
               SuccessOption = DimensionalityFound )
      if ( .not. present ( DimensionalityOption ) &
           .and. .not. DimensionalityFound ) &
      then
        PH % Dimensionality = '3D'
        call Show ( 'Dimensionality not specified, defaulting to 3D', &
                    CONSOLE % INFO_1 )
        call Show ( PH % Dimensionality, 'Dimensionality', CONSOLE % INFO_1 )
      end if
      PH % Name = trim ( Name ) // '_' // trim ( PH % Dimensionality )
    else
      PH % Name = Name
    end if

    Filename = trim ( PH % Name ) // '_Program_Parameters'
    call PH % ParameterStream % Initialize &
           ( Filename, PH % Communicator % Rank, &
             IgnorabilityOption = CONSOLE % INFO_1 )

    call Show ( 'Setting Timer parameters', CONSOLE % INFO_1 )
    allocate ( PH % Timer ( MAX_TIMERS ) )

    PH % TimerLevel = 8
    call PROGRAM_HEADER % GetParameter ( PH % TimerLevel, 'TimerLevel' )
    call Show ( PH % TimerLevel, 'TimerLevel', CONSOLE % INFO_1 )

    PH % TimerDisplayFraction = 0.1_KDR
    call PROGRAM_HEADER % GetParameter &
           ( PH % TimerDisplayFraction, 'TimerDisplayFraction' )
    call Show ( PH % TimerDisplayFraction, 'TimerDisplayFraction', &
                CONSOLE % INFO_1 )

    call PH % AddTimer &
           ( 'Execution', Level = 0, Handle = PH % ExecutionTimeHandle )
    call PH % Timer ( PH % ExecutionTimeHandle ) % Start ( )

!    call Show ( 'Initializing PETSc', CONSOLE % INFO_1)
!    call PETSCINITIALIZE ( PETSC_NULL_CHARACTER, Error )

    call Show ( 'Starting the Program', CONSOLE % INFO_1 ) 
    call Show ( PH % Name, 'Name', CONSOLE % INFO_1 ) 

    end associate  !-- P, CLO 

    nullify ( Abort )
    nullify ( PH )

  end subroutine Initialize
  
  
  subroutine GetParameter_0D_Integer &
               ( Value, Name, ParameterStreamOption, IgnorabilityOption, &
                 SuccessOption )

    integer ( KDI ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_0D_Integer


  subroutine GetParameter_0D_Real &
               ( Value, Name, InputUnitOption, ParameterStreamOption, &
                 IgnorabilityOption, SuccessOption )

    real ( KDR ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_0D_Real


  subroutine GetParameter_0D_MeasuredValue &
               ( Value, Name, InputUnitOption, ParameterStreamOption, &
                 IgnorabilityOption, ConvertOption, SuccessOption )

    type ( MeasuredValueForm ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( inout ), optional :: &
      InputUnitOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( in ), optional :: &
      ConvertOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             ConvertOption = ConvertOption, SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             ConvertOption = ConvertOption, SuccessOption = Success_CLO )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_0D_MeasuredValue


  subroutine GetParameter_0D_Logical &
               ( Value, Name, ParameterStreamOption, IgnorabilityOption, &
                 SuccessOption )

    logical ( KDL ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_0D_Logical


  subroutine GetParameter_0D_Character &
               ( Value, Name, ParameterStreamOption, IgnorabilityOption, &
                 SuccessOption )

    character ( * ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_0D_Character


  subroutine GetParameter_1D_Integer &
               ( Value, Name, nValuesOption, ParameterStreamOption, &
                 IgnorabilityOption, SuccessOption )

    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_Integer


  subroutine GetParameter_1D_Real &
               ( Value, Name, InputUnitOption, nValuesOption, &
                 ParameterStreamOption, IgnorabilityOption, SuccessOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_Real


  subroutine GetParameter_1D_MeasuredValue &
               ( Value, Name, InputUnitOption, nValuesOption, &
                 ParameterStreamOption, IgnorabilityOption, SuccessOption )

    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( inout ), optional :: &
      InputUnitOption
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_MeasuredValue


  subroutine GetParameter_1D_Logical &
               ( Value, Name, nValuesOption, ParameterStreamOption, &
                 IgnorabilityOption, SuccessOption )

    logical ( KDL ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_Logical


  subroutine GetParameter_1D_Character &
               ( Value, Name, nValuesOption, ParameterStreamOption, &
                 IgnorabilityOption, SuccessOption )

    character ( * ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ), optional :: &
      nValuesOption
    type ( ParameterStreamForm ), intent ( in ), target, optional :: &
      ParameterStreamOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption

    integer ( KDI ) :: &
      Ignorability
    logical ( KDL ) :: &
      Success_PS, &
      Success_CLO
    type ( ParameterStreamForm ), pointer :: &
      PS

    Ignorability = CONSOLE % INFO_4
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show ( 'Parameter ' // trim ( Name ) // ' default value', &
                Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = Ignorability, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_Character


  subroutine AddTimer ( Name, Handle, Level )

    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( inout ) :: &
      Handle
    integer ( KDI ), intent ( in ) :: &
      Level    

    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    
    PH => PROGRAM_HEADER 
      
    if ( Level > PH % TimerLevel ) &
      return

    PH % nTimers = PH % nTimers + 1
    Handle = PH % nTimers

    call PH % Timer ( Handle ) % Initialize ( Name, Level )

    call Show ( 'Adding a Timer', CONSOLE % INFO_2 )
    call Show ( PH % Timer ( Handle ) % Name, 'Name', CONSOLE % INFO_2 )

  end subroutine AddTimer


  function TimerPointer ( Handle ) result ( TP )

    integer ( KDI ), intent ( in ) :: &
      Handle
    type ( TimerForm ), pointer :: &
      TP

    TP => null ( )

    if ( Handle > 0 ) &
      TP => PROGRAM_HEADER % Timer ( Handle )

  end function TimerPointer


  subroutine ShowStatistics &
               ( Ignorability, CommunicatorOption, MaxTimeOption, &
                 MinTimeOption, MeanTimeOption )
  
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( out ), optional :: &
      MaxTimeOption, &
      MinTimeOption, &
      MeanTimeOption
      
    type ( MeasuredValueForm )  :: &
      HighWaterMark, &
      MaxHighWaterMark, &
      MinHighWaterMark, &
      MeanHighWaterMark, &
      ResidentSetSize, &
      MaxResidentSetSize, &
      MinResidentSetSize, &
      MeanResidentSetSize
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
      
    PH => PROGRAM_HEADER 
      
    call Show ( 'Program timing', Ignorability )

    call ReadTimers &
           ( Ignorability, CommunicatorOption, MaxTimeOption, MinTimeOption, &
             MeanTimeOption )

    call Show ( 'Program memory usage', Ignorability )
    
    if ( present ( CommunicatorOption ) ) then
      call GetMemoryUsage &
             ( HighWaterMark, ResidentSetSize, Ignorability, &
               C_Option        = CommunicatorOption, &
               Max_HWM_Option  = MaxHighWaterMark, &
               Min_HWM_Option  = MinHighWaterMark, &
               Mean_HWM_Option = MeanHighWaterMark, &
               Max_RSS_Option  = MaxResidentSetSize, &
               Min_RSS_Option  = MinResidentSetSize, &
               Mean_RSS_Option = MeanResidentSetSize )
    else
      call GetMemoryUsage &
             ( HighWaterMark, ResidentSetSize, Ignorability )
    end if
    
    call Show ( HighWaterMark, 'This process HWM', Ignorability + 1 )
    call Show ( ResidentSetSize, 'This process RSS', Ignorability + 1 )
    
    if ( present ( CommunicatorOption ) ) then
      
      call Show ( MaxHighWaterMark, 'Across processes max HWM', &
                  Ignorability )
      call Show ( MinHighWaterMark, 'Across processes min HWM', &
                  Ignorability + 1 )
      call Show ( MeanHighWaterMark, 'Across processes mean HWM', &
                  Ignorability + 1 )
      
      call Show ( MaxResidentSetSize, 'Across processes max RSS', &
                  Ignorability )
      call Show ( MinResidentSetSize, 'Across processes min RSS', &
                  Ignorability + 1 )
      call Show ( MeanResidentSetSize, 'Across processes mean RSS', &
                  Ignorability + 1 )
    
    end if
    
    nullify ( PH )

  end subroutine ShowStatistics 
  
  
  subroutine Abort_PH ( )
  
    if ( PROGRAM_HEADER % Communicator % Initialized ) then
      call PROGRAM_HEADER % Communicator % Abort ( )
    else
      stop 1
    end if
  
  end subroutine Abort_PH


  impure elemental subroutine Finalize ( PH ) 

    type ( ProgramHeaderSingleton ), intent ( inout ) :: &
      PH
    
    call Show ( 'Finishing the Program', CONSOLE % INFO_1 ) 
    call Show ( PH % Name, 'Name', CONSOLE % INFO_1 )
    
!    call Show ( 'Finalizing PETSc', CONSOLE % INFO_1)
!    call PETSCFINALIZE ( Error )
    
    call PH % ShowStatistics &
           ( CONSOLE % INFO_1, &
             CommunicatorOption = PROGRAM_HEADER % Communicator )

    if ( allocated ( PH % CommandLineOptions ) ) &
      deallocate ( PH % CommandLineOptions ) 
    if ( allocated ( PH % Timer ) ) &
      deallocate ( PH % Timer )
    if ( allocated ( PH % ParameterStream ) ) &
      deallocate ( PH % ParameterStream )
    if ( allocated ( PH % Communicator ) ) &
      deallocate ( PH % Communicator )

  end subroutine Finalize
  
  
  subroutine PrepareAndShow_OMP_Environment ( )
  
    integer ( KDI )  :: &
      Length, &
      Status, &
      OMP_ScheduleChunkSize
    integer ( OMP_SCHED_KIND ) :: &
      OMP_ScheduleKind
    character ( LDL )  :: &
      OMP_ScheduleLabel, &
      OMP_ScheduleLabelPrefix
    character ( LDB ) :: &
      OMP_SetSchedule 
    
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    
    PH => PROGRAM_HEADER 
  
    PH % MaxThreads = OMP_GET_MAX_THREADS ( )
    
    call get_environment_variable &
           ( 'OMP_SCHEDULE', OMP_SetSchedule, Length, Status, .false. )
    
    !-- Set default OMP schedule if not in env. var.
    if ( Length == 0 ) &
      call OMP_SET_SCHEDULE ( OMP_SCHED_GUIDED, -1 )
    
    OMP_ScheduleLabelPrefix = ''
    if ( OffloadEnabled ( ) ) then  !-- hardcoded to (static, 1) with offload
      OMP_ScheduleKind = OMP_SCHED_STATIC
      OMP_ScheduleChunkSize = 1
    else
      call OMP_GET_SCHEDULE ( OMP_ScheduleKind, OMP_ScheduleChunkSize )
      OMP_ScheduleLabelPrefix = 'runtime : '
    end if
    
    select case ( OMP_ScheduleKind )
    case ( OMP_SCHED_STATIC )
      OMP_ScheduleLabel = 'static'
    case ( OMP_SCHED_DYNAMIC )
      OMP_ScheduleLabel = 'dynamic'
    case ( OMP_SCHED_GUIDED )
      OMP_ScheduleLabel = 'guided'
    case ( OMP_SCHED_AUTO )
      OMP_ScheduleLabel = 'auto'
    end select    
    
    call Show ( 'OpenMP environment', CONSOLE % INFO_1 )
    call Show ( PH % MaxThreads,  'MaxThreads', CONSOLE % INFO_1 )
    call Show ( GetNumberOfDevices ( ), 'nDevices', CONSOLE % INFO_1 )
    call Show ( OffloadEnabled ( ), 'Offload enabled', CONSOLE % INFO_1 )
    call Show &
           ( adjustl ( adjustr ( OMP_ScheduleLabelPrefix ) &
                    // ' ' // adjustl ( OMP_ScheduleLabel ) ), &
             'Schedule', CONSOLE % INFO_1 )
    call Show ( OMP_ScheduleChunkSize, 'ChunkSize', CONSOLE % INFO_1 )
  
  end subroutine PrepareAndShow_OMP_Environment
  

  subroutine ReadTimers &
               ( Ignorability, CommunicatorOption, MaxTimeOption, &
                 MinTimeOption, MeanTimeOption )

    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( out ), optional :: &
      MaxTimeOption, &
      MinTimeOption, &
      MeanTimeOption
      
    integer ( KDI ) :: &
      iT
    real ( KDR ) :: &
      ExecutionTime
    logical ( KDL ), dimension ( MAX_TIMERS ) :: &
      Running
    type ( CollectiveOperation_R_Form ) :: &
      CO
    type ( TimerForm ), dimension ( MAX_TIMERS ) :: &
      MaxTimer, &
      MinTimer, &
      MeanTimer
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
   
    PH => PROGRAM_HEADER 

    call Show ( 'Running timer intervals', Ignorability + 2 )
    Running = .false.
    do iT = 1, PH % nTimers
      if ( PH % Timer ( iT ) % Running ) then
        Running ( iT ) = .true.
        call PH % Timer ( iT ) % Stop ( )
        call PH % Timer ( iT ) % ShowInterval ( Ignorability + 2 )
      end if
    end do

    call Show ( 'This process timers', Ignorability + 1 )
    do iT = 1, PH % nTimers
      call PH % Timer ( iT ) % ShowTotal ( Ignorability + 1 )
    end do !-- iT

    if ( present ( CommunicatorOption ) ) then

      call CO % Initialize &
             ( CommunicatorOption, &
               nOutgoing = [ PH % nTimers ], nIncoming = [ PH % nTimers ] )

      do iT = 1, PH % nTimers
        call MaxTimer ( iT ) % Initialize ( PH % Timer ( iT ) )
        call MinTimer ( iT ) % Initialize ( PH % Timer ( iT ) )
        call MeanTimer ( iT ) % Initialize ( PH % Timer ( iT ) )
        CO % Outgoing % Value ( iT ) = PH % Timer ( iT ) % TotalTime
      end do !-- iT

      call Show ( 'Max timers', Ignorability + 1 )
      call CO % Reduce ( REDUCTION % MAX )
      do iT = 1, PH % nTimers
        call MaxTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call MaxTimer ( iT ) % ShowTotal ( Ignorability + 1 )
        if ( present ( MaxTimeOption ) ) &
          MaxTimeOption ( iT ) = MaxTimer ( iT ) % TotalTime
      end do !-- iT
      
      call Show ( 'Min timers', Ignorability + 1 )
      call CO % Reduce ( REDUCTION % MIN )
      do iT = 1, PH % nTimers
        call MinTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call MinTimer ( iT ) % ShowTotal ( Ignorability + 1 )
        if ( present ( MinTimeOption ) ) &
          MinTimeOption ( iT ) = MinTimer ( iT ) % TotalTime
      end do !-- iT
      
      call Show ( 'Mean timers', Ignorability )
      call CO % Reduce ( REDUCTION % SUM )
      do iT = 1, PH % nTimers
        call MeanTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) &
                      / CommunicatorOption % Size )
        if ( iT == 1 ) &
          ExecutionTime = MeanTimer ( iT ) % TotalTime
        if ( MeanTimer ( iT ) % TotalTime / ExecutionTime &
             >= PH % TimerDisplayFraction ) &
        then
          call MeanTimer ( iT ) % ShowTotal ( Ignorability )
        else
          call MeanTimer ( iT ) % ShowTotal ( Ignorability + 1 )
        end if
        if ( present ( MeanTimeOption ) ) &
          MeanTimeOption ( iT ) = MeanTimer ( iT ) % TotalTime
      end do !-- iT
      
    end if !-- present ( CommunicatorOption )

    do iT = 1, PH % nTimers
      if ( Running ( iT ) ) call PH % Timer ( iT ) % Start ( )
    end do

  end subroutine ReadTimers


end module PROGRAM_HEADER_Singleton
