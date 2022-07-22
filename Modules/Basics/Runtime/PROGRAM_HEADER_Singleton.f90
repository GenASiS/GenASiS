!-- PROGRAM_HEADER provides functionalities commonly required by all programs
!   (drivers), including initialization of parallel environment, obtaining
!   program parameters, and displaying basic runtime statistics.

module PROGRAM_HEADER_Singleton

  use iso_fortran_env
  use OMP_LIB
  use Specifiers
  use Devices
  use Display
  use MessagePassing
  use FileSystem
  use InitializeRandomSeed_Command
  use CommandLineOptions_Form
  use Timer_Form
  use Timer_1D__Form
  use MemoryUsage_Form
  !  use petsc

  implicit none
  private

  type, public :: ProgramHeaderSingleton
    integer ( KDI ) :: &
      MaxThreads      = 0, &
      iTimerExecution = 0
    character ( LDL ) :: &
      Dimensionality = ''
    character ( LDF ) :: &
      Name
    type ( CommunicatorForm ), allocatable :: &
      Communicator 
    type ( ParameterStreamForm ), allocatable :: &
      ParameterStream
    type ( CommandLineOptionsForm ), allocatable :: &
      CommandLineOptions
    type ( Timer_1D_Form ), allocatable :: &
      Timer_1D
    type ( MemoryUsageForm ), allocatable :: &
      MemoryUsage
  contains
    procedure, public, nopass :: &
      Initialize
    procedure, private, nopass :: &
      GetParameter_0D_Integer
    procedure, private, nopass :: &
      GetParameter_0D_Real
    procedure, private, nopass :: &
      GetParameter_0D_Quantity
    procedure, private, nopass :: &
      GetParameter_0D_Logical
    procedure, private, nopass :: &
      GetParameter_0D_Character
    procedure, private, nopass :: &
      GetParameter_1D_Integer
    procedure, private, nopass :: &
      GetParameter_1D_Real
    procedure, private, nopass :: &
      GetParameter_1D_Quantity
    procedure, private, nopass :: &
      GetParameter_1D_Logical
    procedure, private, nopass :: &
      GetParameter_1D_Character
    generic :: &
      GetParameter &
        => GetParameter_0D_Integer, GetParameter_0D_Real, &
           GetParameter_0D_Quantity, GetParameter_0D_Logical, &
           GetParameter_0D_Character, &
           GetParameter_1D_Integer, GetParameter_1D_Real, &
           GetParameter_1D_Quantity, GetParameter_1D_Logical, &
           GetParameter_1D_Character
    procedure, public, nopass :: &
      Timer
    procedure, public, nopass :: &
      RecordStatistics
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
      OMP_ScheduleChunkSize, &
      TimerLevelMin
    integer ( OMP_SCHED_KIND ) :: &
      OMP_ScheduleKind
    real ( KDR ) :: &
      TimerDisplayFraction
    logical ( KDL ) :: &
      AppendDimensionality, &
      DimensionalityFound
    character ( 5 ) :: &
      Encoding
    character ( LDL )  :: &
      Verbosity, &
      OMP_ScheduleLabel, &
      OMP_ScheduleLabelPrefix
    character ( LDF ) :: &
      Filename
    type ( TimerForm ), pointer :: &
      T
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
    
    !call InitializeRandomSeed ( PH % Communicator )

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

    TimerDisplayFraction  =  0.1_KDR
    TimerLevelMin  =  2
    call PH % GetParameter ( TimerDisplayFraction, 'TimerDisplayFraction' )
    call PH % GetParameter ( TimerLevelMin, 'TimerLevelMin' )

    allocate ( PH % Timer_1D )
    associate ( T_1D  =>  PH % Timer_1D )
    call T_1D % Initialize ( TimerDisplayFraction, TimerLevelMin )
    call T_1D % Get ( PH % iTimerExecution, 'Execution', 0, T )
    call T % Start ( )
    end associate !-- T_1D

    allocate ( PH % MemoryUsage )

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
    type ( QuantityForm ), intent ( inout ), optional :: &
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


  subroutine GetParameter_0D_Quantity &
               ( Value, Name, InputUnitOption, ParameterStreamOption, &
                 IgnorabilityOption, ConvertOption, SuccessOption )

    type ( QuantityForm ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), intent ( inout ), optional :: &
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

  end subroutine GetParameter_0D_Quantity


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
    type ( QuantityForm ), dimension ( : ), intent ( inout ), optional :: &
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


  subroutine GetParameter_1D_Quantity &
               ( Value, Name, InputUnitOption, nValuesOption, &
                 ParameterStreamOption, IgnorabilityOption, SuccessOption )

    type ( QuantityForm ), dimension ( : ), intent ( inout ) :: &
      Value
    character ( * ), intent ( in ) :: &
      Name
    type ( QuantityForm ), dimension ( : ), intent ( inout ), optional :: &
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

  end subroutine GetParameter_1D_Quantity


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


  function Timer ( Handle, Name, Level ) result ( T )

    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( in ) :: &
      Level    
    type ( TimerForm ), pointer :: &
      T

    type ( ProgramHeaderSingleton ), pointer :: &
      PH

    PH  =>  PROGRAM_HEADER 

    associate ( T_1D  =>  PH % Timer_1D )

    if ( Handle  <=  0 ) then  !-- New

      if ( T_1D % nTimers  ==  T_1D % MAX_TIMERS ) then
        call Show ( 'Maximum number of timers reached', CONSOLE % ERROR )
        call Show ( T_1D % nTimers, 'nTimers', CONSOLE % ERROR )
        call Show ( T_1D % MAX_TIMERS, 'MAX_TIMERS', CONSOLE % ERROR )
        call Show ( Name, 'Requested timer', CONSOLE % ERROR )
        call Show ( PH % Communicator % Rank, 'Rank', CONSOLE % ERROR )
        call PH % Abort ( )
      end if

      if ( Level  >  T_1D % LevelMax ) then
        call Show ( 'Timer level exceeds maximum', CONSOLE % ERROR )
        call Show ( Level, 'Level', CONSOLE % ERROR )
        call Show ( T_1D % LevelMax, 'LevelMax', CONSOLE % ERROR )
        call Show ( Name, 'Requested timer', CONSOLE % ERROR )
        call Show ( PH % Communicator % Rank, 'Rank', CONSOLE % ERROR )
        call PH % Abort ( )
      end if

    endif  !-- New

    call T_1D % Get ( Handle, Name, Level, T )

    end associate !-- T_1D

  end function Timer


  subroutine RecordStatistics ( Ignorability, CommunicatorOption )
  
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
      
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    
    PH => PROGRAM_HEADER 
  
    call Show ( 'Program timing', Ignorability )
    call ReadTimers ( Ignorability, CommunicatorOption )

    call Show ( 'Program memory usage', Ignorability )
    call PH % MemoryUsage % Compute ( Ignorability, CommunicatorOption )

  end subroutine RecordStatistics 
  
  
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
    
    call PH % RecordStatistics &
           ( CONSOLE % INFO_1, &
             CommunicatorOption = PROGRAM_HEADER % Communicator )

    if ( allocated ( PH % MemoryUsage ) ) &
      deallocate ( PH % MemoryUsage )
    if ( allocated ( PH % Timer_1D ) ) &
      deallocate ( PH % Timer_1D )
    if ( allocated ( PH % CommandLineOptions ) ) &
      deallocate ( PH % CommandLineOptions ) 
    if ( allocated ( PH % ParameterStream ) ) &
      deallocate ( PH % ParameterStream )
    if ( allocated ( PH % Communicator ) ) &
      deallocate ( PH % Communicator )

  end subroutine Finalize
  
  
  subroutine PrepareAndShow_OMP_Environment ( )
  
    integer ( KDI )  :: &
      nDevices, &
      iDefaultDevice, &
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
    
    !-- Set default OMP schedule to "guided" if not set in env. var.
    if ( Length == 0 ) &
      call OMP_SET_SCHEDULE ( OMP_SCHED_GUIDED, -1 )
    
    OMP_ScheduleLabelPrefix = ''
    if ( OffloadEnabled ( ) .and. NumberOfDevices ( ) >= 1 ) then  
      !-- per Build/Preprocessor, hardcoded to "(static, 1)" for offload
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
    
    !-- Set default device for offload based on MPI rank 
    nDevices = NumberOfDevices ( )
    if ( nDevices > 0 ) then
      iDefaultDevice = mod ( PH % Communicator % Rank, nDevices )
      call SelectDevice ( iDefaultDevice )
    end if
    
    call Show ( 'OpenMP environment', CONSOLE % INFO_1 )
    call Show ( PH % MaxThreads,  'MaxThreads', CONSOLE % INFO_1 )
    call Show ( nDevices, 'nDevices', CONSOLE % INFO_1 )
    call Show ( OffloadEnabled ( ), 'OffloadEnabled', CONSOLE % INFO_1 )
    call Show ( SelectedDevice ( ), 'Selected device', &
                CONSOLE % INFO_1 )
    call Show &
           ( adjustl ( adjustr ( OMP_ScheduleLabelPrefix ) &
                    // ' ' // adjustl ( OMP_ScheduleLabel ) ), &
             'Schedule', CONSOLE % INFO_1 )
    call Show ( OMP_ScheduleChunkSize, 'ChunkSize', CONSOLE % INFO_1 )
  
  end subroutine PrepareAndShow_OMP_Environment
  

  subroutine ReadTimers ( Ignorability, CommunicatorOption )

    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption

    type ( CollectiveOperation_I_Form ) :: &
      CO
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
   
    PH => PROGRAM_HEADER 

    if ( present ( CommunicatorOption ) ) then

      !-- nTimers sanity check
      call CO % Initialize &
             ( CommunicatorOption, nOutgoing = [ 1 ], &
               nIncoming = [ CommunicatorOption % Size ], &
               RootOption = CONSOLE % DisplayRank )
      CO % Outgoing % Value  =  PH % Timer_1D % nTimers
      call CO % Gather ( )
      if ( CommunicatorOption % Rank  ==  CONSOLE % DisplayRank ) then
        if ( any ( CO % Incoming % Value  /=  CO % Incoming % Value ( 1 ) ) ) &
        then
          call Show ( 'Unequal number of timers across processes', &
                      CONSOLE % ERROR )
          call Show ( 'PROGRAM_HEADER_Singleton', 'module', CONSOLE % ERROR )
          call Show ( 'ReadTimers', 'subroutine', CONSOLE % ERROR )
          call Show ( CO % Incoming % Value, 'nTimers', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end if !-- unequal nTimers
      end if !-- CONSOLE % DisplayRank

    end if !-- CommunicatorOption

    call PH % Timer_1D % Read ( Ignorability, CommunicatorOption )

  end subroutine ReadTimers


end module PROGRAM_HEADER_Singleton
