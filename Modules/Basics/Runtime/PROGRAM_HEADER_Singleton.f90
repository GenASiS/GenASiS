!-- PROGRAM_HEADER provides functionalities commonly required by all programs
!   (drivers), including initialization of parallel environment, obtaining
!   program parameters, and displaying basic runtime statistics.

module PROGRAM_HEADER_Singleton

  use OMP_LIB
  use VariableManagement
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
      MAX_TIMERS = 32

  type, public :: ProgramHeaderSingleton
    integer ( KDI ) :: &
      nTimers = 0, &
      MaxThreads = 0, &
      ExecutionTimeHandle
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
      ShowStatistics
    procedure, public, nopass :: &
      Abort => Abort_PH  !-- avoids conflict with intrinsic "abort"
    final :: &
      Finalize
  end type ProgramHeaderSingleton
  
  type ( ProgramHeaderSingleton ), public, target, allocatable :: &
    PROGRAM_HEADER

    private :: &
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
      DisplayRank
    character ( LDL )  :: &
      Verbosity
    character ( LDF ) :: &
      Filename
    logical ( KDL ) :: &
      AppendDimensionality, &
      DimensionalityFound
    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    procedure ( ), pointer :: &
      Abort
      
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
        call Show &
               ( 'Dimensionality not specified, defaulting to 3D', &
                 CONSOLE % WARNING )
        call Show ( PH % Dimensionality, 'Dimensionality' )
      end if
      PH % Name = trim ( Name ) // '_' // trim ( PH % Dimensionality )
    else
      PH % Name = Name
    end if

    Filename = trim ( PH % Name ) // '_Program_Parameters'
    call PH % ParameterStream % Initialize &
           ( Filename, PH % Communicator % Rank )

    DisplayRank  = CONSOLE % DisplayRank
    call PH % GetParameter ( DisplayRank, 'DisplayRank' )
    call PH % Communicator % Synchronize ( )
    call CONSOLE % SetDisplayRank ( DisplayRank )

    Verbosity = CONSOLE % LABEL ( CONSOLE % Verbosity ) 
    call PH % GetParameter ( Verbosity, 'Verbosity' )
    call CONSOLE % SetVerbosity ( Verbosity )

    PH % MaxThreads = OMP_GET_MAX_THREADS ( )
    call Show ( PH % MaxThreads, 'MaxThreads', CONSOLE % INFO_1 )
    
    allocate ( PH % Timer ( MAX_TIMERS ) )
    call PH % AddTimer ( 'Execution', PH % ExecutionTimeHandle )
    call PH % Timer ( PH % ExecutionTimeHandle ) % Start ( )
    
!    call Show ( 'Initializing PETSc', CONSOLE % INFO_1)
!    call PETSCINITIALIZE ( PETSC_NULL_CHARACTER, Error )

    call Show ( 'Starting the Program', CONSOLE % INFO_1 ) 
    call Show ( PH % Name, 'Name', CONSOLE % INFO_1 ) 
    call Show ( PH % MaxThreads, 'MaxThreads', CONSOLE % INFO_1 )

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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
             ConvertOption = ConvertOption, SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, InputUnitOption = InputUnitOption, &
             IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
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

    Ignorability = CONSOLE % INFO_3
    if ( present ( IgnorabilityOption ) ) Ignorability = IgnorabilityOption

    if ( present ( ParameterStreamOption ) ) then
      PS => ParameterStreamOption
    else
      PS => PROGRAM_HEADER % ParameterStream
    end if

    associate ( CLO => PROGRAM_HEADER % CommandLineOptions )

    call Show &
           ( 'Parameter ' // trim ( Name ) // ' default value', &
             Ignorability )
    call Show ( Value, Name, Ignorability )
    call PS % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_PS, nValuesOption = nValuesOption )
    call CLO % Read &
           ( Value, Name, IgnorabilityOption = IgnorabilityOption, &
             SuccessOption = Success_CLO, nValuesOption = nValuesOption )
    if ( present ( SuccessOption ) ) &
      SuccessOption = Success_PS .or. Success_CLO

    end associate !-- CLO

    nullify ( PS )

  end subroutine GetParameter_1D_Character


  subroutine AddTimer ( Name, Handle )

    character ( * ), intent ( in ) :: &
      Name
    integer ( KDI ), intent ( out ) :: &
      Handle

    type ( ProgramHeaderSingleton ), pointer :: &
      PH
    
    PH => PROGRAM_HEADER 
      
    PH % nTimers = PH % nTimers + 1
    Handle = PH % nTimers

    call PH % Timer ( Handle ) % Initialize ( Name )

    call Show ( 'Adding a Timer', CONSOLE % INFO_1 )
    call Show ( PH % Timer ( Handle ) % Name, 'Name' )

  end subroutine AddTimer


  subroutine ShowStatistics &
               ( Verbosity, CommunicatorOption, MaxTimeOption, MinTimeOption, &
                 MeanTimeOption )
  
    integer ( KDI ), intent ( in ) :: &
      Verbosity
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
      
    if ( Verbosity > CONSOLE % Verbosity ) return

    PH => PROGRAM_HEADER 
      
    call Show ( 'Runtime statistics', Verbosity )

    call ReadTimers &
           ( Verbosity, CommunicatorOption, MaxTimeOption, MinTimeOption, &
             MeanTimeOption )

    call Show ( 'Program memory usage', Verbosity )
    
    if ( present ( CommunicatorOption ) ) then
      call GetMemoryUsage &
             ( HighWaterMark, ResidentSetSize, CommunicatorOption, &
               Max_HWM_Option  = MaxHighWaterMark, &
               Min_HWM_Option  = MinHighWaterMark, &
               Mean_HWM_Option = MeanHighWaterMark, &
               Max_RSS_Option  = MaxResidentSetSize, &
               Min_RSS_Option  = MinResidentSetSize, &
               Mean_RSS_Option = MeanResidentSetSize )
    else
      call GetMemoryUsage &
             ( HighWaterMark, ResidentSetSize )
    end if
    
    call Show ( HighWaterMark, 'This process HWM', Verbosity )
    call Show ( ResidentSetSize, 'This process RSS', Verbosity )
    
    if ( present ( CommunicatorOption ) ) then
      
      call Show ( MaxHighWaterMark, 'Across processes max HWM', &
                  Verbosity )
      call Show ( MinHighWaterMark, 'Across processes min HWM', &
                  Verbosity )
      call Show ( MeanHighWaterMark, 'Across processes mean HWM', &
                  Verbosity )
      
      call Show ( MaxResidentSetSize, 'Across processes max RSS', &
                  Verbosity )
      call Show ( MinResidentSetSize, 'Across processes min RSS', &
                  Verbosity )
      call Show ( MeanResidentSetSize, 'Across processes mean RSS', &
                  Verbosity )
    
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
  

  subroutine ReadTimers &
               ( Verbosity, CommunicatorOption, MaxTimeOption, MinTimeOption, &
                 MeanTimeOption )

    integer ( KDI ), intent ( in ) :: &
      Verbosity
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    real ( KDR ), dimension ( : ), intent ( out ), optional :: &
      MaxTimeOption, &
      MinTimeOption, &
      MeanTimeOption
      
    integer ( KDI ) :: &
      iT
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

    call Show ( 'Running timer intervals', Verbosity )
    Running = .false.
    do iT = 1, PH % nTimers
      if ( PH % Timer ( iT ) % Running ) then
        Running ( iT ) = .true.
        call PH % Timer ( iT ) % Stop ( )
        call PH % Timer ( iT ) % ShowInterval ( Verbosity )
      end if
    end do

    call Show ( 'This process timers', Verbosity )
    do iT = 1, PH % nTimers
      call PH % Timer ( iT ) % ShowTotal ( Verbosity )
    end do !-- iT

    if ( present ( CommunicatorOption ) ) then

      call CO % Initialize &
             ( CommunicatorOption, &
               nOutgoing = [ PH % nTimers ], nIncoming = [ PH % nTimers ], &
               RootOption = CONSOLE % DisplayRank )

      do iT = 1, PH % nTimers
        call MaxTimer ( iT ) % Initialize ( PH % Timer ( iT ) % Name )
        call MinTimer ( iT ) % Initialize ( PH % Timer ( iT ) % Name )
        call MeanTimer ( iT ) % Initialize ( PH % Timer ( iT ) % Name )
        CO % Outgoing % Value ( iT ) = PH % Timer ( iT ) % TotalTime
      end do !-- iT

      call Show ( 'Max timers', Verbosity )
      call CO % Reduce ( REDUCTION % MAX )
      do iT = 1, PH % nTimers
        call MaxTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call MaxTimer ( iT ) % ShowTotal ( Verbosity )
        if ( present ( MaxTimeOption ) ) &
          MaxTimeOption ( iT ) = MaxTimer ( iT ) % TotalTime
      end do !-- iT
      
      call Show ( 'Min timers', Verbosity )
      call CO % Reduce ( REDUCTION % MIN )
      do iT = 1, PH % nTimers
        call MinTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) )
        call MinTimer ( iT ) % ShowTotal ( Verbosity )
        if ( present ( MinTimeOption ) ) &
          MinTimeOption ( iT ) = MinTimer ( iT ) % TotalTime
      end do !-- iT
      
      call Show ( 'Mean timers', Verbosity )
      call CO % Reduce ( REDUCTION % SUM )
      do iT = 1, PH % nTimers
        call MeanTimer ( iT ) % TotalTime % Initialize &
               ( 's', CO % Incoming % Value ( iT ) / CommunicatorOption % Size )
        call MeanTimer ( iT ) % ShowTotal ( Verbosity )
        if ( present ( MeanTimeOption ) ) &
          MeanTimeOption ( iT ) = MeanTimer ( iT ) % TotalTime
      end do !-- iT
      
    end if !-- present ( CommunicatorOption )

    do iT = 1, PH % nTimers
      if ( Running ( iT ) ) call PH % Timer ( iT ) % Start ( )
    end do

  end subroutine ReadTimers


end module PROGRAM_HEADER_Singleton
