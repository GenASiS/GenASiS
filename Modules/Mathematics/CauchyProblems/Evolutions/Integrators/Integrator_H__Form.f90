module Integrator_H__Form

  !-- Integrator_Header_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

  type, public :: Integrator_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimer_E = 0, &   !-- Evolution
      iTimer_AC = 0, &  !-- AdministerCheckpoint
      ! iTimerCycle = 0, &
      ! iTimerNewTime = 0, &
      ! iTimerTally = 0, &
      ! iTimerAnalyze = 0, &
      ! iTimerWriteSeries = 0, &
      iCycle, &
      iCheckpoint, &
      nRampCycles, &
      FinishCycle, &
      nWrite, &
      n_dT_Candidates, &
      CheckpointDisplayInterval
    real ( KDR ) :: &
      T_Start, &
      T_Finish, &
      T_CheckpointInterval, &
      T_Checkpoint, &
      T
    type ( MeasuredValueForm ) :: &
      Unit_T
    real ( KDR ), dimension ( : ), allocatable :: &
      dT_Candidate
    logical ( KDL ) :: &
      Start, &
      Restart, &
      NoWrite, &
      AllWrite, &
      T_CheckpointExact, &
      CheckpointDue
    character ( LDL ), dimension ( : ), allocatable :: &
      dT_Label
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream
    class ( Atlas_H_Form ), allocatable :: &
      X_A
    type ( Stream_A_Form ), allocatable :: &
      Checkpoint_X_A
    class ( Geometry_F_A_Form ), allocatable :: &
      Geometry_X_A
    procedure ( SI ), public, pointer :: &
      SetInitial => null ( )
    procedure ( RI ), public, pointer :: &
      ResetInitial => null ( )
    procedure ( W ), public, pointer :: &
      Write => null ( )
    procedure ( R ), public, pointer :: &
      Read => null ( )
    procedure ( SCTI ), pointer :: &
      Set_T_CheckpointInterval => null ( )
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H
    generic, public :: &           
      Initialize => Initialize_H
    procedure, private, pass :: &  !-- 1
      Show_I
    generic, public :: &
      Show => Show_I
    procedure, public, pass :: &   !-- 1
      Evolve
    final :: &                     !-- 1
      Finalize
    procedure, private, pass :: &  !-- 2
      ShowManifold
    procedure, private, pass :: &  !-- 2
      ShowFields
    procedure, private, pass :: &  !-- 2
      ShowCheckpoint
    procedure, public, pass :: &   !-- 2
      PrepareInitial
    procedure, public, pass :: &   !-- 2
      PrepareEvolution
    procedure, private, pass :: &  !-- 2
      AdministerCheckpoint
    procedure, public, pass :: &   !-- 3
      SetInitial_H
    procedure, public, pass :: &   !-- 3
      ResetInitial_H
    procedure, public, pass :: &  !-- 3
      UpdateHost => UpdateHost_H
    procedure, public, pass :: &   !-- 3
      Write_H
    procedure, public, pass :: &   !-- 3
      Read_H
  end type Integrator_H_Form

  interface

    subroutine SI ( I )
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
    end subroutine SI

    subroutine RI ( I, RestartFrom, T_Restart )
      use Basics
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ) :: &
        RestartFrom
      type ( MeasuredValueForm ), intent ( out ) :: &
        T_Restart
    end subroutine RI

    subroutine W ( I, TimerLevelOption )
      use Basics
      import Integrator_H_Form
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ), optional :: &
        TimerLevelOption
    end subroutine W

    subroutine R ( I, ReadFrom, T, CycleNumber )
      use Basics
      import Integrator_H_Form
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ) :: &
        ReadFrom
      type ( MeasuredValueForm ), intent ( out ) :: &
        T
      integer ( KDI ), intent ( out ) :: &
        CycleNumber
    end subroutine R

    subroutine SCTI ( I )
      import Integrator_H_Form
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
    end subroutine SCTI

  end interface

  
    private :: &
      Set_T_CheckpointInterval
          

contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, &
                 Unit_T_Option, T_FinishOption, nWriteOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      Unit_T_Option
    real ( KDR ), intent ( in ), optional :: &
      T_FinishOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    character ( LDF ) :: &
      OutputDirectory

    I % IGNORABILITY  =  CONSOLE % INFO_1

    if ( I % Type == '' ) &
      I % Type = 'an Integrator' 

    I % Name = 'Integrator'
    if ( present ( NameOption ) ) &
      I % Name  =  NameOption

    !-- Communicator

    if ( present ( CommunicatorOption ) ) then
      I % Communicator  =>  CommunicatorOption
    else
      I % Communicator  =>  PROGRAM_HEADER % Communicator
    end if

    ! !-- Device

    ! I % DeviceMemory        =  OffloadEnabled ( )  &
    !                            .and.  GetNumberOfDevices ( ) >= 1 
    ! I % PinnedMemory        =  OffloadEnabled ( )  &
    !                            .and.  GetNumberOfDevices ( ) >= 1 
    ! I % DevicesCommunicate  =  OffloadEnabled ( )  &
    !                            .and.  GetNumberOfDevices ( ) >= 1
    ! call PROGRAM_HEADER % GetParameter &
    !        ( I % DeviceMemory, 'DeviceMemory' )
    ! call PROGRAM_HEADER % GetParameter &
    !        ( I % PinnedMemory, 'PinnedMemory' )
    ! call PROGRAM_HEADER % GetParameter &
    !        ( I % DevicesCommunicate, 'DevicesCommunicate' )

    !-- Atlas, if necessary

    if ( .not. allocated ( I % X_A ) ) then
      allocate ( Atlas_SCG_Form :: I % X_A )
      select type ( A  =>  I % X_A )
        class is ( Atlas_SCG_Form )
      call A % Initialize &
             ( CommunicatorOption = I % Communicator, &
               NameOption = 'X', &
               PeriodicOption = [ .true., .true., .true. ] )
      end select !-- A
    end if

    !-- Geometry, if necessary

    if ( .not. allocated ( I % Geometry_X_A ) ) then
      allocate ( I % Geometry_X_A )
      associate ( GA  =>  I % Geometry_X_A )
      call GA % Initialize &
             ( I % X_A, &
               DeviceMemoryOption = DeviceMemoryOption, &
               PinnedMemoryOption = PinnedMemoryOption, &
               DevicesCommunicateOption = DevicesCommunicateOption )
      end associate !-- GA
    end if

    !-- Integration parameters

    call Show ( 'Initializing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    I % T_Start   =  0.0_KDR
    I % T_Finish  =  1.0_KDR
    I % Unit_T    =  UNIT % IDENTITY
    if ( present ( T_FinishOption ) ) &
      I % T_Finish  =  T_FinishOption
    if ( present ( Unit_T_Option ) ) &
      I % Unit_T  =  Unit_T_Option
    call PROGRAM_HEADER % GetParameter &
           ( I % T_Start, 'T_Start', InputUnitOption = I % Unit_T )    
    call PROGRAM_HEADER % GetParameter &
           ( I % T_Finish, 'T_Finish', InputUnitOption = I % Unit_T )
           
    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'Candidate'
    end if
    I % n_dT_Candidates  =  size ( I % dT_Label )
    allocate ( I % dT_Candidate ( I % n_dT_Candidates ) )

    I % iCycle = 0
    I % iCheckpoint = 0
    I % nRampCycles = 100
    call PROGRAM_HEADER % GetParameter ( I % nRampCycles, 'nRampCycles' )

    I % FinishCycle = huge ( 1 )
    call PROGRAM_HEADER % GetParameter ( I % FinishCycle, 'FinishCycle' )

    !-- Checkpointing

    I % nWrite  =  100
    if ( present ( nWriteOption ) ) &
      I % nWrite  =  nWriteOption
    I %  NoWrite  =  .false.
    I % AllWrite  =  .false.
    call PROGRAM_HEADER % GetParameter ( I %   nWrite,   'nWrite' )
    call PROGRAM_HEADER % GetParameter ( I %  NoWrite,  'NoWrite' )
    call PROGRAM_HEADER % GetParameter ( I % AllWrite, 'AllWrite' )

    I % CheckpointDisplayInterval  =  100
    I % T_CheckpointExact  =  .false.
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval' )
    call PROGRAM_HEADER % GetParameter &
           ( I % T_CheckpointExact, 'T_CheckpointExact' )

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )

    allocate ( I % GridImageStream )
    associate ( GIS => I % GridImageStream )
    call GIS % Initialize &
           ( PROGRAM_HEADER % Name, &
             CommunicatorOption = I % Communicator, &
             WorkingDirectoryOption = OutputDirectory )
    end associate !-- GIS

    allocate ( I % Checkpoint_X_A )
    associate &
      (  SA  =>  I % Checkpoint_X_A, &
          A  =>  I % X_A, &
        GIS  =>  I % GridImageStream, &
         GA  =>  I % Geometry_X_A )
    call SA % Initialize ( A, GIS, NameOption = 'Checkpoint' )
    call GA % SetStream ( SA )
    end associate !--SA
    
  end subroutine Initialize_H


  subroutine Evolve ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

!     real ( KDR ) :: &
!       TimeStepRatio
    type ( TimerForm ), pointer :: &
      T

!     call I % OpenManifoldStreams ( )
!     call I % InitializeTimers ( )
!     call I % InitializeTimeSeries ( )

    associate ( iT  =>  I % iTimer_E )
    if ( iT == 0 ) then
      call PROGRAM_HEADER % AddTimer ( 'Evolve', iT, Level = 1 )
    end if
    end associate !-- iT
    T  =>  PROGRAM_HEADER % TimerPointer ( I % iTimer_E )

    call T % Start ( )

    call I % PrepareInitial ( )
    call I % PrepareEvolution ( )
    call I % AdministerCheckpoint ( ComputeChangeOption = .false. )

    call Show ( 'Starting evolution', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

!    do while ( I % T  <  I % T_Finish .and. I % iCycle  <  I % FinishCycle )
!      call Show ( 'Computing a cycle', I % IGNORABILITY + 1 )

!       call I % ComputeCycle ( )

!       call Show ( 'Cycle computed', I % IGNORABILITY + 1 )
!       call Show ( I % iCycle, 'iCycle', I % IGNORABILITY + 1 )
!       call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY + 1 )

!       TimeStepRatio  &
!         =  minval ( I % TimeStepCandidate ) &
!              / max ( I % T_CheckpointInterval, sqrt ( tiny ( 0.0_KDR ) ) )
!       if ( TimeStepRatio  <  1.0e-6  *  I % nWrite ) then
!         call I % AdministerCheckpoint ( )
!         call Show ( 'TimeStepRatio too small', CONSOLE % WARNING )
!         call Show ( TimeStepRatio, 'TimeStepRatio', CONSOLE % WARNING )
!         exit
!       end if

! !call I % Write ( )
!       if ( I % IsT_Checkpoint ) &
!         call I % AdministerCheckpoint ( )

!    end do !-- T  <  T_Finish

    call T % Stop ( )   

  end subroutine Evolve


  subroutine Show_I ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( I % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    call Show ( I % Communicator % Name, 'Communicator', I % IGNORABILITY )

    call Show ( I % T_Start, I % Unit_T, 'T_Start', I % IGNORABILITY )
    call Show ( I % T_Finish, I % Unit_T, 'T_Finish', I % IGNORABILITY )

    call Show ( I % n_dT_Candidates, 'n_dT_Candidates', I % IGNORABILITY )
    call Show ( I % dT_Label, 'dT_Label', I % IGNORABILITY )

    call Show ( I % nRampCycles, 'nRampCycles', I % IGNORABILITY )
    call Show ( I % FinishCycle, 'FinishCycle', I % IGNORABILITY )

    call I % ShowManifold ( )
    call I % ShowFields ( )
    call I % ShowCheckpoint ( )

  end subroutine Show_I


  impure elemental subroutine Finalize ( I )

    type ( Integrator_H_Form ), intent ( inout ) :: &
      I

    if ( I % Name == '' ) &
      return

    if ( allocated ( I % Geometry_X_A ) ) &
      deallocate ( I % Geometry_X_A )
    if ( allocated ( I % X_A ) ) &
      deallocate ( I % X_A )
    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )

    nullify ( I % Communicator )

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

  end subroutine Finalize


  subroutine ShowManifold ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % X_A % Show ( )

  end subroutine ShowManifold


  subroutine ShowFields ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % Geometry_X_A % Show ( )

  end subroutine ShowFields


  subroutine ShowCheckpoint ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call Show ( 'Checkpointing', I % IGNORABILITY )

    call Show ( I %   nWrite,   'nWrite', I % IGNORABILITY )
    call Show ( I %  NoWrite,  'NoWrite', I % IGNORABILITY )
    call Show ( I % AllWrite, 'AllWrite', I % IGNORABILITY )
    
    call Show ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval', &
                I % IGNORABILITY )
    call Show ( I % T_CheckpointExact, 'T_CheckpointExact', &
                I % IGNORABILITY )

    call Show ( I % GridImageStream % Name, 'GridImageStream', &
                I % IGNORABILITY )

    associate ( SA  =>  I % Checkpoint_X_A )
    call SA % Show ( )
    end associate !-- SA
  end subroutine ShowCheckpoint


  subroutine PrepareInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      RestartFrom
    type ( MeasuredValueForm ) :: &
      T_Restart

    if ( .not. associated ( I % SetInitial ) ) then
      call Show ( 'SetInitial unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % SetInitial  =>  SetInitial_H
    end if

    if ( .not. associated ( I % ResetInitial ) ) then
      call Show ( 'ResetInitial unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % ResetInitial  =>  ResetInitial_H
    end if

    if ( .not. associated ( I % Write ) ) then
      call Show ( 'Write unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Write  =>  Write_H
    end if

    if ( .not. associated ( I % Read ) ) then
      call Show ( 'Read unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Read  =>  Read_H
    end if

    if ( .not. associated ( I % Set_T_CheckpointInterval ) ) then
      call Show ( 'Set_T_CheckpointInterval unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end if

    RestartFrom  =  - huge ( 1 )
    call PROGRAM_HEADER % GetParameter ( RestartFrom, 'RestartFrom' )

    if ( RestartFrom >= 0 ) then
      call I % ResetInitial ( RestartFrom, T_Restart )
      I % Start    =  .false.
      I % Restart  =  .true.
      I % T        =  T_Restart
    else !-- no restart
      call I % SetInitial ( )
      I % Start    =  .true.
      I % Restart  =  .false.
      I % T        =  I % T_Start
    end if !-- restart

  end subroutine PrepareInitial


  subroutine PrepareEvolution ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine PrepareEvolution


  subroutine AdministerCheckpoint ( I, ComputeChangeOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    integer ( KDI ) :: &
      iTSC, &  !-- iTimeStepCandidates
!       TallyIgnorability, &
      StatisticsIgnorability
    real ( KDR ), dimension ( : ), allocatable :: &
      MaxTime, &
      MinTime, &
      MeanTime
!     logical ( KDL ) :: &
!       WriteSeries
    type ( TimerForm ), pointer :: &
      T_AC!, &
!       Timer_T, &
!       Timer_A, &
!       Timer_WS
    
!     Timer_T  => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
!     Timer_A  => PROGRAM_HEADER % TimerPointer ( I % iTimerAnalyze )
!     Timer_WS => PROGRAM_HEADER % TimerPointer ( I % iTimerWriteSeries )

    associate ( iT  =>  I % iTimer_AC )
    if ( iT == 0 ) then
      call PROGRAM_HEADER % AddTimer ( 'AdministerCheckpoint', iT, Level = 2 )
    end if
    end associate !-- iT
    T_AC  =>  PROGRAM_HEADER % TimerPointer ( I % iTimer_AC )

    call T_AC % Start ( )   

    call Show ( 'Checkpoint reached', I % IGNORABILITY )
    call Show ( I % iCheckpoint, 'iCheckpoint', I % IGNORABILITY )
    call Show ( I % iCycle, 'iCycle', I % IGNORABILITY )
    call Show ( I % T, I % Unit_T, 'T', I % IGNORABILITY )
    if ( .not. I % Start .and. .not. I % Restart ) then
      do iTSC = 1, I % n_dT_Candidates
        call Show ( I % dT_Candidate ( iTSC ), I % Unit_T, &
                    trim ( I % dT_Label ( iTSC ) ) // ' dT', &
                    I % IGNORABILITY )
      end do !-- iTSC
    end if

    call I % UpdateHost ( TimerLevelOption  =  T_AC % Level  +  1 )

!     WriteSeries = .true.

    if ( .not. I % Start & !.and. .not. I % Restart &
         .and. I % T  <  I % T_Finish &
         .and. mod ( I % iCheckpoint, I % CheckpointDisplayInterval ) > 0 ) &
    then
!       TallyIgnorability       =  I % IGNORABILITY + 2
      StatisticsIgnorability  =  I % IGNORABILITY + 2
! !      WriteSeries = .false.
    else
!       TallyIgnorability       =  CONSOLE % INFO_1
      StatisticsIgnorability  =  CONSOLE % INFO_1
! !      WriteSeries = .true.
    end if

!     if ( associated ( Timer_T ) ) call Timer_T % Start ( )   
!     call I % ComputeTally &
!            ( ComputeChangeOption = ComputeChangeOption, &
!              IgnorabilityOption  = TallyIgnorability )
!     if ( associated ( Timer_T ) ) call Timer_T % Stop ( )   

!     if ( associated ( Timer_A ) ) call Timer_A % Start ( )   
!     call I % Analyze ( )
!     if ( associated ( Timer_A ) ) call Timer_A % Stop ( )   

    if ( .not. I % NoWrite .and. .not. I % Restart ) &
      call I % Write ( TimerLevelOption  =  T_AC % Level  +  1 )

    associate ( nT  =>  PROGRAM_HEADER % nTimers )
    allocate ( MaxTime ( nT ), MinTime ( nT ), MeanTime ( nT ) )
    call PROGRAM_HEADER % ShowStatistics &
           ( StatisticsIgnorability, &
             CommunicatorOption = PROGRAM_HEADER % Communicator, &
             MaxTimeOption = MaxTime, MinTimeOption = MinTime, &
             MeanTimeOption = MeanTime )
    end associate !-- nT

!     if ( .not. I % Restart ) &
!       call I % RecordTimeSeries ( MaxTime, MinTime, MeanTime )

!     if ( associated ( Timer_WS ) ) call Timer_WS % Start ( )   
!     if ( WriteSeries .and. .not. I % NoWrite .and. .not. I % Restart ) &
!       call I % WriteTimeSeries ( )
!     if ( associated ( Timer_WS ) ) call Timer_WS % Stop ( )   

    I % CheckpointDue  =  .false.
    if ( I % T  <  I % T_Finish ) then
      call I % Set_T_CheckpointInterval ( )
      I % T_Checkpoint &
        =  min ( I % T  +  I % T_CheckpointInterval, I % T_Finish )
      if ( I % T_Checkpoint  ==  I % T_Finish ) &
        I % T_CheckpointExact  =  .true.
      call Show ( I % T_CheckpointInterval, I % Unit_T, &
                  'T_CheckpointInterval', &
                  I % IGNORABILITY )
      call Show ( I % T_Checkpoint, I % Unit_T, 'Next T_Checkpoint', &
                  I % IGNORABILITY )
    else 
      call Show ( 'T_Finish reached', I % IGNORABILITY )
    end if  !-- T  <  T_Finish

    I % iCheckpoint  =  I % iCheckpoint + 1
    
    I % Start    =  .false.
    I % Restart  =  .false.

    call T_AC % Stop ( )

  end subroutine AdministerCheckpoint


  subroutine SetInitial_H ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine SetInitial_H


  subroutine ResetInitial_H ( I, RestartFrom, T_Restart )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      RestartFrom
    type ( MeasuredValueForm ), intent ( out ) :: &
      T_Restart

    integer ( KDI ) :: &
      CycleNumber
    ! real ( KDR ), dimension ( PROGRAM_HEADER % nTimers ) :: &
    !   MaxTime, &
    !   MinTime, &
    !   MeanTime

    call I % Read ( RestartFrom, T_Restart, CycleNumber )

    I % iCheckpoint  =  RestartFrom
    I % iCycle       =  CycleNumber  !-- needed by RestoreTimeSeries

    ! call I % ReadTimeSeries ( nSeries = RestartFrom + 1 )

    ! call I % RestoreTimeSeries ( MaxTime, MinTime, MeanTime )

    ! call PROGRAM_HEADER % RestoreStatistics &
    !        ( Ignorability = CONSOLE % INFO_1, &
    !          CommunicatorOption = PROGRAM_HEADER % Communicator, &
    !          MeanTimeOption = MeanTime )
 
  end subroutine ResetInitial_H


  subroutine UpdateHost_H ( I, TimerLevelOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    associate ( GA  =>  I % Geometry_X_A )
    call GA % UpdateHost ( TimerLevelOption )
    end associate !-- GA

  end subroutine UpdateHost_H


  subroutine Write_H ( I, TimerLevelOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % MarkFibersWritten ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    associate ( GIS => I % GridImageStream )
    call GIS % Open ( GIS % ACCESS_CREATE )

    associate ( SA  =>  I % Checkpoint_X_A )
    call SA % Write &
           ( TimeOption  =  I % T  /  I % Unit_T, &
             CycleNumberOption  =  I % iCycle, &
             TimerLevelOption  =  TimerLevelOption )
    end associate !-- SA

    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call GIS % Close ( )

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % Write &
    !            ( iStream = iS, TimeOption = I % Time / I % TimeUnit, &
    !              CycleNumberOption = I % iCycle )
    !   class default
    !     call Show ( 'Bundle type not found', CONSOLE % ERROR )
    !     call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
    !     call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
    !     call PROGRAM_HEADER % Abort ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    end associate !-- GIS

  end subroutine Write_H


  subroutine Read_H ( I, ReadFrom, T, CycleNumber )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      ReadFrom
    type ( MeasuredValueForm ), intent ( out ) :: &
      T
    integer ( KDI ), intent ( out ) :: &
      CycleNumber

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % MarkFibersWritten ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    associate ( GIS => I % GridImageStream )
    call GIS % Open ( GIS % ACCESS_READ, NumberOption = ReadFrom )

    associate ( SA  =>  I % Checkpoint_X_A )
    call SA % Read &
           ( TimeOption = T, &
             CycleNumberOption = CycleNumber )
    T  =  T  *  I % Unit_T
    end associate !-- SA

    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call GIS % Close ( )

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % Write &
    !            ( iStream = iS, TimeOption = I % T / I % TimeUnit, &
    !              CycleNumberOption = I % iCycle )
    !   class default
    !     call Show ( 'Bundle type not found', CONSOLE % ERROR )
    !     call Show ( 'Integrator_Template', 'module', CONSOLE % ERROR )
    !     call Show ( 'Write', 'subroutine', CONSOLE % ERROR ) 
    !     call PROGRAM_HEADER % Abort ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    end associate !-- GIS

  end subroutine Read_H


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    I % T_CheckpointInterval &
      =  ( I % T_Finish  -  I % T_Start )  /  I % nWrite

  end subroutine Set_T_CheckpointInterval


end module Integrator_H__Form
