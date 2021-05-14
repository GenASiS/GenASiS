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
      ! iTimerEvolve = 0, &
      ! iTimerCycle = 0, &
      ! iTimerNewTime = 0, &
      ! iTimerCheckpoint = 0, &
      ! iTimerTally = 0, &
      ! iTimerAnalyze = 0, &
      ! iTimerWrite = 0, &
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
    !   CheckpointTimeInterval, &
    !   CheckpointTime, &
      Time
    type ( MeasuredValueForm ) :: &
      Unit_T
    real ( KDR ), dimension ( : ), allocatable :: &
      dT_Candidate
    logical ( KDL ) :: &
    !   Start, &
    !   Restart, &
    !   IsCheckpointTime, &
      NoWrite, &
      AllWrite, &
      CheckpointTimeExact
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
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H
    generic, public :: &           
      Initialize => Initialize_H
    procedure, public, pass :: &   !-- 1
      Evolve
    procedure, private, pass :: &  !-- 1
      Show_I
    generic, public :: &
      Show => Show_I
    final :: &                     !-- 1
      Finalize
    procedure, private, pass :: &   !-- 2
      ShowManifold
    procedure, private, pass :: &   !-- 2
      ShowFields
    procedure, private, pass :: &   !-- 2
      ShowCheckpoint
  end type Integrator_H_Form

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
      I % dT_Label ( 1 ) = 'dT_Candidate'
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
    I % CheckpointTimeExact  =  .false.
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval' )
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointTimeExact, 'CheckpointTimeExact' )

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
!     type ( TimerForm ), pointer :: &
!       Timer

!     call I % OpenManifoldStreams ( )
!     call I % InitializeTimers ( )
!     call I % InitializeTimeSeries ( )

!     Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerEvolve )
!     if ( associated ( Timer ) ) call Timer % Start ( )   

!     call I % PrepareInitial ( )
!     call I % PrepareEvolution ( )
!     call I % AdministerCheckpoint ( ComputeChangeOption = .false. )

!     call Show ( 'Starting evolution', I % IGNORABILITY )
!     call Show ( I % Name, 'Name', I % IGNORABILITY )

!     do while ( I % Time < I % FinishTime .and. I % iCycle < I % FinishCycle )
!       call Show ( 'Computing a cycle', I % IGNORABILITY + 1 )

!       call I % ComputeCycle ( )

!       call Show ( 'Cycle computed', I % IGNORABILITY + 1 )
!       call Show ( I % iCycle, 'iCycle', I % IGNORABILITY + 1 )
!       call Show ( I % Time, I % TimeUnit, 'Time', I % IGNORABILITY + 1 )

!       TimeStepRatio  &
!         =  minval ( I % TimeStepCandidate ) &
!              / max ( I % CheckpointTimeInterval, sqrt ( tiny ( 0.0_KDR ) ) )
!       if ( TimeStepRatio  <  1.0e-6  *  I % nWrite ) then
!         call I % AdministerCheckpoint ( )
!         call Show ( 'TimeStepRatio too small', CONSOLE % WARNING )
!         call Show ( TimeStepRatio, 'TimeStepRatio', CONSOLE % WARNING )
!         exit
!       end if

! !call I % Write ( )
!       if ( I % IsCheckpointTime ) &
!         call I % AdministerCheckpoint ( )

!     end do !-- Time < FinishTime 

!     if ( associated ( Timer ) ) call Timer % Stop ( )   

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
    call Show ( I % CheckpointTimeExact, 'CheckpointTimeExact', &
                I % IGNORABILITY )

    call Show ( I % GridImageStream % Name, 'GridImageStream', &
                I % IGNORABILITY )

    associate ( SA  =>  I % Checkpoint_X_A )
    call SA % Show ( )
    end associate !-- SA
  end subroutine ShowCheckpoint


end module Integrator_H__Form
