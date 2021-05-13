module Integrator_H__Form

  !-- Integrator_Header_Form

  use Basics
  use Manifolds

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
      CheckpointTimeExact
    character ( LDL ), dimension ( : ), allocatable :: &
      dT_Label
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    ! type ( GridImageStreamForm ), allocatable :: &
    !   GridImageStream
    class ( Atlas_H_Form ), allocatable :: &
      X
  contains
    procedure, private, pass :: &
      Initialize_H
    generic, public :: &
      Initialize => Initialize_H
    procedure, private, pass :: &
      Show_I
    generic, public :: &
      Show => Show_I
    final :: &
      Finalize
  end type Integrator_H_Form

contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, Unit_T_Option, &
                 T_FinishOption, nWriteOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      Unit_T_Option
    real ( KDR ), intent ( in ), optional :: &
      T_FinishOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    I % IGNORABILITY  =  CONSOLE % INFO_1

    if ( I % Type == '' ) &
      I % Type = 'an Integrator' 

    I % Name = 'Integrator'
    if ( present ( NameOption ) ) &
      I % Name  =  NameOption

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

    I % nWrite  =  100
    if ( present ( nWriteOption ) ) &
      I % nWrite  =  nWriteOption
    I % NoWrite  =  .false.
    call PROGRAM_HEADER % GetParameter ( I % nWrite, 'nWrite' )
    call PROGRAM_HEADER % GetParameter ( I % NoWrite, 'NoWrite' )

    I % CheckpointDisplayInterval  =  100
    I % CheckpointTimeExact  =  .false.
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval' )
    call PROGRAM_HEADER % GetParameter &
           ( I % CheckpointTimeExact, 'CheckpointTimeExact' )

    if ( present ( CommunicatorOption ) ) then
      I % Communicator  =>  CommunicatorOption
    else
      I % Communicator  =>  PROGRAM_HEADER % Communicator
    end if

    if ( .not. allocated ( I % X ) ) then
      allocate ( Atlas_SCG_Form :: I % X )
      select type ( A  =>  I % X )
        class is ( Atlas_SCG_Form )
      call A % Initialize &
             ( CommunicatorOption = I % Communicator, &
               NameOption = 'X', &
               PeriodicOption = [ .true., .true., .true. ] )
      end select !-- A
    end if

  end subroutine Initialize_H


  subroutine Show_I ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( I % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    call Show ( I % T_Start, I % Unit_T, 'T_Start', I % IGNORABILITY )
    call Show ( I % T_Finish, I % Unit_T, 'T_Finish', I % IGNORABILITY )

    call Show ( I % n_dT_Candidates, 'n_dT_Candidates', I % IGNORABILITY )
    call Show ( I % dT_Label, 'dT_Label', I % IGNORABILITY )

    call Show ( I % nRampCycles, 'nRampCycles', I % IGNORABILITY )
    call Show ( I % FinishCycle, 'FinishCycle', I % IGNORABILITY )

    call Show ( I % nWrite, 'nWrite', I % IGNORABILITY )
    call Show ( I % NoWrite, 'NoWrite', I % IGNORABILITY )
    
    call Show ( I % CheckpointDisplayInterval, 'CheckpointDisplayInterval', &
                I % IGNORABILITY )
    call Show ( I % CheckpointTimeExact, 'CheckpointTimeExact', &
                I % IGNORABILITY )

    call Show ( I % Communicator % Name, 'Communicator', I % IGNORABILITY )

    call I % X % Show ( )

 end subroutine Show_I


  impure elemental subroutine Finalize ( I )

    type ( Integrator_H_Form ), intent ( inout ) :: &
      I

    if ( I % Name == '' ) &
      return

    if ( allocated ( I % X ) ) &
      deallocate ( I % X )

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

!    if ( allocated ( I % GridImageStream ) ) &
!      deallocate ( I % GridImageStream )

!    nullify ( I % Communicator )

  end subroutine Finalize


end module Integrator_H__Form
