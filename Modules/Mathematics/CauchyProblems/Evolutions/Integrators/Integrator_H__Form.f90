module Integrator_H__Form

  !-- Integrator_Header_Form

  use Basics
  use Manifolds
  use Fields
  use Steps

  implicit none
  private

  type, public :: Integrator_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iCycle, &
      iCheckpoint, &
      nRampCycles, &
      FinishCycle, &
      nWrite, &
      n_dT_Candidates, &
      CheckpointDisplayInterval
    integer ( KDI ) :: &
      iTimer_E   = 0, &   !-- Evolution
      iTimer_AC  = 0, &   !-- AdministerCheckpoint
      iTimer_UH  = 0, &   !-- UpdateHost
      iTimer_A   = 0, &   !-- Analyze
      iTimer_W   = 0, &   !-- Write
      iTimer_CC  = 0, &   !-- ComputeCycle
      iTimer_PC  = 0, &   !-- PrepareCycle
      iTimer_CTN = 0!, &  !-- Compute_T_New
      ! iTimerTally = 0, &
      ! iTimerWriteSeries = 0, &
    real ( KDR ) :: &
      T_Start               = 0.0_KDR, &
      T_Finish              = 1.0_KDR, &
      T_CheckpointInterval  = 0.0_KDR, &
      T_Checkpoint          = 0.0_KDR, & 
      T                     = 0.0_KDR
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
      X
    type ( StreamForm ), allocatable :: &
      Checkpoint_X
    class ( Geometry_F_Form ), allocatable :: &
      Geometry_X
    class ( Step_RK_H_Form ), allocatable :: &
      Step_X
    class ( * ), pointer :: &
      System => null ( )
    procedure ( SI ), public, pointer :: &
      SetInitial => null ( )
    procedure ( RI ), public, pointer :: &
      ResetInitial => null ( )
    !-- FIXME: renamed SS to SS_I to work around CCE bug
    !procedure ( SS ), public, pointer :: &
    procedure ( SS_I ), public, pointer :: &
      ShowSystem => null ( )
    procedure ( A ), public, pointer :: &
      Analyze => null ( )
    procedure ( W ), public, pointer :: &
      Write => null ( )
    procedure ( R ), public, pointer :: &
      Read => null ( )
    procedure ( SR ), public, pointer :: &
      SetReference => null ( )
    procedure ( STCI ), pointer :: &
      Set_T_CheckpointInterval => null ( )
    procedure ( C_dT_L ), pointer :: &
      Compute_dT_Local => null ( )
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H
    generic, public :: &           
      Initialize => Initialize_H
    procedure, public, pass :: &
      Show => Show_I
    procedure, public, pass :: &   !-- 1
      Evolve
    final :: &                     !-- 1
      Finalize
    procedure, public, pass :: &  !-- 2
      ShowParameters
    procedure, public, pass :: &  !-- 2
      ShowManifold
    procedure, public, pass :: &  !-- 2
      ShowFields
    procedure, public, pass :: &  !-- 2
      ShowSteps
    procedure, public, pass :: &  !-- 2
      ShowCheckpointing
    procedure, private, pass :: &
      Timer_E
    procedure, private, pass :: &
      Timer_AC
    procedure, private, pass :: &
      Timer_UH
    procedure, public, pass :: &
      Timer_A
    procedure, public, pass :: &
      Timer_W
    procedure, private, pass :: &
      Timer_CC
    procedure, private, pass :: &
      Timer_PC
    procedure, private, pass :: &
      Timer_CTN
    procedure, private, pass :: &   !-- 2
      PrepareInitial
    procedure, private, pass :: &   !-- 2
      PrepareEvolution
    procedure, private, pass :: &  !-- 2
      AdministerCheckpoint
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &   !-- 3
      SetInitial_H
    procedure, public, pass :: &   !-- 3
      ResetInitial_H
    procedure, private, pass :: &  !-- 3
      ShowSystem_H
    procedure, public, pass :: &   !-- 3
      UpdateHost => UpdateHost_H
    procedure, public, pass :: &   !-- 3
      Analyze_H
    procedure, private, pass :: &   !-- 3
      Write_H
    procedure, private, pass :: &   !-- 3
      Read_H
    procedure, private, pass :: &   !-- 3
      PrepareCycle
    procedure, private, pass :: &   !-- 3
      Compute_T_New
    procedure, private, pass :: &   !-- 4
      Compute_dT
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

    subroutine SS_I ( I )
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( in ) :: &
        I
    end subroutine SS_I

    subroutine A ( I, T_A )
      use Basics
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      type ( TimerForm ), intent ( in ) :: &
        T_A
    end subroutine A
    
    subroutine W ( I, T_W )
      use Basics
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      type ( TimerForm ), intent ( in ) :: &
        T_W
    end subroutine W

    subroutine R ( I, ReadFrom, T, CycleNumber )
      use Basics
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
      integer ( KDI ), intent ( in ) :: &
        ReadFrom
      type ( MeasuredValueForm ), intent ( out ) :: &
        T
      integer ( KDI ), intent ( out ) :: &
        CycleNumber
    end subroutine R

    subroutine SR ( I )
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
    end subroutine SR
    
    subroutine STCI ( I )
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ) :: &
        I
    end subroutine STCI

    subroutine C_dT_L ( I, dT_Candidate, iC, T_Option )
      use Basics
      import Integrator_H_Form
      implicit none
      class ( Integrator_H_Form ), intent ( inout ), target :: &
        I
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        dT_Candidate
      integer ( KDI ), intent ( in ) :: &
        iC
      type ( TimerForm ), intent ( in ), optional :: &
        T_Option
    end subroutine C_dT_L

  end interface

  
    private :: &
      Set_T_CheckpointInterval, &
      Compute_dT_Local
          

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

    !-- Atlas, if necessary

    if ( .not. allocated ( I % X ) ) then
      allocate ( Atlas_SCG_Form :: I % X )
      select type ( A_X  =>  I % X )
        class is ( Atlas_SCG_Form )
      call A_X % Initialize &
             ( CommunicatorOption = I % Communicator, &
               NameOption = 'X' )
      end select !-- A_X
    end if

    !-- Geometry, if necessary

    if ( .not. allocated ( I % Geometry_X ) ) then
      allocate ( I % Geometry_X )
      associate ( G_X  =>  I % Geometry_X )
      call G_X % Initialize &
             ( I % X, &
               DeviceMemoryOption = DeviceMemoryOption, &
               PinnedMemoryOption = PinnedMemoryOption, &
               DevicesCommunicateOption = DevicesCommunicateOption )
      end associate !-- G_X
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

    !-- Step, if necessary (only for Integrator_H_Form test)

    if ( .not. allocated ( I % Step_X ) ) then
      call Show ( 'Step_X not allocated', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      allocate ( I % Step_X )
      associate ( S  =>  I % Step_X )
      call S % Initialize_H ( I % X )
      end associate !-- S
    end if

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

    allocate ( I % Checkpoint_X )
    associate &
      ( S_X  =>  I % Checkpoint_X, &
        A_X  =>  I % X, &
        GIS  =>  I % GridImageStream )
    call S_X % Initialize ( A_X, GIS, NameOption = A_X % Name )
    call I % Geometry_X % SetStream ( S_X )
    end associate !--SA
    
  end subroutine Initialize_H


  subroutine Evolve ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    real ( KDR ) :: &
      dT_Ratio
    type ( TimerForm ), pointer :: &
      T_E, &
      T_AC, &
      T_CC, &
      T_A, &
      T_W

!     call I % InitializeTimeSeries ( )

    call I % PrepareInitial ( )
    call I % PrepareEvolution ( )

    call Show ( 'Starting evolution', I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

    T_E   =>  I % Timer_E  ( LevelOption = 1 )
    call T_E % Start ( )

    T_AC  =>  I % Timer_AC ( LevelOption = T_E % Level + 1 )
    call T_AC % Start ( )
    call I % AdministerCheckpoint ( T_AC, ComputeChangeOption = .false. )
    call T_AC % Stop ( )

    do while ( I % T  <  I % T_Finish .and. I % iCycle  <  I % FinishCycle )
      call Show ( 'Computing a cycle', I % IGNORABILITY + 1 )

      T_CC  =>  I % Timer_CC ( LevelOption = T_E % Level + 1 )
      call T_CC % Start ( )
      call I % ComputeCycle ( T_CC )
      call T_CC % Stop ( )

      call Show ( 'Cycle computed', I % IGNORABILITY + 1 )
      call Show ( I % iCycle, 'iCycle', I % IGNORABILITY + 1 )
      call Show ( I % T, I % Unit_T, 'T', I % IGNORABILITY + 1 )

      dT_Ratio  &
        =  minval ( I % dT_Candidate ) &
             /  max ( I % T_CheckpointInterval, sqrt ( tiny ( 0.0_KDR ) ) )
      if ( dT_Ratio  <  1.0e-6  *  I % nWrite ) then
        call T_AC % Start ( )
        call I % AdministerCheckpoint ( T_AC )
        call T_AC % Stop ( )
        call Show ( '*** dT_Ratio too small', CONSOLE % WARNING, &
                    nLeadingLinesOption = 2 )
        call Show ( dT_Ratio, 'dT_Ratio', CONSOLE % WARNING, &
                    nTrailingLinesOption = 2 )
        exit
      end if

      if ( I % AllWrite  .and. .not. I % CheckpointDue ) then

        T_A   =>  I % Timer_A  ( )
        call T_A % Start ( )
        call I % Analyze ( T_A )
        call T_A % Stop ( )

        T_W   =>  I % Timer_W ( )
        call T_W % Start ( )
        call I % Write ( T_W )
        call T_W % Stop ( )

      end if

      if ( I % CheckpointDue ) then
        call T_AC % Start ( )
        call I % AdministerCheckpoint ( T_AC )
        call T_AC % Stop ( )
      end if

    end do !-- T  <  T_Finish

    call T_E % Stop ( )   

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

    call I % ShowParameters ( )
    call I % ShowManifold ( )
    call I % ShowFields ( )
    call I % ShowSteps ( )
    call I % ShowCheckpointing ( )

  end subroutine Show_I


  impure elemental subroutine Finalize ( I )

    type ( Integrator_H_Form ), intent ( inout ) :: &
      I

    if ( I % Name == '' ) &
      return

    if ( allocated ( I % Step_X ) ) &
      deallocate ( I % Step_X )
    if ( allocated ( I % Geometry_X ) ) &
      deallocate ( I % Geometry_X )
    if ( allocated ( I % Checkpoint_X ) ) &
      deallocate ( I % Checkpoint_X )
    if ( allocated ( I % X ) ) &
      deallocate ( I % X )
    if ( allocated ( I % GridImageStream ) ) &
      deallocate ( I % GridImageStream )

    nullify ( I % System )
    nullify ( I % Communicator )

    call Show ( 'Finalizing ' // trim ( I % Type ), I % IGNORABILITY )
    call Show ( I % Name, 'Name', I % IGNORABILITY )

  end subroutine Finalize


  subroutine ShowParameters ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call Show ( I % T_Start, I % Unit_T, 'T_Start', I % IGNORABILITY )
    call Show ( I % T_Finish, I % Unit_T, 'T_Finish', I % IGNORABILITY )

    call Show ( I % n_dT_Candidates, 'n_dT_Candidates', I % IGNORABILITY )
    call Show ( I % dT_Label, 'dT_Label', I % IGNORABILITY )

    call Show ( I % nRampCycles, 'nRampCycles', I % IGNORABILITY )
    call Show ( I % FinishCycle, 'FinishCycle', I % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowManifold ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % X % Show ( )

  end subroutine ShowManifold


  subroutine ShowFields ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % Geometry_X % Show ( )

  end subroutine ShowFields


  subroutine ShowSteps ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % Step_X % Show ( )

  end subroutine ShowSteps


  subroutine ShowCheckpointing ( I )

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

    associate ( SA  =>  I % Checkpoint_X )
    call SA % Show ( )
    end associate !-- SA

  end subroutine ShowCheckpointing


  function Timer_E ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_E )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_Evltn'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_E


  function Timer_AC ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_AC )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_Chckpnt'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_AC


  function Timer_UH ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_UH )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_UpdtHst'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_UH


  function Timer_A ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_A )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_Anlz'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_A


  function Timer_W ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_W )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_Wrt'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_W


  function Timer_CC ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_CC )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_CmptCcl'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_CC


  function Timer_PC ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_PC )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_PrprCcl'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_PC


  function Timer_CTN ( I, LevelOption ) result ( T )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  I % iTimer_CTN )

    if ( iT == 0 ) then
      TimerName  =  trim ( I % Name ) // '_CmptTmStp'
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_CTN


  subroutine PrepareInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      RestartFrom
    type ( MeasuredValueForm ) :: &
      T_Restart

    if ( .not. associated ( I % SetInitial ) ) then
      call Show ( 'SetInitial method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % SetInitial  =>  SetInitial_H
    end if

    if ( .not. associated ( I % ResetInitial ) ) then
      call Show ( 'ResetInitial method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % ResetInitial  =>  ResetInitial_H
    end if

    if ( .not. associated ( I % ShowSystem ) ) then
      call Show ( 'ShowSystem method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % ShowSystem  =>  ShowSystem_H
    end if

    if ( .not. associated ( I % Analyze ) ) then
      call Show ( 'Analyze method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Analyze  =>  Analyze_H
    end if

    if ( .not. associated ( I % Write ) ) then
      call Show ( 'Write method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Write  =>  Write_H
    end if

    if ( .not. associated ( I % Read ) ) then
      call Show ( 'Read method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Read  =>  Read_H
    end if

    if ( .not. associated ( I % Set_T_CheckpointInterval ) ) then
      call Show ( 'Set_T_CheckpointInterval method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    end if

    if ( .not. associated ( I % Compute_dT_Local ) ) then
      call Show ( 'Compute_dT_Local method unset', CONSOLE % WARNING )
      call Show ( 'Integrator_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'PrepareInitial', 'subroutine', CONSOLE % WARNING )
      I % Compute_dT_Local  =>  Compute_dT_Local
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

    call I % ShowSystem ( )

  end subroutine PrepareInitial


  subroutine PrepareEvolution ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine PrepareEvolution


  subroutine AdministerCheckpoint ( I, T_AC, ComputeChangeOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ), optional :: &
      T_AC
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
      T_UH, &
      T_A, &
      T_W
!       Timer_T, &
!       Timer_WS
    
!     Timer_T  => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
!     Timer_WS => PROGRAM_HEADER % TimerPointer ( I % iTimerWriteSeries )

      T_UH  =>  I % Timer_UH ( LevelOption = T_AC % Level + 1 )
      T_A   =>  I % Timer_A  ( LevelOption = T_AC % Level + 1 )
      T_W   =>  I % Timer_W  ( LevelOption = T_AC % Level + 1 )

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

    call T_UH % Start ( )   
    call I % UpdateHost ( )
    call T_UH % Stop ( )

!     WriteSeries = .true.

    if ( .not. I % Start & !.and. .not. I % Restart &
         .and. I % T  <  I % T_Finish &
         .and. mod ( I % iCheckpoint, I % CheckpointDisplayInterval ) > 0 ) &
    then
!       TallyIgnorability       =  I % IGNORABILITY + 2
      StatisticsIgnorability  =  I % IGNORABILITY + 2
!      WriteSeries = .false.
    else
!       TallyIgnorability       =  CONSOLE % INFO_1
      StatisticsIgnorability  =  CONSOLE % INFO_1
!      WriteSeries = .true.
    end if

!     if ( associated ( Timer_T ) ) call Timer_T % Start ( )   
!     call I % ComputeTally &
!            ( ComputeChangeOption = ComputeChangeOption, &
!              IgnorabilityOption  = TallyIgnorability )
!     if ( associated ( Timer_T ) ) call Timer_T % Stop ( )   

    call T_A % Start ( )   
    call I % Analyze ( T_A )
    call T_A % Stop ( )   

    if ( .not. I % NoWrite .and. .not. I % Restart ) then
      call T_W % Start ( )
      call I % Write ( T_W )
      call T_W % Stop ( )
    end if

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
      if ( abs ( I % T_Finish  -  I % T_Checkpoint )  /  I % T_Finish  &
                 <  1.0e-14 ) &
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

  end subroutine AdministerCheckpoint


  subroutine ComputeCycle ( I, T_CC )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_CC

    real ( KDR ) :: &
      T_New  !-- Use of Compute_T_New is a relic of past AMR evolution
    type ( TimerForm ), pointer :: &
      T_PC, &
      T_CTN, &
      T_S

    T_PC  =>  I % Timer_PC ( LevelOption = T_CC % Level + 1 )
    call T_PC % Start ( )
    call I % PrepareCycle ( )
    call T_PC % Stop ( )

    T_CTN  =>  I % Timer_CTN ( LevelOption = T_CC % Level + 1 )
    call T_CTN % Start ( )
    call I % Compute_T_New ( T_New )
    call T_CTN % Stop ( )

    associate (  S  =>  I % Step_X )
    associate ( dT  =>  T_New  -  I % T )    

    ! select type ( Chart => PS % Chart )
    ! class is ( Chart_SLD_Form )

    T_S  =>  S % Timer ( LevelOption = T_CC % Level + 1 )
    call T_S % Start ( )
    call S % Compute ( I % T, dT, T_Option = T_S )
    call T_S % Stop ( )

    ! class default
    !   call Show ( 'Chart type not found', CONSOLE % ERROR )
    !   call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % ERROR )
    !   call Show ( 'ComputeCycle_ASC', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end select !-- C

    I % iCycle  =  I % iCycle  +   1
    I % T       =  I % T       +  dT

    if ( I % T_CheckpointExact ) then
      if ( abs ( I % T_Checkpoint  -  I % T )  /  I % T_Checkpoint  &
           <  1.0e-14 ) &
        I % CheckpointDue  =  .true.
    else 
      if ( I % T  >  I % T_Checkpoint &
           .or. abs ( I % T_Checkpoint  -  I % T )  <  0.5_KDR * dT ) &
        I % CheckpointDue  =  .true.
    end if

    end associate !-- dT
    end associate !--  S

  end subroutine ComputeCycle


  subroutine SetInitial_H ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine SetInitial_H


  subroutine ResetInitial_H ( I, RestartFrom, T_Restart  )

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
    
    call Show ( 'Restarting', I % IGNORABILITY )
    call Show ( RestartFrom, 'RestartFrom', I % IGNORABILITY )
    call Show ( T_Restart, I % Unit_T, 'T_Restart', I % IGNORABILITY )
    call Show ( CycleNumber, 'CycleNumber', I % IGNORABILITY )

    I % iCheckpoint  =  RestartFrom
    I % iCycle       =  CycleNumber  !-- needed by RestoreTimeSeries

    ! call I % ReadTimeSeries ( nSeries = RestartFrom + 1 )

    ! call I % RestoreTimeSeries ( MaxTime, MinTime, MeanTime )

    ! call PROGRAM_HEADER % RestoreStatistics &
    !        ( Ignorability = CONSOLE % INFO_1, &
    !          CommunicatorOption = PROGRAM_HEADER % Communicator, &
    !          MeanTimeOption = MeanTime )
 
  end subroutine ResetInitial_H


  subroutine ShowSystem_H ( I )

    class ( Integrator_H_Form ), intent ( in ) :: &
      I

    call I % Show ( )

  end subroutine ShowSystem_H


  subroutine UpdateHost_H ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    associate ( GA  =>  I % Geometry_X )
    call GA % UpdateHost ( )
    end associate !-- GA

  end subroutine UpdateHost_H


  subroutine Analyze_H ( I, T_A )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_A

    type ( TimerForm ), pointer :: &
      T_SS

    if ( associated ( I % SetReference ) ) &
      call I % SetReference ( )

    associate ( Sp_X  =>  I % Step_X )
    T_SS  =>  Sp_X % TimerSlopeSum ( LevelOption = T_A % Level + 1 )
    call T_SS % Start ( )
    call Sp_X % ComputeSlopeSum ( )
    call T_SS % Stop ( )
    end associate !-- Sp_X

  end subroutine Analyze_H


  subroutine Write_H ( I, T_W )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_W

    type ( TimerForm ), pointer :: &
      T_X

    ! if ( allocated ( I % MomentumSpace ) ) then
    !   select type ( MS => I % MomentumSpace )
    !   class is ( Bundle_SLL_ASC_CSLD_Form )
    !     call MS % MarkFibersWritten ( )
    !   end select !-- MS
    ! end if !-- MomentumSpace

    associate ( GIS => I % GridImageStream )
    associate ( S_X  =>  I % Checkpoint_X )
    T_X  =>  S_X % TimerWrite ( LevelOption = T_W % Level + 1 )
    call T_X % Start ( )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call S_X % Write &
           ( TimeOption  =  I % T  /  I % Unit_T, &
             CycleNumberOption  =  I % iCycle )
    !-- Base's GIS must be closed before call to Bundle % Write ( ).
    call GIS % Close ( )

    call T_X % Stop ( )
    end associate !-- S_X
    end associate !-- GIS

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

    associate ( SA  =>  I % Checkpoint_X )
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


  subroutine PrepareCycle ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

  end subroutine PrepareCycle


  subroutine Compute_T_New ( I, T_New, HoldCheckpointSolveOption )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      T_New
    logical ( KDL ), dimension ( : ), intent ( inout ), optional :: &
      HoldCheckpointSolveOption

    real ( KDR ) :: &
      dT

!    call Show ( 'Computing T_New', I % IGNORABILITY )
!    call Show ( I % Name, 'Name', I % IGNORABILITY )

    ! associate ( CFC => I % ConservedFields % Chart ( 1 ) % Element )
    ! select type ( C => I % Atlas % Chart ( 1 ) % Element )
    ! class is ( Chart_SL_Template )
    !   call I % Compute_dT_CSL ( CFC, C, dT )
    ! ! class is ( ChartMultiLevelTemplate )
    ! !   call I % Compute_dT_CML ( CFC, C, dT )
    ! end select !-- C
    ! end associate !-- CFC

    call I % Compute_dT ( dT )

!    associate ( C => I % Atlas % Chart ( 1 ) % Element )

    if ( I % T_CheckpointExact ) then
      if ( I % T  +  dT  >  I % T_Checkpoint ) then
        call Show ( 'T_Checkpoint encountered', I % IGNORABILITY ) 
!        if ( present ( HoldCheckpointSolveOption ) ) &
!          HoldCheckpointSolveOption ( 2 : C % nLevels ) = .true.
        dT  =  I % T_Checkpoint  -  I % T
        call Show ( dT, I % Unit_T, 'Modified dT', I % IGNORABILITY + 1 )
      end if
    end if

    T_New  =  I % T  +  dT
    call Show ( T_New, I % Unit_T, 'T_New', I % IGNORABILITY + 1 )

!    end associate !-- C

  end subroutine Compute_T_New


  subroutine Compute_dT ( I, dT, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( out ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iTSC   !-- iTimeStepCandidate
    real ( KDR ) :: &
      RampFactor
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO

    I % dT_Candidate  =  huge ( 0.0_KDR )

    associate ( A  =>  I % X )
    do iC  =  1,  A % nCharts

      call I % Compute_dT_Local ( I % dT_Candidate, iC, T_Option )

      select type ( C  =>  A % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      allocate ( CO )
      call CO % Initialize &
             ( C % Communicator, nOutgoing = [ I % n_dT_Candidates ], &
               nIncoming = [ I % n_dT_Candidates ] )

      CO % Outgoing % Value  =  I % dT_Candidate
      call CO % Reduce ( REDUCTION % MIN )
      I % dT_Candidate  =  CO % Incoming % Value

      deallocate ( CO )

      class default 
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Integrator_H_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute_dT', 'subroutine', CONSOLE % ERROR )
      end select !-- C

    end do !-- iC
    end associate !-- A

    do iTSC  =  1,  I % n_dT_Candidates
      call Show ( I % dT_Candidate ( iTSC ), I % Unit_T, &
                  trim ( I % dT_Label ( iTSC ) ) // ' dT', &
                  I % IGNORABILITY + 1 )
    end do !-- iTSC

    dT  =  minval ( I % dT_Candidate )

    RampFactor &
      =  min ( real ( I % iCycle + 1, KDR ) / I % nRampCycles, 1.0_KDR )
    if ( RampFactor  <  1.0_KDR ) then
      dT  =  RampFactor * dT
      call Show ( dT, I % Unit_T, 'Ramped dT', I % IGNORABILITY + 1 )
    end if

  end subroutine Compute_dT


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    I % T_CheckpointInterval &
      =  ( I % T_Finish  -  I % T_Start )  /  I % nWrite

  end subroutine Set_T_CheckpointInterval


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    !-- Very arbitrary for Integrator_H_Form test

    dT_Candidate ( 1 )  =  I % T_CheckpointInterval  /  10
    
  end subroutine Compute_dT_Local


end module Integrator_H__Form
