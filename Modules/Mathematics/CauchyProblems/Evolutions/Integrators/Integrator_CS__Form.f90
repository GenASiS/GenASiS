module Integrator_CS__Form

  !-- Integrator_CurrentSet_Form

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Series_CS__Form
  use Integrator_H__Form

  implicit none
  private

  type, public, extends ( Integrator_H_Form ) :: Integrator_CS_Form
    integer ( KDI ) :: &
      iTimer_CT  = 0  !-- ComputeTally
    real ( KDR ) :: &
      CourantFactor
    class ( CurrentSetForm ), allocatable :: &
      CurrentSet_X
    class ( EigenspeedSet_F_Form ), dimension ( : ), allocatable :: &
      EigenspeedSet_X
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    procedure, public, pass :: &   !-- 2
      ShowParameters
    procedure, public, pass :: &   !-- 2
      ShowFields
    final :: &
      Finalize
    procedure, private, pass :: &   !-- 2
      PrepareEvolution
    procedure, public, pass :: &   !-- 3
      UpdateHost => UpdateHost_CS
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, public, nopass :: &   !-- 3
      Analyze_CS
    procedure, public, pass :: &
      Compute_dT_CS_CGS
  end type Integrator_CS_Form

    private :: &
      Compute_dT_Local, &
      InitializeSeries

      private :: &
        Compute_dT_CS_CGS_Kernel

    interface
    
      module subroutine Compute_dT_CS_CGS_Kernel &
               ( dT, ProperCell, &
                 FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, dX_1, dX_2, dX_3, &
                 nDimensions, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          FEP_1, FEP_2, FEP_3, &
          FEM_1, FEM_2, FEM_3, &
          dX_1, dX_2, dX_3
        integer ( KDI ), intent ( in ) :: &
          nDimensions
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_CS_CGS_Kernel

    end interface


contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, &
                 Unit_T_Option, T_FinishOption, nWriteOption )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    type ( QuantityForm ), intent ( in ), optional :: &
      Unit_T_Option
    real ( KDR ), intent ( in ), optional :: &
      T_FinishOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iD
    logical ( KDL ) :: &
      InitializeStep
    character ( 1 ) :: &
      Dimension

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_CS'

    InitializeStep  =  .false.
    if ( .not. allocated ( I % Step_X ) ) then
      allocate ( Step_RK_CS_Form :: I % Step_X )
      InitializeStep  =  .true.
    end if

    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'FastEigenspeed'
    end if

    call I % Integrator_H_Form % Initialize &
           ( CommunicatorOption, NameOption, DeviceMemoryOption, &
             PinnedMemoryOption, DevicesCommunicateOption, &
             Unit_T_Option, T_FinishOption, nWriteOption )

    !-- Integrator methods

    I % Compute_dT_Local  =>  Compute_dT_Local
    I % InitializeSeries  =>  InitializeSeries
    I % Analyze           =>  Analyze_CS

    !-- CurrentSet, if necessary

    if ( .not. allocated ( I % CurrentSet_X ) ) then
      allocate ( I % CurrentSet_X )
      associate &
        ( CS  =>  I % CurrentSet_X, &
           G  =>  I % Geometry_X )
      call CS % Initialize ( G )
      end associate !-- CS, etc.
    end if

    !-- Step, if necessary

    if ( InitializeStep ) then
      select type ( S  =>  I % Step_X )
        class is ( Step_RK_CS_Form )
      associate &
        ( CS  =>  I % CurrentSet_X )
      call S % Initialize ( CS )
      end associate !-- CS, etc.
      end select !-- S
    end if

    !-- EigenspeedSet

    allocate ( I % EigenspeedSet_X ( 3 ) )
    do iD  =  1, 3
      associate &
        ( ES  =>  I % EigenspeedSet_X ( iD ), &
          CS  =>  I % CurrentSet_X )
      write ( Dimension, fmt = '(i1.1)' ) iD
      call ES % Initialize ( CS, CS, SuffixOption = Dimension ) 
      end associate !-- ES, etc.
    end do !-- iD

    !-- Stream

    select type ( S  =>  I % Step_X )
      class is ( Step_RK_CS_Form )
    associate &
      ( CS_X  =>  I % CurrentSet_X, &
         S_X  =>  I % Checkpoint_X )
    call CS_X % SetStream ( S_X )
    call  S   % SetStream ( S_X )
    end associate !-- CS_X, etc.
    end select !-- S

    !-- Courant factor

    I % CourantFactor  =  0.95_KDR
    call PROGRAM_HEADER % GetParameter ( I % CourantFactor, 'CourantFactor' )

  end subroutine Initialize_H


  subroutine ShowParameters ( I )

    class ( Integrator_CS_Form ), intent ( in ) :: &
      I

    call I % Integrator_H_Form % ShowParameters ( )

    call Show ( I % CourantFactor, 'CourantFactor', I % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowFields ( I )

    class ( Integrator_CS_Form ), intent ( in ) :: &
      I

    integer ( KDI ) :: &
      iD

    call I % Integrator_H_Form % ShowFields ( )

    call I % CurrentSet_X % Show ( )

    do iD  =  1, 3
      call I % EigenspeedSet_X ( iD ) % Show ( )
    end do !-- iD

  end subroutine ShowFields


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % EigenspeedSet_X ) ) &
      deallocate ( I % EigenspeedSet_X )
    if ( allocated ( I % CurrentSet_X ) ) &
      deallocate ( I % CurrentSet_X )

  end subroutine Finalize


  subroutine PrepareEvolution ( I )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I

    associate ( G  =>  I % Geometry_X )
    call G % UpdateDevice ( )
    end associate !-- G

    if ( .not. allocated ( I % CurrentSet_X ) ) &
      return

    associate ( CS  =>  I % CurrentSet_X )
    call CS % UpdateDevice ( )
    call CS % ExchangeGhostData ( )
    call CS % ComputeFromInitial ( )
    call CS % ApplyBoundaryConditions ( )
    call CS % UpdateHost ( )
    end associate !-- CS

  end subroutine PrepareEvolution


  subroutine UpdateHost_CS ( I )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I

    call I % Integrator_H_Form % UpdateHost ( )

    associate ( CS  =>  I % CurrentSet_X )
    call CS % UpdateHost ( )
    end associate !-- CS

  end subroutine UpdateHost_CS


  subroutine ComputeTally ( I, ChangeOption, IgnorabilityOption )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( .not. allocated ( I % CurrentSet_X ) ) &
      return

    associate ( CS => I % CurrentSet_X )
    call CS % ComputeTally &
           ( ChangeOption = ChangeOption, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CS

  end subroutine ComputeTally


  subroutine Analyze_CS ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_CT

    select type ( I )
      class is ( Integrator_CS_Form )

    !-- Tally

    if ( present ( T_Option ) ) then
      T_CT  =>  PROGRAM_HEADER % Timer &
                  ( Handle = I % iTimer_CT, &
                    Name = trim ( I % Name ) // '_CmptTlly', &
                    Level = T_Option % Level + 1 )
    else
      T_CT  =>  null ( )
    end if
    if ( associated ( T_CT ) ) call T_CT % Start ( )   
    call I % ComputeTally &
           ( ChangeOption = .not. I % Start .and. .not. I % Restart, &
             IgnorabilityOption  = Ignorability )
    if ( associated ( T_CT ) ) call T_CT % Stop ( )   

    !-- Analyze_H

    call I % Analyze_H ( Ignorability, T_Option )

    end select !-- I

  end subroutine Analyze_CS


  subroutine Compute_dT_CS_CGS ( I, dT, iC, T_Option )

    class ( Integrator_CS_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( inout ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    associate &
      ( ES_1  =>  I % EigenspeedSet_X ( 1 ), &
        ES_2  =>  I % EigenspeedSet_X ( 2 ), &
        ES_3  =>  I % EigenspeedSet_X ( 3 ), &
         G    =>  I % Geometry_X )

    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    call ES_1 % Compute ( iC = 1, iD = 1 )
    call ES_2 % Compute ( iC = 1, iD = 2 )
    call ES_3 % Compute ( iC = 1, iD = 3 )

    associate &
      ( EV_1  =>  ES_1 % Storage ( 1 ) % Value, &
        EV_2  =>  ES_2 % Storage ( 1 ) % Value, &
        EV_3  =>  ES_3 % Storage ( 1 ) % Value, &
        GV    =>   G   % Storage ( 1 ) % Value )

    call Compute_dT_CS_CGS_Kernel &
           ( dT, C % ProperCell, &
             EV_1 ( :, ES_1 % EIGENSPEED_FAST_PLUS_U ), &
             EV_2 ( :, ES_2 % EIGENSPEED_FAST_PLUS_U ), &
             EV_3 ( :, ES_3 % EIGENSPEED_FAST_PLUS_U ), &
             EV_1 ( :, ES_1 % EIGENSPEED_FAST_MINUS_U ), &
             EV_2 ( :, ES_2 % EIGENSPEED_FAST_MINUS_U ), &
             EV_3 ( :, ES_3 % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, G % WIDTH_U_1 ), &
             GV ( :, G % WIDTH_U_2 ), &
             GV ( :, G % WIDTH_U_3 ), &
             C % nDimensions, &
             UseDeviceOption = G % DeviceMemory )

    end associate !-- EV, etc.
    end associate !-- C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Integrator_CS_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_CS_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end associate !-- ES_1, etc.

    dT  =  I % CourantFactor  *  dT
    
  end subroutine Compute_dT_CS_CGS


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( I )
    class is ( Integrator_CS_Form )
      call I % Compute_dT_CS_CGS ( dT_Candidate ( 1 ), iC, T_Option )
    end select !-- I

  end subroutine Compute_dT_Local


  subroutine InitializeSeries ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    allocate ( Series_CS_Form :: I % Series )

    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( S  =>  I % Series )
      class is ( Series_CS_Form )
    call S % Initialize &
      ( I % CurrentSet_X, I % GridImageStream, I % dT_Label, I % Unit_T, &
        I % dT_Candidate, I % T, I % Communicator % Rank, I % nWrite, &
        I % iCycle )
    end select !-- S
    end select !-- I

  end subroutine InitializeSeries


end module Integrator_CS__Form
