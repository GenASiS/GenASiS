module Integrator_CSA__Form

  !-- Integrator_CurrentSetAtlas_Form

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_H__Form

  implicit none
  private

  type, public, extends ( Integrator_H_Form ) :: Integrator_CSA_Form
    real ( KDR ) :: &
      CourantFactor
    class ( CurrentSet_A_Form ), allocatable :: &
      CurrentSet_X_A
    class ( FieldSet_A_Element ), dimension ( : ), allocatable :: &
      Eigenspeeds_X_A
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    procedure, public, pass :: &  !-- 2
      ShowParameters
    procedure, public, pass :: &  !-- 2
      ShowFields
    final :: &
      Finalize
    procedure, private, pass :: &   !-- 2
      PrepareEvolution
    procedure, public, pass :: &   !-- 3
      UpdateHost => UpdateHost_CSA
    procedure, public, pass :: &
      Compute_dT_CSC
  end type Integrator_CSA_Form

    private :: &
      Compute_dT_Local

    private :: &
      Compute_dT_CSC_Kernel

  interface
  
    module subroutine Compute_dT_CSC_Kernel &
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
    end subroutine Compute_dT_CSC_Kernel

  end interface

contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, &
                 Unit_T_Option, T_FinishOption, nWriteOption )

    class ( Integrator_CSA_Form ), intent ( inout ) :: &
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

    integer ( KDI ) :: &
      iD
    logical ( KDL ) :: &
      InitializeStep
    character ( 1 ) :: &
      Dimension

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_CSA'

    InitializeStep  =  .false.
    if ( .not. allocated ( I % Step_X_A ) ) then
      allocate ( Step_RK_CSA_Form :: I % Step_X_A )
      InitializeStep  =  .true.
    end if

    if ( .not. allocated ( I % dT_Label ) ) then
      allocate ( I % dT_Label ( 1 ) )
      I % dT_Label ( 1 ) = 'FastEigenspeed'
    end if

    I % Compute_dT_Local  =>  Compute_dT_Local

    call I % Integrator_H_Form % Initialize &
           ( CommunicatorOption, NameOption, DeviceMemoryOption, &
             PinnedMemoryOption, DevicesCommunicateOption, &
             Unit_T_Option, T_FinishOption, nWriteOption )

    !-- CurrentSet, if necessary

    if ( .not. allocated ( I % CurrentSet_X_A ) ) then
      allocate ( I % CurrentSet_X_A )
      associate &
        ( CSA  =>  I % CurrentSet_X_A, &
           GA  =>  I % Geometry_X_A, &
           SA  =>  I % Checkpoint_X_A )
      call CSA % Initialize( GA )
      call CSA % SetStream ( SA )
      end associate !-- CSA, etc.
    end if

    !-- Step, if necessary

    if ( InitializeStep ) then
      select type ( S  =>  I % Step_X_A )
        class is ( Step_RK_CSA_Form )
      associate &
        ( CSA  =>  I % CurrentSet_X_A, &
           SA  =>  I % Checkpoint_X_A )
      call S % Initialize ( CSA )
      call S % SetStream ( SA )
      end associate !-- CSA, etc.
      end select !-- S
    end if

    !-- Eigenspeeds

    allocate ( I % Eigenspeeds_X_A ( 3 ) )
    do iD  =  1, 3
      associate &
        ( EAE  =>  I % Eigenspeeds_X_A ( iD ), &
          CSA  =>  I % CurrentSet_X_A )
      allocate ( Eigenspeeds_F_A_Form :: EAE % Element )
      select type ( EA  =>  EAE % Element )
        class is ( Eigenspeeds_F_A_Form )
      write ( Dimension, fmt = '(i1.1)' ) iD
      call EA % Initialize &
             ( CSA, &
               NameOption = 'E_' // Dimension // '_' // trim ( CSA % Name ) ) 
      end select !-- EA
      end associate !-- EAE, CSA
    end do !-- iD

    !-- Courant factor

    I % CourantFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter ( I % CourantFactor, 'CourantFactor' )

  end subroutine Initialize_H


  subroutine ShowParameters ( I )

    class ( Integrator_CSA_Form ), intent ( in ) :: &
      I

    call I % Integrator_H_Form % ShowParameters ( )

    call Show ( I % CourantFactor, 'CourantFactor', I % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowFields ( I )

    class ( Integrator_CSA_Form ), intent ( in ) :: &
      I

    integer ( KDI ) :: &
      iD

    call I % Integrator_H_Form % ShowFields ( )

    call I % CurrentSet_X_A % Show ( )

    do iD  =  1, 3
      associate ( EA  =>  I % Eigenspeeds_X_A ( iD ) % Element )
      call EA % Show ( )
      end associate !-- EA
    end do !-- iD

  end subroutine ShowFields


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CSA_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % Eigenspeeds_X_A ) ) &
      deallocate ( I % Eigenspeeds_X_A )
    if ( allocated ( I % CurrentSet_X_A ) ) &
      deallocate ( I % CurrentSet_X_A )

  end subroutine Finalize


  subroutine PrepareEvolution ( I )

    class ( Integrator_CSA_Form ), intent ( inout ) :: &
      I

    associate ( GA  =>  I % Geometry_X_A )
    call GA % UpdateDevice ( )
    end associate !-- GA

    if ( .not. allocated ( I % CurrentSet_X_A ) ) &
      return

    associate ( CSA  =>  I % CurrentSet_X_A )
    call CSA % UpdateDevice ( )
    call CSA % ExchangeGhostData ( )
    call CSA % ComputeFromInitial ( )
    call CSA % UpdateHost ( )
    end associate !-- CSA

    ! call I % ComputeConstraints ( )

  end subroutine PrepareEvolution


  subroutine UpdateHost_CSA ( I, TimerLevelOption )

    class ( Integrator_CSA_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call I % Integrator_H_Form % UpdateHost ( TimerLevelOption )

    associate ( CSA  =>  I % CurrentSet_X_A )
    call CSA % UpdateHost ( TimerLevelOption )
    end associate !-- CSA

  end subroutine UpdateHost_CSA


  subroutine Compute_dT_CSC ( I, dT, iC, TimerLevelOption )

    class ( Integrator_CSA_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( inout ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    associate &
      ( EAE  =>  I % Eigenspeeds_X_A, &
        CSA  =>  I % CurrentSet_X_A )
    select type ( EC_1  =>  EAE ( 1 ) % Element % FieldSet_C ( iC ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( EC_2  =>  EAE ( 2 ) % Element % FieldSet_C ( iC ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( EC_3  =>  EAE ( 3 ) % Element % FieldSet_C ( iC ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( GC  =>  CSA % Geometry_A % FieldSet_C ( iC ) % Element )
      class is ( Geometry_F_C_Form )
    select type ( C  =>  GC % Chart )
      class is ( Chart_GS_Form )

    call EC_1 % Compute ( iD = 1, TimerLevelOption = TimerLevelOption )
    call EC_2 % Compute ( iD = 2, TimerLevelOption = TimerLevelOption )
    call EC_3 % Compute ( iD = 3, TimerLevelOption = TimerLevelOption )

    associate &
      ( EV_1  =>  EC_1 % Storage_FSC % Storage % Value, &
        EV_2  =>  EC_2 % Storage_FSC % Storage % Value, &
        EV_3  =>  EC_3 % Storage_FSC % Storage % Value, &
        GV    =>  GC   % Storage_FSC % Storage % Value, &
        DeviceMemory  =>  GC % Storage_FSC % DeviceMemory )

    call Compute_dT_CSC_Kernel &
           ( dT, C % ProperCell, &
             EV_1 ( :, EC_1 % EIGENSPEED_FAST_PLUS_U ), &
             EV_2 ( :, EC_2 % EIGENSPEED_FAST_PLUS_U ), &
             EV_3 ( :, EC_3 % EIGENSPEED_FAST_PLUS_U ), &
             EV_1 ( :, EC_1 % EIGENSPEED_FAST_MINUS_U ), &
             EV_2 ( :, EC_2 % EIGENSPEED_FAST_MINUS_U ), &
             EV_3 ( :, EC_3 % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, GC % WIDTH_U_1 ), &
             GV ( :, GC % WIDTH_U_2 ), &
             GV ( :, GC % WIDTH_U_3 ), &
             C % nDimensions, &
             UseDeviceOption = DeviceMemory )

    end associate !-- EV, etc.
    end select !-- C
    end select !-- GC
    end select !-- EC_3
    end select !-- EC_2
    end select !-- EC_1
    end associate !-- EAE, etc.

    dT  =  I % CourantFactor  *  dT
    
  end subroutine Compute_dT_CSC


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, TimerLevelOption )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    select type ( I )
      class is ( Integrator_CSA_Form )

    call I % Compute_dT_CSC ( dT_Candidate ( 1 ), iC, TimerLevelOption )

    end select !-- I

  end subroutine Compute_dT_Local


end module Integrator_CSA__Form
