!-- Integrator_CS_1D_BM_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets on the base manifold, and an additional 
!   conserved current set on the base manifold.

module Integrator_CS_1D_BM_CS__Form

  !-- Integrator_CurrentSet_1D_BaseManifold_CurrentSet__Form

  use Basics
  use Fields
  use Steps
  use Integrator_H__Form
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_BM_CS_Form
      type ( CommunicatorForm ), allocatable :: &
        Communicator_X_1D
      class ( CurrentSetForm ), allocatable :: &
        CurrentSet_X_1D
      class ( EigenspeedSet_F_Form ), dimension ( : ), allocatable :: &
        EigenspeedSet_X_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
    procedure, public, pass :: &  !-- 2
      ShowParameters
    procedure, public, pass :: &   !-- 2
      ShowFields
    procedure, public, pass :: &   !-- 2
      PrepareEvolution
    procedure, private, pass :: &   !-- 2
      SetCommunicator_1D
    procedure, public, pass :: &   !-- 3
      UpdateHost => UpdateHost_CS_1D
    procedure, private, pass :: &
      ComputeTally_1D
  end type Integrator_CS_1D_BM_CS_Form


contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, &
                 Unit_T_Option, T_FinishOption, nWriteOption )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
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
    character ( 1 ) :: &
      Suffix

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_CS_1D_BM_CS'

    call I % Integrator_CS_1D_CS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             Unit_T_Option = Unit_T_Option, & 
             T_FinishOption = T_FinishOption, &
             nWriteOption = nWriteOption )

    !-- CurrentSet, if necessary

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) then
      allocate ( I % CurrentSet_X_1D )
      associate &
        ( CS  =>  I % CurrentSet_X_1D, &
           G  =>  I % Geometry_X )
      write ( Suffix, fmt = '(i1.1)' ) I % iCurrentSet
      call CS % Initialize ( G, NameOption = 'CurrentSet_X_1D_' // Suffix )
      end associate !-- CS, etc.
    end if

    !-- EigenspeedSet

    allocate ( I % EigenspeedSet_X_1D ( 3 ) )
    do iD  =  1, 3
      associate &
        ( ES  =>  I % EigenspeedSet_X_1D ( iD ), &
          CS  =>  I % CurrentSet_X_1D )
      write ( Suffix, fmt = '(i1.1)' ) iD
      call ES % Initialize ( CS, CS, SuffixOption = Suffix ) 
      end associate !-- ES, etc.
    end do !-- iD

    !-- Stream

    associate &
      ( CS_X_1D  =>  I % CurrentSet_X_1D, &
         S_X     =>  I % Checkpoint_X )

    call CS_X_1D % SetStream ( S_X )

    if ( allocated ( I % Step_1D ) ) then
      associate ( S  =>  I % Step_1D )
      call S % SetStream ( S_X )
      end associate !-- S
    end if

    end associate !-- CS_X, etc.

    !-- Communicator

    call I % SetCommunicator_1D ( )

    !-- Series

    ! select type ( TS  =>  I % TimeSeries )
    ! type is ( TimeSeries_C_1D_C_Form )
    !   I % InitializeTimeSeries  =>  InitializeTimeSeries_C_1D_PS
    ! end select !-- TS
    
  end subroutine Initialize_H


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

    if ( allocated ( I % EigenspeedSet_X_1D ) ) &
      deallocate ( I % EigenspeedSet_X_1D )
    if ( allocated ( I % CurrentSet_X_1D ) ) &
      deallocate ( I % CurrentSet_X_1D )
    if ( allocated ( I % Communicator_X_1D ) ) &
      deallocate ( I % Communicator_X_1D )

  end subroutine Finalize


  subroutine ShowParameters ( I )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( in ) :: &
      I

    call I % Integrator_CS_1D_CS_Form % ShowParameters ( )

    call I % Communicator_X_1D % Show ( I % IGNORABILITY )

  end subroutine ShowParameters


  subroutine ShowFields ( I )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( in ) :: &
      I

    integer ( KDI ) :: &
      iD

    call I % Integrator_CS_Form % ShowFields ( )

    call I % CurrentSet_X_1D % Show ( )

    do iD  =  1, 3
      call I % EigenspeedSet_X_1D ( iD ) % Show ( )
    end do !-- iD

  end subroutine ShowFields


  subroutine PrepareEvolution ( I )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

    call I % Integrator_CS_Form % PrepareEvolution ( )

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) &
      return

    associate ( CS  =>  I % CurrentSet_X_1D )
    call CS % UpdateDevice ( )
    call CS % ExchangeGhostData ( )
    call CS % ComputeFromInitial ( )
    call CS % ApplyBoundaryConditions ( )
    call CS % UpdateHost ( )
    end associate !-- CS

  end subroutine PrepareEvolution
  
  
  subroutine SetCommunicator_1D ( I )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iR  !-- iRadiation
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank

    allocate ( I % Communicator_X_1D )

    associate &
      ( C     =>  I % Communicator, &
        C_P   =>  I % Communicator % Parent, &
        C_1D  =>  I % Communicator_X_1D )

    allocate ( Rank, &
               source = [ ( iR, iR = C % Rank, C_P % Size - 1, C % Size ) ] )

    !-- C_1D % Size will be I % N_CURRENT_SETS_1D

    call C_1D % Initialize ( C_P, Rank, NameOption = 'Communicator_X_1D' )

    end associate !-- C, etc.

  end subroutine SetCommunicator_1D


  subroutine UpdateHost_CS_1D ( I )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

    call I % Integrator_CS_Form % UpdateHost ( )

    associate ( CS  =>  I % CurrentSet_X_1D )
    call CS % UpdateHost ( )
    end associate !-- CS

  end subroutine UpdateHost_CS_1D


  subroutine ComputeTally_1D ( I, ChangeOption, IgnorabilityOption )

    class ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( .not. allocated ( I % CurrentSet_X_1D ) ) &
      return

    associate ( CS => I % CurrentSet_X_1D )
    call CS % ComputeTally &
           ( ChangeOption = ChangeOption, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CS

  end subroutine ComputeTally_1D


end module Integrator_CS_1D_BM_CS__Form
