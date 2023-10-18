!-- Integrator_CS_1D_BM_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets on the base manifold, and an additional 
!   conserved current set on the base manifold.

module Integrator_CS_1D_BM_CS__Form

  !-- Integrator_CurrentSet_1D_BaseManifold_CurrentSet__Form

  use Basics
  use Fields
  use Steps
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_BM_CS_Form
      class ( CurrentSetForm ), allocatable :: &
        CurrentSet_X_1D
      class ( EigenspeedSet_F_Form ), dimension ( : ), allocatable :: &
        EigenspeedSet_X_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
    procedure, public, pass :: &   !-- 2
      ShowFields
    procedure, public, pass :: &   !-- 2
      PrepareEvolution
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

    !-- CurrentSet, if necessary. Member iCurrentSet should already be set.

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

    select type ( S  =>  I % Step_1D )
      class is ( Step_RK_CS_Form )
    associate &
      ( CS_X_1D  =>  I % CurrentSet_X_1D, &
         S_X     =>  I % Checkpoint_X )
    call CS_X_1D % SetStream ( S_X )
    call  S      % SetStream ( S_X )
    end associate !-- CS_X_1D, etc.
    end select !-- S

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

  end subroutine Finalize


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


end module Integrator_CS_1D_BM_CS__Form
