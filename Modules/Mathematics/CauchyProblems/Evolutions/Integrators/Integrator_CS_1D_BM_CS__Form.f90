!-- Integrator_CS_1D_BM_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets on the base manifold, and an additional 
!   conserved current set on the base manifold.

module Integrator_CS_1D_BM_CS__Form

  !-- Integrator_CurrentSet_1D_BaseManifold_CurrentSet__Form

  use Basics
  use Fields
  use Integrator_CS_1D_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_1D_CS_Form ) :: &
    Integrator_CS_1D_BM_CS_Form
      class ( CurrentSetForm ), allocatable :: &
        CurrentSet_X_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
    procedure, public, pass :: &   !-- 2
      ShowFields
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

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_CS_1D_BM_CS'

    ! if ( .not. allocated ( I % Current_ASC_1D ) ) then
    !   call Show ( 'Current_ASC_1D not allocated by an extension', &
    !               CONSOLE % WARNING )
    !   call Show ( 'Integrator_C_1D_PS_C_PS__Form', 'module', &
    !               CONSOLE % WARNING )
    !   call Show ( 'Initialize', 'subroutine', &
    !               CONSOLE % WARNING )
    ! end if

    call I % Integrator_CS_1D_CS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             Unit_T_Option = Unit_T_Option, & 
             T_FinishOption = T_FinishOption, &
             nWriteOption = nWriteOption )

    !-- Stream

    ! select type ( S  =>  I % Step_X )
    !   class is ( Step_RK_CS_Form )
    associate &
      ( CS_X_1D  =>  I % CurrentSet_X_1D, &
         S_X     =>  I % Checkpoint_X )
    call CS_X_1D % SetStream ( S_X )
    ! call  S   % SetStream ( S_X )
    end associate !-- CS_X_1D, etc.
    ! end select !-- S

    ! select type ( TS  =>  I % TimeSeries )
    ! type is ( TimeSeries_C_1D_C_Form )
    !   I % InitializeTimeSeries  =>  InitializeTimeSeries_C_1D_PS
    ! end select !-- TS

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_BM_CS_Form ), intent ( inout ) :: &
      I

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

!    do iD  =  1, 3
!      call I % EigenspeedSet_X ( iD ) % Show ( )
!    end do !-- iD

  end subroutine ShowFields


end module Integrator_CS_1D_BM_CS__Form
