!-- Integrator_CS_1D_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets, and an additional conserved
!   current set on position space.

module Integrator_CS_1D_CS__Form

  !-- Integrator_CurrentSet_1D_CurrentSet__Form

  use Basics
  use Integrator_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_Form ) :: Integrator_CS_1D_CS_Form
    integer ( KDI ) :: &
      N_CURRENT_SETS_1D = 0, &
      iCurrentSet = 0
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
  end type Integrator_CS_1D_CS_Form


contains


  subroutine Initialize_H &
               ( I, CommunicatorOption, NameOption, DeviceMemoryOption, &
                 PinnedMemoryOption, DevicesCommunicateOption, &
                 Unit_T_Option, T_FinishOption, nWriteOption )

    class ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
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
      I % Type = 'an Integrator_CS_1D_CS'

    ! if ( .not. allocated ( I % Step_1D ) ) then
    !   call Show ( 'Step_1D not allocated by an extension', &
    !               CONSOLE % WARNING )
    !   call Show ( 'Integrator_C_1D_C_PS__Template', 'module', &
    !               CONSOLE % WARNING )
    !   call Show ( 'InitializeTemplate_C_1D_C_PS', 'subroutine', &
    !               CONSOLE % WARNING )
    ! end if

    ! if ( .not. allocated ( I % TimeSeries ) ) then
    !   allocate ( TimeSeries_C_1D_C_Form :: I % TimeSeries )
    !   !-- Initialized in MS or PS extension of this class
    ! end if

    ! if ( .not. associated ( I % ComputeTimeStepLocal ) ) &
    !   I % ComputeTimeStepLocal => ComputeTimeStepLocal

    call I % Integrator_CS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             Unit_T_Option = Unit_T_Option, & 
             T_FinishOption = T_FinishOption, &
             nWriteOption = nWriteOption )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( I )

    type ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Integrator_CS_1D_CS__Form
