!-- Integrator_CS_1D_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets, and an additional conserved
!   current set on position space.

module Integrator_CS_1D_CS__Form

  !-- Integrator_CurrentSet_1D_CurrentSet__Form

  use Basics
  use Steps
  use Integrator_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_Form ) :: Integrator_CS_1D_CS_Form
    integer ( KDI ) :: &
      N_CURRENT_SETS_1D = 0, &
      iCurrentSet = 0
    class ( Step_RK_H_Form ), allocatable :: &
      Step_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
    procedure, private, pass :: &  !-- 2
      ComputeCycle
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

    if ( .not. allocated ( I % Step_1D ) ) then
      call Show ( 'Step_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_CS_1D_CS__Form', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'Initialize_H', 'subroutine', &
                  CONSOLE % WARNING )
    end if

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

    ! nullify ( I % PrepareStep )
    ! nullify ( I % PrepareStep_1D )
    ! nullify ( I % SeriesChangeGrandTotal )
    ! nullify ( I % iTime )

    if ( allocated ( I % Step_1D ) ) &
      deallocate ( I % Step_1D )

  end subroutine Finalize


  subroutine ComputeCycle ( I, T_CC )

    class ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_CC

    real ( KDR ) :: &
      T_New  !-- Use of Compute_T_New is a relic of past AMR evolution
    type ( TimerForm ), pointer :: &
      T_CTN, &
      T_S

    T_CTN  =>  PROGRAM_HEADER % Timer &
                ( Handle = I % iTimer_CTN, &
                  Name = trim ( I % Name ) // '_CmptDt', &
                  Level = T_CC % Level + 1 )
    call T_CTN % Start ( )
    call I % Compute_T_New ( T_New )
    call T_CTN % Stop ( )

    associate ( dT  =>  T_New  -  I % T )    

    ! select type ( Chart => PS % Chart )
    ! class is ( Chart_SLD_Form )

    if ( allocated ( I % Step_X ) ) then
      associate ( S  =>  I % Step_X )
      T_S  =>  S % Timer ( Level = T_CC % Level + 1 )
      call T_S % Start ( )
      call S % Compute ( I % T, dT, T_Option = T_S )
      call T_S % Stop ( )
      end associate !--  S
    end if !-- allocated Step

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

  end subroutine ComputeCycle


end module Integrator_CS_1D_CS__Form
