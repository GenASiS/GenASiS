!-- Integrator_CS_1D_CS_PS is a parent for time evolution of multiple
!   similar conserved current sets, and an additional conserved
!   current set on position space.

module Integrator_CS_1D_CS__Form

  !-- Integrator_CurrentSet_1D_CurrentSet__Form

  use Basics
  use Fields
  use Steps
  use Integrator_CS__Form

  implicit none
  private

  type, public, extends ( Integrator_CS_Form ) :: Integrator_CS_1D_CS_Form
    integer ( KDI ) :: &
      nCurrentSets = 0, &
      iCurrentSet = 0
    real ( KDR ) :: &
      CourantFactor_1D
    class ( Step_RK_H_Form ), allocatable :: &
      Step_1D
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_H      
    final :: &
      Finalize
    procedure, public, pass :: &  !-- 2
      ShowParameters
    procedure, public, pass :: &  !-- 2
      ShowSteps
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeTally_1D
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
      call Show ( 'Step_1D not allocated', CONSOLE % WARNING )
      call Show ( 'Integrator_CS_1D_CS__Form', 'module', CONSOLE % WARNING )
      call Show ( 'Initialize_H', 'subroutine', CONSOLE % WARNING )
    end if

    ! if ( .not. allocated ( I % TimeSeries ) ) then
    !   allocate ( TimeSeries_C_1D_C_Form :: I % TimeSeries )
    !   !-- Initialized in MS or PS extension of this class
    ! end if

    call I % Integrator_CS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             Unit_T_Option = Unit_T_Option, & 
             T_FinishOption = T_FinishOption, &
             nWriteOption = nWriteOption )

    !-- Courant factor

    I % CourantFactor_1D  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter &
           ( I % CourantFactor_1D, 'CourantFactor_1D' )

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


  subroutine ShowParameters ( I )

    class ( Integrator_CS_1D_CS_Form ), intent ( in ) :: &
      I

    call I % Integrator_CS_Form % ShowParameters ( )

    call Show ( I % nCurrentSets, 'nCurrentSets', I % IGNORABILITY )
    call Show ( I % iCurrentSet,  'iCurrentSet',  I % IGNORABILITY )
    call Show ( I % CourantFactor_1D, 'CourantFactor_1D', I % IGNORABILITY )


  end subroutine ShowParameters


  subroutine ShowSteps ( I )

    class ( Integrator_CS_1D_CS_Form ), intent ( in ) :: &
      I

    call I % Integrator_CS_Form % ShowSteps ( )

    if ( allocated ( I % Step_1D ) ) &
      call I % Step_1D % Show ( )

  end subroutine ShowSteps


  subroutine ComputeCycle ( I, T_CC )

    class ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I
    type ( TimerForm ), intent ( in ) :: &
      T_CC

    real ( KDR ) :: &
      T_New  !-- Use of Compute_T_New is a relic of past AMR evolution
    type ( TimerForm ), pointer :: &
      T_CTN, &
      T_S, &
      T_S_1D

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

    if ( allocated ( I % Step_1D ) ) then
      associate ( S_1D  =>  I % Step_1D )
      T_S_1D  =>  S_1D % Timer ( Level = T_CC % Level + 1 )
      call T_S_1D % Start ( )
      call S_1D % Compute ( I % T, dT, T_Option = T_S_1D )
      call T_S_1D % Stop ( )
      end associate !-- S_1D
    end if !-- allocated Step_1D

    if ( allocated ( I % Step_X ) ) then
      associate ( S  =>  I % Step_X )
      T_S  =>  S % Timer ( Level = T_CC % Level + 1 )
      call T_S % Start ( )
      call S % Compute ( I % T, dT, T_Option = T_S )
      call T_S % Stop ( )
      end associate !-- S
    end if !-- allocated Step_X

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


  subroutine ComputeTally ( I, ChangeOption, IgnorabilityOption )

    class ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call I % ComputeTally_1D &
           ( ChangeOption = ChangeOption, &
             IgnorabilityOption  = IgnorabilityOption )

    if ( allocated ( I % CurrentSet_X ) ) then
      associate ( CS => I % CurrentSet_X )
      call CS % ComputeTally &
             ( ChangeOption = ChangeOption, &
               IgnorabilityOption = IgnorabilityOption )
      end associate !-- CS
    end if

    ! if ( associated ( I % SeriesChangeGrandTotal ) ) then
    !   associate ( SCGT => I % SeriesChangeGrandTotal )
    !   call Show ( 'Change in Grand Total Tally', IgnorabilityOption )
    !   do iV = 1, SCGT % nVariables
    !     iS = SCGT % iaSelected ( iV )
    !     call Show ( SCGT % Value ( I % iTime, iS ), SCGT % Unit ( iS ), &
    !                 SCGT % Variable ( iS ), IgnorabilityOption )
    !   end do !-- iV
    !   end associate !-- SCGT, etc.
    ! end if

  end subroutine ComputeTally


  subroutine ComputeTally_1D ( I, ChangeOption, IgnorabilityOption )

    class ( Integrator_CS_1D_CS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    !-- To be filled in by extension

  end subroutine ComputeTally_1D


end module Integrator_CS_1D_CS__Form
