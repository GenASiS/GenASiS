!-- Integrator_C_PS_C_PS is a template for operator-split time evolution of two
!   conserved currents on position space.

module Integrator_C_PS_C_PS__Template

  !-- Integrator_Current_PositionSpace_Current_PositionSpace__Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    Integrator_C_PS_C_PS_Template
      class ( Current_ASC_Template ), allocatable :: &
        Current_ASC_1, &
        Current_ASC_2
      class ( Step_RK_C_ASC_Template ), allocatable :: &
        Step_1, &
        Step_2
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_PS_C_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_PS_C_PS
    procedure, private, pass :: &  !-- 3
      InitializeStepTimers
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeCycle_ASC
    procedure, private, pass :: &
      PrepareStep_1
    procedure, private, pass :: &
      PrepareStep_2
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type Integrator_C_PS_C_PS_Template

contains


  subroutine InitializeTemplate_C_PS_C_PS &
               ( I, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_C_PS_C_PS'

    if ( .not. allocated ( I % Current_ASC_1 ) ) then
      call Show ( 'Current_ASC_1 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Current_ASC_2 ) ) then
      call Show ( 'Current_ASC_2 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Step_1 ) ) then
      call Show ( 'Step_1 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_C_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Step_2 ) ) then
      call Show ( 'Step_2 not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS_C_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine InitializeTemplate_C_PS_C_PS



  impure elemental subroutine FinalizeTemplate_C_PS_C_PS ( I )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % TimeSeries ) ) &
     deallocate ( I % TimeSeries )
   if ( allocated ( I % Step_2 ) ) &
     deallocate ( I % Step_2 )
   if ( allocated ( I % Step_1 ) ) &
     deallocate ( I % Step_1 )
   if ( allocated ( I % Current_ASC_2 ) ) &
     deallocate ( I % Current_ASC_2 )
   if ( allocated ( I % Current_ASC_1 ) ) &
     deallocate ( I % Current_ASC_1 )

    call I % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C_PS_C_PS


  subroutine InitializeStepTimers ( I, BaseLevel )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( allocated ( I % Step_1 ) ) &
      call I % Step_1 % InitializeTimers ( BaseLevel )
    if ( allocated ( I % Step_2 ) ) &
      call I % Step_2 % InitializeTimers ( BaseLevel )

  end subroutine InitializeStepTimers


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % Current_ASC_1 ) ) &
      return
    if ( .not. allocated ( I % Current_ASC_2 ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    associate ( CA => I % Current_ASC_1 )
    call CA % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CA

    associate ( CA => I % Current_ASC_2 )
    call CA % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CA

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine ComputeCycle_ASC ( I, PS )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I
    class ( Atlas_SC_Form ), intent ( inout ) :: &
      PS

    real ( KDR ) :: &
      TimeNew

    select type ( Chart => PS % Chart )
    class is ( Chart_SLD_Form )

    call I % ComputeNewTime ( TimeNew )
    associate ( TimeStep => TimeNew - I % Time )    

    associate &
      ( S_1 => I % Step_1, &
        S_2 => I % Step_2 )

    call I % PrepareStep_1 ( )
    call S_1 % Compute ( I % Time, TimeStep )
    
    call I % PrepareStep_2 ( )
    call S_2 % Compute ( I % Time, TimeStep )

    end associate !-- S_1, S_2

    I % iCycle = I % iCycle + 1
    I % Time = I % Time + TimeStep

    if ( I % WriteTimeExact ) then
      if ( I % Time == I % WriteTime ) &
        I % IsCheckpointTime = .true.
    else 
      if ( I % Time  >  I % WriteTime &
           .or. abs ( I % Time - I % WriteTime )  <  0.5_KDR * TimeStep ) &
        I % IsCheckpointTime = .true.
    end if

    end associate !-- TimeStep

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_C_PS_C_PS__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeCycle_ASC', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeCycle_ASC


  subroutine PrepareStep_1 ( I )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep_1


  subroutine PrepareStep_2 ( I )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep_2


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_PS_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    call I % ComputeTimeStep_C_ASC &
           ( TimeStepCandidate ( 1 ), I % Current_ASC_1 )
    call I % ComputeTimeStep_C_ASC &
           ( TimeStepCandidate ( 2 ), I % Current_ASC_2 )

  end subroutine ComputeTimeStepLocal


end module Integrator_C_PS_C_PS__Template
