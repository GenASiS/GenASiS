!-- Integrator_C_1D_PS is a template for time evolution of multiple conserved 
!   currents on position space.

module Integrator_C_1D_PS__Template

  !-- Integrator_Current_1D_PositionSpace__Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    Integrator_C_1D_PS_Template
      integer ( KDI ) :: &
        N_CURRENTS_1D = 0
      type ( Current_ASC_ElementForm ), dimension ( : ), allocatable :: &
        Current_ASC_1D
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_1D_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_1D_PS
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeCycle_ASC
    procedure, private, pass :: &
      ComputeTimeStepLocal
    procedure, public, pass :: &
      ComputeTimeStepLocalTemplate
  end type Integrator_C_1D_PS_Template

contains


  subroutine InitializeTemplate_C_1D_PS &
               ( I, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_1D_PS_Template ), intent ( inout ) :: &
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
      I % Type = 'an Integrator_C_1D_PS'

    if ( I % N_CURRENTS_1D <= 0 ) then
      call Show ( 'I % N_CURRENTS_1D not set to a positive integer', &
                  CONSOLE % WARNING )
      call Show ( I % N_CURRENTS_1D, 'I % N_CURRENTS_1D', CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_PS__Template', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_1D_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Current_ASC_1D ) ) then
      call Show ( 'Current_ASC_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_PS__Template', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_1D_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

  end subroutine InitializeTemplate_C_1D_PS


  impure elemental subroutine FinalizeTemplate_C_1D_PS ( I )

    class ( Integrator_C_1D_PS_Template ), intent ( inout ) :: &
      I

    if ( allocated ( I % Current_ASC_1D ) ) &
      deallocate ( I % Current_ASC_1D )
   
    call I % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_C_1D_PS


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_1D_PS_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % Current_ASC_1D ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, I % N_CURRENTS_1D
      associate ( CA => I % Current_ASC_1D ( iC ) % Element )
      call CA % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption = IgnorabilityOption )
      end associate !-- CA
    end do !-- iC

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine ComputeCycle_ASC ( I, PS )

    class ( Integrator_C_1D_PS_Template ), intent ( inout ) :: &
      I
    class ( Atlas_SC_Form ), intent ( inout ) :: &
      PS

    real ( KDR ) :: &
      TimeNew

    call I % PrepareCycle ( )
    call I % ComputeNewTime ( TimeNew )

    select type ( S => I % Step )
    class is ( Step_RK_C_ASC_1D_Template )

    associate &
      ( CA_1D => I % Current_ASC_1D, &
        TimeStep => TimeNew - I % Time )    

    select type ( Chart => PS % Chart )
    class is ( Chart_SLD_Form )

      call S % Compute ( I % Time, TimeStep )

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_C_1D_PS__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeCycle_ASC', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

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

! if ( I % iCycle > I % nRampCycles &
!      .and. TimeStep < 2.0e-6_KDR * UNIT % SECOND &
!      .and. mod ( I % iCycle, 1000 ) == 0 ) &
! then
!if ( I % iCycle > 73993 ) then
!  I % IsCheckpointTime = .true.
!end if

! if ( TimeStep < 1.0e-12_KDR * UNIT % SECOND ) then
!   call Show ( I % iCycle, '>>> iCycle' )
!   call Show ( TimeStep, I % TimeUnit, '>>> TimeStep too small', CONSOLE % ERROR )
!   call PROGRAM_HEADER % Abort ( )
! end if

    end associate !-- CA_1D, etc.
    end select !-- S

  end subroutine ComputeCycle_ASC


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_1D_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    call I % ComputeTimeStepLocalTemplate ( TimeStepCandidate )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStepLocalTemplate ( I, TimeStepCandidate )

    class ( Integrator_C_1D_PS_Template ), intent ( in ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC  !-- iCurrent
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C

    if ( .not. allocated ( I % Current_ASC_1D ) ) &
      return

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    G => CSL % Geometry ( )

    do iC = 1, I % N_CURRENTS_1D
      associate ( CA => I % Current_ASC_1D ( iC ) % Element )
      C => CA % Current ( )

      call I % ComputeTimeStepKernel_CSL &
             ( CSL % IsProperCell, &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
               C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_LEFT_U ( 2 ) ), & 
               G % Value ( :, G % WIDTH_LEFT_U ( 3 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
               G % Value ( :, G % WIDTH_RIGHT_U ( 2 ) ), & 
               G % Value ( :, G % WIDTH_RIGHT_U ( 3 ) ), &
               G % Value ( :, G % COARSENING ( 2 ) ), &
               G % Value ( :, G % COARSENING ( 3 ) ), &
               CSL % nDimensions, TimeStepCandidate ( iC ) )

      end associate !-- CA
    end do !-- iC

    end select !-- CSL
    end select !-- PS

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( C, G )

  end subroutine ComputeTimeStepLocalTemplate


end module Integrator_C_1D_PS__Template
