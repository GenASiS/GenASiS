!-- Integrator_C_1D_MS_C_PS is a template for time evolution of multiple
!   conserved currents on momentum space (and position space), and a conserved
!   current on position space.

module Integrator_C_1D_MS_C_PS__Template

  !-- Integrator_Current_1D_MomentumSpace_C_PositionSpace__Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_C_PS__Template

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ), abstract :: &
    Integrator_C_1D_MS_C_PS_Template
      integer ( KDI ) :: &
        N_CURRENTS_MS = 0
      logical ( KDL ) :: &
        UseLimiterParameter_S, &
        UseLimiterParameter_F
      type ( Current_BSLL_ASC_CSLD_ElementForm ), dimension ( : ), &
        allocatable :: &
          Current_BSLL_ASC_CSLD_1D
      class ( Step_RK_C_ASC_Template ), allocatable :: &
        Step_PS
      class ( Step_RK_C_BSLL_ASC_CSLD_1D_Template ), allocatable :: &
        Step_MS
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C_1D_MS_C_PS
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C_1D_MS_C_PS
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      InitializeStepTimers
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &
      ComputeCycle_BSLL_ASC_CSLD_1D
    procedure, private, pass :: &
      PrepareStep_MS
    procedure, private, pass :: &
      PrepareStep_PS
    procedure, private, pass :: &
      ComputeTimeStepLocal
    procedure, public, pass :: &
      ComputeTimeStepLocalTemplate
  end type Integrator_C_1D_MS_C_PS_Template

contains


  subroutine InitializeTemplate_C_1D_MS_C_PS &
               ( I, Name, UseLimiterParameter_S_Option, &
                 UseLimiterParameter_F_Option, TimeUnitOption, &
                 FinishTimeOption, CourantFactorOption, nWriteOption )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameter_S_Option, &
      UseLimiterParameter_F_Option
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( I % Type == '' ) &
      I % Type = 'an Integrator_C_1D_MS_C_PS'

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) then
      call Show ( 'Current_BSLL_ASC_CSLD_1D not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_1D_MS_C_PS__Template', 'module', &
                  CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_1D_MS_C_PS', 'subroutine', &
                  CONSOLE % WARNING )
    end if

    call I % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )

    I % UseLimiterParameter_S = .false.
    I % UseLimiterParameter_F = .false.
    if ( present ( UseLimiterParameter_S_Option ) ) &
      I % UseLimiterParameter_S = UseLimiterParameter_S_Option
    if ( present ( UseLimiterParameter_F_Option ) ) &
      I % UseLimiterParameter_F = UseLimiterParameter_F_Option
    call PROGRAM_HEADER % GetParameter &
           ( I % UseLimiterParameter_S, 'UseLimiterParameter_S' )
    call PROGRAM_HEADER % GetParameter &
           ( I % UseLimiterParameter_F, 'UseLimiterParameter_F' )

  end subroutine InitializeTemplate_C_1D_MS_C_PS


  subroutine FinalizeTemplate_C_1D_MS_C_PS ( I )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % Step_MS ) ) &
     deallocate ( I % Step_MS )
   if ( allocated ( I % Step_PS ) ) &
     deallocate ( I % Step_PS )
   if ( allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) &
     deallocate ( I % Current_BSLL_ASC_CSLD_1D )

    call I % FinalizeTemplate_C_PS ( )

  end subroutine FinalizeTemplate_C_1D_MS_C_PS


  subroutine ComputeCycle ( I )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerCycle )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
      call I % ComputeCycle_BSLL_ASC_CSLD_1D ( MS )
    end select

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeCycle


  subroutine InitializeStepTimers ( I, BaseLevel )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( allocated ( I % Step_MS ) ) &
      call I % Step_MS % InitializeTimers ( BaseLevel )
    if ( allocated ( I % Step_PS ) ) &
      call I % Step_PS % InitializeTimers ( BaseLevel )

  end subroutine InitializeStepTimers


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    integer ( KDI ) :: &
      iC  !-- iCurrent
    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    do iC = 1, I % N_CURRENTS_MS
      associate ( CB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
      call CB % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption  = IgnorabilityOption )
      end associate !-- CB
    end do !-- iC

    if ( allocated ( I % Current_ASC ) ) then 
      associate ( CA => I % Current_ASC )
      call CA % ComputeTally &
             ( ComputeChangeOption = ComputeChangeOption, &
               IgnorabilityOption  = IgnorabilityOption )
      end associate !-- CA
    end if

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine ComputeCycle_BSLL_ASC_CSLD_1D ( I, MS )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I
    class ( Bundle_SLL_ASC_CSLD_Form ), intent ( inout ) :: &
      MS

    integer ( KDI ) :: &
      iS  !-- iSection
    real ( KDR ) :: &
      TimeNew

    call I % PrepareCycle ( )
    call I % ComputeNewTime ( TimeNew )

    associate ( TimeStep => TimeNew - I % Time )    

    call I % PrepareStep_MS ( )
    call I % Step_MS % Compute ( I % Time, TimeStep )

    if ( allocated ( I % Step_PS ) ) then
      call I % PrepareStep_PS ( )
      call I % Step_PS % Compute ( I % Time, TimeStep )
    end if
      
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

  end subroutine ComputeCycle_BSLL_ASC_CSLD_1D


  subroutine PrepareStep_MS ( I )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep_MS


  subroutine PrepareStep_PS ( I )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ) :: &
      I

  end subroutine PrepareStep_PS


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    call I % ComputeTimeStepLocalTemplate ( TimeStepCandidate )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeTimeStepLocalTemplate ( I, TimeStepCandidate )

    class ( Integrator_C_1D_MS_C_PS_Template ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC, &  !-- iCurrent
      N_PS
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C
    class ( Current_ASC_Template ), pointer :: &
      CA

    if ( .not. allocated ( I % Current_BSLL_ASC_CSLD_1D ) ) &
      return

    N_PS  =  I % N_CURRENTS_MS + 1

    select type ( PS  =>  I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL  =>  PS % Chart )
    class is ( Chart_SL_Template )

    G  =>  CSL % Geometry ( )

    do iC = 1, N_PS

      if ( iC == N_PS ) then
        if ( allocated ( I % Current_ASC ) ) then
          CA  =>  I % Current_ASC
        else
          CA  =>  null ( )
        end if
      else
        associate ( CB  =>  I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        select type ( CAE => CB % Section % Atlas ( 1 ) % Element )
        class is ( Current_ASC_Template )
          CA => CAE
        end select !-- CAE
        end associate !-- CB
      end if

      if ( associated ( CA ) ) then

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
      else
        TimeStepCandidate ( iC )  =  huge ( 1.0_KDR )
      end if

    end do !-- iC

    end select !-- CSL
    end select !-- PS

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( G, C, CA )

  end subroutine ComputeTimeStepLocalTemplate


end module Integrator_C_1D_MS_C_PS__Template
