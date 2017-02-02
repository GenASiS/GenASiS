!-- Integrator_C is a template for time evolution of a conserved current.

module Integrator_C__Template

  !-- Integrator_Current_Template

  use Basics
  use Manifolds
  use Fields
  use Steps
  use Integrator_Template
  use TimeSeries_C__Form

  implicit none
  private

  type, public, extends ( IntegratorTemplate ), abstract :: &
    Integrator_C_Template
      real ( KDR ) :: &
        CourantFactor
      logical ( KDL ) :: &
        UseLimiterParameter
      class ( Current_ASC_Template ), allocatable :: &
        Current_ASC
      class ( Step_RK_Template ), allocatable :: &
        Step
      class ( TimeSeries_C_Form ), allocatable :: &
        TimeSeries
  contains
    procedure, public, pass :: &  !-- 1
      InitializeTemplate_C
    procedure, public, pass :: &  !-- 1
      FinalizeTemplate_C
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &  !-- 3
      RecordTimeSeries
    procedure, private, pass :: &  !-- 3
      WriteTimeSeries
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type Integrator_C_Template

    private :: &
      ComputeCycle_ASC

      private :: &
        ComputeCycle_ASC_CSL

      private :: &
        ComputeTimeStepKernel_CSL

contains


  subroutine InitializeTemplate_C &
               ( I, Name, UseLimiterParameterOption, TimeUnitOption, &
                 FinishTimeOption, nWriteOption )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I
    character ( * ), intent ( in )  :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      UseLimiterParameterOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    if ( .not. allocated ( I % PositionSpace ) ) then
      call Show ( 'PositionSpace must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( .not. allocated ( I % Current_ASC ) ) then
      call Show ( 'Current must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( .not. allocated ( I % Step ) ) then
      call Show ( 'Step must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    I % CourantFactor = 0.7
    call PROGRAM_HEADER % GetParameter &
           ( I % CourantFactor, 'CourantFactor' )

    I % UseLimiterParameter = .true.
    if ( present ( UseLimiterParameterOption ) ) &
      I % UseLimiterParameter = UseLimiterParameterOption
    call PROGRAM_HEADER % GetParameter &
           ( I % UseLimiterParameter, 'UseLimiterParameter' )

    call I % InitializeTemplate &
           ( Name, TimeUnitOption, FinishTimeOption, nWriteOption )

    allocate ( I % TimeSeries )
    associate &
      ( TS => I % TimeSeries, &
        CA => I % Current_ASC )
    call TS % Initialize &
           ( I, CA % TallyInterior, &
             CA % TallyBoundaryGlobal ( 1 ) % Element, &
             CA % TallyTotal, CA % TallyChange )
    end associate !-- TS, etc.

  end subroutine InitializeTemplate_C


  subroutine FinalizeTemplate_C ( I )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I

   if ( allocated ( I % TimeSeries ) ) &
     deallocate ( I % TimeSeries )
   if ( allocated ( I % Step ) ) &
     deallocate ( I % Step )
   if ( allocated ( I % Current_ASC ) ) &
     deallocate ( I % Current_ASC )

    call I % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_C


  subroutine ComputeCycle ( I )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I

    associate ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeCycle ) )
    call Timer % Start ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call ComputeCycle_ASC ( I, PS )
    end select

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeCycle


  subroutine ComputeTally ( I, ComputeChangeOption )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption

    associate ( Timer => PROGRAM_HEADER % Timer ( I % iTimerComputeTally ) )
    call Timer % Start ( )

    associate ( CA => I % Current_ASC )
    call CA % ComputeTally ( ComputeChangeOption = ComputeChangeOption )
    end associate !-- CA

    call Timer % Stop ( )
    end associate !-- Timer

  end subroutine ComputeTally


  subroutine RecordTimeSeries ( I, MaxTime, MinTime, MeanTime )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iT  !-- iTimer
    real ( KDR ) :: &
      ReconstructionImbalance

    call I % TimeSeries % Record ( MaxTime, MinTime, MeanTime )

    if ( I % TimeSeries % iTime <= 1 ) &
      return  !-- no Reconstruction yet
    if ( PROGRAM_HEADER % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return  !-- only DisplayRank has MinTime, MaxTime

    do iT = 1, PROGRAM_HEADER % nTimers
      if ( PROGRAM_HEADER % Timer ( iT ) % Name == 'ComputeReconstruction' ) &
        exit
    end do !-- iT

    ReconstructionImbalance &
      = ( MaxTime ( iT ) - MinTime ( iT ) ) / MinTime ( iT )
    call Show ( ReconstructionImbalance, 'ReconstructionImbalance' )
    if ( ReconstructionImbalance > 0.5_KDR .and. MaxTime ( iT ) > 10.0 ) then
      call Show ( 'ReconstructionBalance > 0.5', CONSOLE % ERROR )
      call Show ( 'Integrator_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'RecordTimeSeries', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine RecordTimeSeries


  subroutine WriteTimeSeries ( I )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I

    call I % TimeSeries % Write ( )

  end subroutine WriteTimeSeries


  subroutine ComputeTimeStepLocal ( I, TimeStep )

    class ( Integrator_C_Template ), intent ( in ) :: &
      I
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    associate ( CA => I % Current_ASC )

    G => CSL % Geometry ( )
    C => CA % Current ( )
    
    call ComputeTimeStepKernel_CSL &
           ( CSL % IsProperCell, &
             C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 1 ) ), &
             C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 2 ) ), &
             C % Value ( :, C % FAST_EIGENSPEED_PLUS ( 3 ) ), &
             C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 1 ) ), &
             C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 2 ) ), &
             C % Value ( :, C % FAST_EIGENSPEED_MINUS ( 3 ) ), &
             G % Value ( :, G % WIDTH ( 1 ) ), &
             G % Value ( :, G % WIDTH ( 2 ) ), & 
             G % Value ( :, G % WIDTH ( 3 ) ), &
             CSL % nDimensions, TimeStep )

    end associate !-- CA
    end select !-- CSL
    end select !-- PS

    TimeStep = I % CourantFactor * TimeStep

    nullify ( C, G )

  end subroutine ComputeTimeStepLocal


  subroutine ComputeCycle_ASC ( I, PS )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I
    class ( Atlas_SC_Form ), intent ( inout ) :: &
      PS

    select type ( Chart => PS % Chart )
    class is ( Chart_SLD_Form )
      call ComputeCycle_ASC_CSL ( I, I % Current_ASC, Chart )
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_C__Template', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeCycle_ASC', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine ComputeCycle_ASC


  subroutine ComputeCycle_ASC_CSL ( I, CA, CSL )

    class ( Integrator_C_Template ), intent ( inout ) :: &
      I
    class ( Current_ASC_Template ), intent ( inout ) :: &
      CA
    class ( Chart_SLD_Form ), intent ( inout ) :: &
      CSL

    real ( KDR ) :: &
      TimeNew

    select type ( S => I % Step )
    class is ( Step_RK_C_Template )

    select type ( FC => CA % Chart )
    class is ( Field_CSL_Template )

    select type ( C => FC % Field )
    class is ( CurrentTemplate )

    call I % ComputeNewTime ( TimeNew )
    associate ( TimeStep => TimeNew - I % Time )    

    call S % Compute &
           ( C, CSL, I % Time, TimeStep, &
             UseLimiterParameterOption = I % UseLimiterParameter )

    I % iCycle = I % iCycle + 1
    I % Time = I % Time + TimeStep
    if ( I % Time == I % WriteTime ) &
      I % IsCheckpointTime = .true.

    call CA % AccumulateBoundaryTally ( S % BoundaryFluence_CSL )

    end associate !-- TimeStep
    end select !-- C
    end select !-- FC
    end select !-- S

  end subroutine ComputeCycle_ASC_CSL


  subroutine ComputeTimeStepKernel_CSL &
               ( IsProperCell, FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 dX_1, dX_2, dX_3, nDimensions, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      dX_1, dX_2, dX_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    real ( KDR ), intent ( inout ) :: &
      TimeStep

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      TimeStepInverse

    nV = size ( FEP_1 )

    select case ( nDimensions )
    case ( 1 )
      TimeStepInverse &
        = maxval ( max ( FEP_1, -FEM_1 ) / dX_1, mask = IsProperCell )
    case ( 2 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
                   + max ( FEP_2, -FEM_2 ) / dX_2, &
                   mask = IsProperCell )
    case ( 3 )
      ! TimeStepInverse &
      !   = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
      !              + max ( FEP_2, -FEM_2 ) / dX_2 &
      !              + max ( FEP_3, -FEM_3 ) / dX_3, &
      !              mask = IsProperCell )
      TimeStepInverse = - huge ( 0.0_KDR )
      !$OMP parallel do private ( iV ) &
      !$OMP reduction ( max : TimeStepInverse )
      do iV = 1, nV
        if ( IsProperCell ( iV ) ) &
          TimeStepInverse &
            = max ( TimeStepInverse, &
                      max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                      / dX_1 ( iV ) &
                    + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                      / dX_2 ( iV ) &
                    + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                      / dX_3 ( iV ) )
      end do
      !$OMP end parallel do
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStepKernel_CSL


end module Integrator_C__Template
