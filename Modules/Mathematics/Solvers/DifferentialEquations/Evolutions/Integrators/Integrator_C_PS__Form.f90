!-- Integrator_C_PS is a class for time evolution of a conserved current
!   on position space.

module Integrator_C_PS__Form

  !-- Integrator_Current_PositionSpace_Template

  use Basics
  use Manifolds
  use Fields
  use EvolutionBasics
  use Steps
  use Integrator_Template
  use TimeSeries_C__Form

  implicit none
  private

  type, public, extends ( IntegratorTemplate ) :: Integrator_C_PS_Form
    real ( KDR ) :: &
      CourantFactor
    class ( Current_ASC_Template ), allocatable :: &
      Current_ASC
    class ( Step_RK_C_ASC_Template ), allocatable :: &
      Step
    class ( TimeSeries_C_Form ), allocatable :: &
      TimeSeries
  contains
    procedure, private, pass :: &  !-- 1
      Initialize_I
    generic, public :: &
      Initialize => Initialize_I
    final :: &  !-- 1
      Finalize
    procedure, private, pass :: &  !-- 2
      ComputeConstraints
    procedure, private, pass :: &  !-- 2
      ComputeCycle
    procedure, private, pass :: &  !-- 3
      InitializeStepTimers
    procedure, private, pass :: &  !-- 3
      ComputeTally
    procedure, private, pass :: &  !-- 3
      RecordTimeSeries
    procedure, private, pass :: &  !-- 3
      WriteTimeSeries
    procedure, private, pass :: &
      ComputeCycle_ASC
    procedure, public, pass :: &
      ComputeTimeStep_C_ASC
    procedure, public, nopass :: &
      ComputeTimeStepKernel_CSL
  end type Integrator_C_PS_Form

    private :: &
      ComputeTimeStepLocal

contains


  subroutine Initialize_I &
               ( I, U, Name, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nWriteOption )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I
    class ( UniverseHeaderForm ), intent ( in ) :: &
      U
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
      I % Type = 'an Integrator_C_PS'

    if ( .not. allocated ( I % PositionSpace ) ) then
      call Show ( 'PositionSpace must be allocated by an extension', &
                  CONSOLE % ERROR )
      call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C_PS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( .not. allocated ( I % Current_ASC ) ) then
      call Show ( 'Current_ASC not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS', 'subroutine', CONSOLE % WARNING )
    end if

    if ( .not. allocated ( I % Step ) ) then
      call Show ( 'Step not allocated by an extension', &
                  CONSOLE % WARNING )
      call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % WARNING )
      call Show ( 'InitializeTemplate_C_PS', 'subroutine', CONSOLE % WARNING )
    end if

    I % CourantFactor = 0.7_KDR
    if ( present ( CourantFactorOption ) ) &
      I % CourantFactor = CourantFactorOption
    call PROGRAM_HEADER % GetParameter &
           ( I % CourantFactor, 'CourantFactor' )

    call I % InitializeTemplate &
           ( U, Name, TimeUnitOption, FinishTimeOption, nWriteOption )
    call Show ( I % CourantFactor, 'CourantFactor', I % IGNORABILITY )

    I % ComputeTimeStepLocal  =>  ComputeTimeStepLocal

    if ( allocated ( I % Current_ASC ) ) then
      allocate ( I % TimeSeries )
      associate &
        ( TS => I % TimeSeries, &
          CA => I % Current_ASC )
      call TS % Initialize &
             ( I, CA % TallyInterior, &
               CA % TallyBoundaryGlobal ( 1 ) % Element, &
               CA % TallyTotal, CA % TallyChange )
      end associate !-- TS, etc.
    end if

  end subroutine Initialize_I


  impure elemental subroutine Finalize ( I )

    type ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I

   if ( allocated ( I % TimeSeries ) ) &
     deallocate ( I % TimeSeries )
   if ( allocated ( I % Step ) ) &
     deallocate ( I % Step )
   if ( allocated ( I % Current_ASC ) ) &
     deallocate ( I % Current_ASC )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeConstraints ( I )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I

    if ( .not. allocated ( I % Step ) ) &
      return

    select type ( S => I % Step )
    class is ( Step_RK_C_ASC_Template )

    if ( associated ( S % ComputeConstraints % Pointer ) ) &
      call S % ComputeConstraints % Pointer ( S )

    end select !-- S

  end subroutine ComputeConstraints


  subroutine ComputeCycle ( I )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I

    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerCycle )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
      call I % ComputeCycle_ASC ( PS )
    end select

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeCycle


  subroutine InitializeStepTimers ( I, BaseLevel )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      BaseLevel

    if ( allocated ( I % Step ) ) &
      call I % Step % InitializeTimers ( BaseLevel )

  end subroutine InitializeStepTimers


  subroutine ComputeTally ( I, ComputeChangeOption, IgnorabilityOption )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I
    logical ( KDL ), intent ( in ), optional :: &
      ComputeChangeOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    type ( TimerForm ), pointer :: &
      Timer

    if ( .not. allocated ( I % Current_ASC ) ) &
      return

    Timer => PROGRAM_HEADER % TimerPointer ( I % iTimerTally )
    if ( associated ( Timer ) ) call Timer % Start ( )

    associate ( CA => I % Current_ASC )
    call CA % ComputeTally &
           ( ComputeChangeOption = ComputeChangeOption, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CA

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ComputeTally


  subroutine RecordTimeSeries ( I, MaxTime, MinTime, MeanTime )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iT  !-- iTimer
    real ( KDR ) :: &
      ReconstructionImbalance

    if ( .not. allocated ( I % TimeSeries ) ) &
      return

    call I % TimeSeries % Record ( MaxTime, MinTime, MeanTime )

    if ( I % TimeSeries % iTime <= 1 ) &
      return  !-- no Reconstruction yet
    if ( PROGRAM_HEADER % Communicator % Rank /= CONSOLE % DisplayRank ) &
      return  !-- only DisplayRank has MinTime, MaxTime

    ! do iT = 1, PROGRAM_HEADER % nTimers
    !   if ( PROGRAM_HEADER % Timer ( iT ) % Name == 'ComputeReconstruction' ) &
    !     exit
    ! end do !-- iT

    ! ReconstructionImbalance &
    !   = ( MaxTime ( iT ) - MinTime ( iT ) ) / MinTime ( iT )
    ! call Show ( ReconstructionImbalance, 'ReconstructionImbalance', &
    !             I % IGNORABILITY + 2 )
    ! if ( ReconstructionImbalance > 0.5_KDR .and. MaxTime ( iT ) > 10.0 ) then
    !   call Show ( 'ReconstructionBalance > 0.5', CONSOLE % ERROR )
    !   call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % ERROR )
    !   call Show ( 'RecordTimeSeries', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end if

  end subroutine RecordTimeSeries


  subroutine WriteTimeSeries ( I )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I

    if ( .not. allocated ( I % TimeSeries ) ) &
      return

    call I % TimeSeries % Write ( )

  end subroutine WriteTimeSeries


  subroutine ComputeCycle_ASC ( I, PS )

    class ( Integrator_C_PS_Form ), intent ( inout ) :: &
      I
    class ( Atlas_SC_Form ), intent ( inout ) :: &
      PS

    real ( KDR ) :: &
      TimeNew

    call I % PrepareCycle ( )
    call I % ComputeNewTime ( TimeNew )

    associate ( S => I % Step )
    associate ( TimeStep => TimeNew - I % Time )    

    select type ( Chart => PS % Chart )
    class is ( Chart_SLD_Form )

      call S % Compute ( I % Time, TimeStep )

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'Integrator_C_PS__Template', 'module', CONSOLE % ERROR )
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

    end associate !-- TimeStep
    end associate !-- S

  end subroutine ComputeCycle_ASC


  subroutine ComputeTimeStep_C_ASC ( I, TimeStepCandidate, CA )

    class ( Integrator_C_PS_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), intent ( inout ) :: &
      TimeStepCandidate
    class ( Current_ASC_Template ), intent ( in ) :: &
      CA

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( CurrentTemplate ), pointer :: &
      C

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    G => CSL % Geometry ( )
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
             CSL % nDimensions, TimeStepCandidate )

    end select !-- CSL
    end select !-- PS

    TimeStepCandidate = I % CourantFactor * TimeStepCandidate

    nullify ( C, G )

  end subroutine ComputeTimeStep_C_ASC


  subroutine ComputeTimeStepKernel_CSL &
               ( IsProperCell, FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 dXL_1, dXL_2, dXL_3, dXR_1, dXR_2, dXR_3, Crsn_2, Crsn_3, &
                 nDimensions, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      dXL_1, dXL_2, dXL_3, &
      dXR_1, dXR_2, dXR_3, &
      Crsn_2, Crsn_3
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
        = maxval ( max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ), &
                   mask = IsProperCell )
    case ( 2 )
      TimeStepInverse &
        = maxval (   max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ) &
                   + max ( FEP_2, -FEM_2 ) / ( Crsn_2 * ( dXL_2 + dXR_2 ) ), &
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
                      / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                    + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                      / ( Crsn_2 ( iV ) * ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) &
                    + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                      / ( Crsn_3 ( iV ) * ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
      end do
      !$OMP end parallel do
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStepKernel_CSL


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( IntegratorTemplate ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    select type ( I )
    class is ( Integrator_C_PS_Form )

    call I % ComputeTimeStep_C_ASC &
           ( TimeStepCandidate ( 1 ), I % Current_ASC )

    end select !-- I

  end subroutine ComputeTimeStepLocal


end module Integrator_C_PS__Form
