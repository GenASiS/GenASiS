module Relaxation_RM_BSLL_ASC_CSLD__Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use RadiationMoments_BSLL_ASC_CSLD__Form
  use Relaxation_RM__Template

  implicit none
  private
  
  type, public, extends ( Relaxation_RM_Template ) :: &
    Relaxation_RM_BSLL_ASC_CSLD_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, nopass :: &
      ApplySubroutine
    final :: &
      Finalize
  end type Relaxation_RM_BSLL_ASC_CSLD_Form

    type ( Relaxation_RM_BSLL_ASC_CSLD_Form ), private, pointer :: &
      Relaxation => null ( )

contains


  subroutine Initialize ( R, RMB, Name )

    class ( Relaxation_RM_BSLL_ASC_CSLD_Form ), intent ( inout ), target :: &
      R
    class ( RadiationMoments_BSLL_ASC_CSLD_Form ), intent ( in ) :: &
      RMB
    character ( * ), intent ( in ) :: &
      Name

    if ( R % Type == '' ) &
      R % Type  =  'a Relaxation_RM_BSLL_ASC_CSLD'

    R % Apply  =>  ApplySubroutine
    call R % InitializeTemplate ( Name )

    select type ( RMA  =>  RMB % Fiber % Atlas ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call R % SetIndices ( RMA )
    end select !-- RMA

    allocate ( R % LinearEquations ( 1 ) )
    associate ( LE => R % LinearEquations ( 1 ) )
    call LE % Initialize ( nEquations = 4, nSolutions = 1 )
    end associate !-- LE

    Relaxation => R

  end subroutine Initialize


  subroutine ApplySubroutine &
               ( S, RadiationMoments, Sources_RM, Increment, Chart, &
                 TimeStep, iStage, GeometryOption, iStrgeometryValueOption )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( CurrentTemplate ), intent ( inout ) :: &
      RadiationMoments
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( StorageForm ), intent ( inout ) :: &
      Increment
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage
    class ( GeometryFlatForm ), intent ( in ), optional :: &
      GeometryOption
    integer ( KDI ), intent ( in ), optional :: &
      iStrgeometryValueOption

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues

    call Show ( 'Relaxation_RM_BSLL_ASC_CSLD % Apply', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    if ( .not. present ( GeometryOption ) &
         .or. .not. present ( iStrgeometryValueOption ) ) &
    then
      call Show ( 'GeometryOption must be present', CONSOLE % ERROR )
      call Show ( 'iStrgeometryValueOption must be present', CONSOLE % ERROR )
      call Show ( 'Apply', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Relaxation_RM_BSLL_ASC_CSLD__Form', 'module', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iStage == 1 ) &
      call Clear ( Sources_RM % Value )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    associate ( I => RM % Interactions )
    call I % Compute ( RM )

    select type ( Chart )
    class is ( Chart_SL_Template )

    nV = size ( Increment % Value, dim = 1 )

    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

      call Relaxation % ComputeIncrements &
             ( S, Sources_RM, Relaxation % LinearEquations ( 1 ), &
               Increment, RadiationMoments, TimeStep, iStage, iV, &
               GeometryOption, iStrgeometryValueOption )

    end do !-- iV

    end select !-- Chart
    end associate !-- I
    end select !-- RM

  end subroutine ApplySubroutine


  subroutine Finalize ( R )

    type ( Relaxation_RM_BSLL_ASC_CSLD_Form ), intent ( inout ) :: &
      R

    call R % FinalizeTemplate ( )

  end subroutine Finalize


end module Relaxation_RM_BSLL_ASC_CSLD__Form
