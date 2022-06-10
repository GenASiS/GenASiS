module Relaxation_RM_ASC__Form

  use OMP_LIB
  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use Relaxation_RM__Template

  implicit none
  private
  
  type, public, extends ( Relaxation_RM_Template ) :: Relaxation_RM_ASC_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, nopass :: &
      ApplySubroutine
    final :: &
      Finalize
  end type Relaxation_RM_ASC_Form

    type ( Relaxation_RM_ASC_Form ), private, pointer :: &
      Relaxation => null ( )

contains


  subroutine Initialize ( R, RMA, Name )

    class ( Relaxation_RM_ASC_Form ), intent ( inout ), target :: &
      R
    class ( RadiationMoments_ASC_Form ), intent ( in ) :: &
      RMA
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iT  !-- iThread

    if ( R % Type == '' ) &
      R % Type  =  'a Relaxation_RM_ASC'

    R % Apply  =>  ApplySubroutine
    call R % InitializeTemplate ( Name )

    call R % SetIndices ( RMA )

    allocate ( R % LinearEquations ( 0 : PROGRAM_HEADER % MaxThreads - 1 ) )
    do iT = 0, size ( R % LinearEquations ) - 1
      associate ( LE => R % LinearEquations ( iT ) )
      call LE % Initialize ( nEquations = 4, nSolutions = 1 )
      end associate !-- LE
    end do !-- iT

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
      nV, &  !-- nValues
      iThread
    class ( GeometryFlatForm ), pointer :: &
      G

    call Show ( 'Relaxation_RM_ASC % Apply', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    if ( iStage == 1 ) &
      call Clear ( Sources_RM % Value )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    associate ( I => RM % Interactions )
    call I % Compute ( RM )

    select type ( Chart )
    class is ( Chart_SL_Template )
    G => Chart % Geometry ( )

    nV = size ( Increment % Value, dim = 1 )

    !$OMP parallel do private ( iV )
    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

      iThread = OMP_GET_THREAD_NUM ( )

      call Relaxation % ComputeIncrements &
             ( S, Sources_RM, Relaxation % LinearEquations ( iThread ), &
               Increment, RadiationMoments, TimeStep, iStage, iV, G )

    end do !-- iV
    !$OMP end parallel do


    end select !-- Chart
    end associate !-- I
    end select !-- RM
    nullify ( G )

  end subroutine ApplySubroutine


  subroutine Finalize ( R )

    type ( Relaxation_RM_ASC_Form ), intent ( inout ) :: &
      R

    call R % FinalizeTemplate ( )

  end subroutine Finalize


end module Relaxation_RM_ASC__Form
