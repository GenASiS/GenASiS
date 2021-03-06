module DensityWaveStep_Form

  use Basics
  use ProtoFields
  use EvolutionBasics
  use ProtoIncrements
  use Step_RK2_C_ASC__Form

  implicit none
  private

  type, public, extends ( DensityWaveIncrementForm ) :: DensityWaveStepForm
    type ( Step_RK2_C_ASC_Form ), allocatable :: &
      Step
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type DensityWaveStepForm

contains


  subroutine Initialize ( DW, Name )

    class ( DensityWaveStepForm ), intent ( inout ) :: &
      DW
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      CourantFactor, &
      TimeStep
    type ( IntegratorHeaderForm ) :: &
      I

    !-- Initialize parent

    call DW % DensityWaveForm % Initialize ( Name )
    call DW % ProtoCurrent % ComputeTally ( ComputeChangeOption = .false. )
    call DW % Write ( )

    !-- Initialize values needed for Step

    select case ( DW % Atlas % nDimensions )
    case ( 1 )
      CourantFactor = 0.7
    case ( 2 ) 
      CourantFactor = 0.4
    case ( 3 )
      CourantFactor = 0.25
    end select
    call PROGRAM_HEADER % GetParameter &
           ( CourantFactor, 'CourantFactor' )

    call DW % ComputeTimeStep ( TimeStep )
    TimeStep = CourantFactor * TimeStep
    call Show ( TimeStep, 'TimeStep', DW % IGNORABILITY )

    !-- Step

    allocate ( DW % Step )
    associate ( S => DW % Step )

    call S % Initialize ( I, DW % ProtoCurrent, Name )

    call S % Compute ( DW % Time, TimeStep ) 
    call DW % ProtoCurrent % ComputeTally ( )
    call DW % Write ( )

    end associate !-- S

  end subroutine Initialize


  subroutine Finalize ( DW )

    type ( DensityWaveStepForm ), intent ( inout ) :: &
      DW

    if ( allocated ( DW % Step ) ) &
      deallocate ( DW % Step )

  end subroutine Finalize


end module DensityWaveStep_Form
