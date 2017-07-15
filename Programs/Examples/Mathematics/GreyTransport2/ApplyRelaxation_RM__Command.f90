module ApplyRelaxation_RM__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form
  use PrepareRelaxation_RM__Command

  implicit none
  private

  public :: &
    ApplyRelaxation_RM

contains

  
  subroutine ApplyRelaxation_RM &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 RadiationMoments, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      RadiationMoments
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    associate ( ID => S % IncrementDamping )
    call PrepareRelaxation_RM &
           ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
             RadiationMoments, Chart, TimeStep, iStage )
    call ID % Compute &
           ( IncrementExplicit, RadiationMoments, IncrementExplicit, &
             DampingCoefficient, TimeStep )
    end associate !-- iD

  end subroutine ApplyRelaxation_RM


end module ApplyRelaxation_RM__Command
