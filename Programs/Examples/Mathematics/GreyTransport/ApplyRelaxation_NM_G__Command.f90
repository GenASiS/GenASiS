module ApplyRelaxation_NM_G__Command

  use Basics
  use Mathematics
  use NeutrinoMoments_G__Form
  use Sources_RM__Form
  use PrepareRelaxation_NM_G__Command
  use ApplyRelaxation_RM__Command

  implicit none
  private

  public :: &
    ApplyRelaxation_NM_G

contains

  
  subroutine ApplyRelaxation_NM_G &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 NeutrinoMoments_G, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      NeutrinoMoments_G
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iNumber

    associate ( ID => S % IncrementDamping )

    !-- Number

    call PrepareRelaxation_NM_G_Number &
           ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
             NeutrinoMoments_G, Chart, TimeStep, iStage, iNumber )
    call ID % Compute &
           ( IncrementExplicit, NeutrinoMoments_G, IncrementExplicit, &
             DampingCoefficient, TimeStep, [ iNumber ] )

    !-- Energy and Momentum

    call ApplyRelaxation_RM &
           ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
             NeutrinoMoments_G, Chart, TimeStep, iStage )

    end associate !-- iD

  end subroutine ApplyRelaxation_NM_G


end module ApplyRelaxation_NM_G__Command
