module NeutrinoMoments_S__Form

  !-- NeutrinoMoments_Spectral__Form
  !-- Adds nothing but a name to RadiationMoments_Form

  use Basics
  use Mathematics
  use StressEnergyBasics
  use RadiationBasics

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: NeutrinoMoments_S_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    final :: &
      Finalize
  end type NeutrinoMoments_S_Form

contains


  subroutine InitializeAllocate_RM &
               ( RM, RadiationMomentsType, RiemannSolverType, &
                 ReconstructedType, UseLimiter, Units, LimiterParameter, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, PinnedOption, UnitOption, VectorIndicesOption )

    class ( NeutrinoMoments_S_Form ), intent ( inout ) :: &
      RM
    character ( * ), intent ( in ) :: &
      RadiationMomentsType, &
      RiemannSolverType, &
      ReconstructedType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      Units
    real ( KDR ), intent ( in ) :: &
      LimiterParameter
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption, &
      VectorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption, &
      PinnedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    if ( RM % Type == '' ) &
      RM % Type = 'NeutrinoMoments_S'

    call RM % RadiationMomentsForm % Initialize &
           ( RadiationMomentsType, RiemannSolverType, ReconstructedType, &
             UseLimiter, Units, LimiterParameter, nValues, &
             VariableOption = VariableOption, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             PinnedOption = PinnedOption, UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_RM


  impure elemental subroutine Finalize ( NM )

    type ( NeutrinoMoments_S_Form ), intent ( inout ) :: &
      NM

    !-- Trigger finalization in parent

  end subroutine Finalize


end module NeutrinoMoments_S__Form
