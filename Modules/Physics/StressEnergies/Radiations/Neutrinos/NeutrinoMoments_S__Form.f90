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
    procedure, public, pass :: &
      InitializeAllocate_NM
    final :: &
      Finalize
  end type NeutrinoMoments_S_Form

contains


  subroutine InitializeAllocate_NM &
               ( RM, NeutrinoType, RiemannSolverType, ReconstructedType, &
                 UseLimiter, Units, LimiterParameter, nValues, VariableOption, &
                 VectorOption, NameOption, ClearOption, UnitOption, &
                 VectorIndicesOption )

    class ( NeutrinoMoments_S_Form ), intent ( inout ) :: &
      RM
    character ( * ), intent ( in ) :: &
      NeutrinoType, &
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
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ), optional ::&
      VectorIndicesOption

    RM % Type = NeutrinoType

    call RM % RadiationMomentsForm % InitializeAllocate_RM &
           ( RiemannSolverType, ReconstructedType, UseLimiter, Units, &
             LimiterParameter, nValues, VariableOption = VariableOption, &
             VectorOption = VectorOption, NameOption = NameOption, &
             ClearOption = ClearOption, UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_NM


  impure elemental subroutine Finalize ( PM )

    type ( NeutrinoMoments_S_Form ), intent ( inout ) :: &
      PM

    !-- Trigger finalization in parent

  end subroutine Finalize


end module NeutrinoMoments_S__Form
