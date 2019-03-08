module PhotonMoments_S__Form

  !-- PhotonMoments_Spectral__Form
  !-- Adds nothing but a name to RadiationMoments_Form

  use Basics
  use Mathematics
  use RadiationBasics

  implicit none
  private

  type, public, extends ( RadiationMomentsForm ) :: PhotonMoments_S_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM
    final :: &
      Finalize
  end type PhotonMoments_S_Form

contains


  subroutine InitializeAllocate_RM &
               ( RM, RiemannSolverType, UseLimiter, Velocity_U_Unit, &
                 MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
                 EnergyDensityUnit, LimiterParameter, &
                 nValues, VariableOption, VectorOption, NameOption, &
                 ClearOption, UnitOption, VectorIndicesOption )

    class ( PhotonMoments_S_Form ), intent ( inout ) :: &
      RM
    character ( * ), intent ( in ) :: &
      RiemannSolverType
    logical ( KDL ), intent ( in ) :: &
      UseLimiter
    type ( MeasuredValueForm ), dimension ( 3 ), intent ( in ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyDensityUnit
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

    if ( RM % Type == '' ) &
      RM % Type = 'PhotonMoments_S'

    call RM % RadiationMomentsForm % Initialize &
           ( RiemannSolverType, UseLimiter, Velocity_U_Unit, &
             MomentumDensity_U_Unit, MomentumDensity_D_Unit, &
             EnergyDensityUnit, LimiterParameter, nValues, &
             VariableOption = VariableOption, VectorOption = VectorOption, &
             NameOption = NameOption, ClearOption = ClearOption, &
             UnitOption = UnitOption, &
             VectorIndicesOption = VectorIndicesOption )

  end subroutine InitializeAllocate_RM


  impure elemental subroutine Finalize ( PM )

    type ( PhotonMoments_S_Form ), intent ( inout ) :: &
      PM

    !-- Trigger finalization in parent

  end subroutine Finalize


end module PhotonMoments_S__Form
