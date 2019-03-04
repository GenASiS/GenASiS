module StressEnergyUnits_Form
  
  use Basics

  implicit none
  private

  type, public :: StressEnergyUnitsForm
    !-- Spacetime
    type ( MeasuredValueForm ) :: &
      Length, &
      Time
    !-- Local
    type ( MeasuredValueForm ) :: &
      BaryonMassUnit, &
      NumberDensityUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U_Unit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit
    !-- Global
    type ( MeasuredValueForm ) :: &
      NumberUnit, &
      EnergyUnit, &
      MomentumUnit, &
      AngularMomentumUnit
  end type StressEnergyUnitsForm

end module StressEnergyUnits_Form
