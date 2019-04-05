module StressEnergyUnits_Form
  
  use Basics

  implicit none
  private

  type, public :: StressEnergyUnitsForm
    !-- Phase space 
    type ( MeasuredValueForm ) :: &
      Time, &
      Length
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Coordinate_PS, &
      Coordinate_MS
    !-- Local
    type ( MeasuredValueForm ) :: &
      BaryonMass, &
      NumberDensity, &
      EnergyDensity, &
      Temperature
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U, &
      MomentumDensity_U, &
      MomentumDensity_D
    !-- Global
    type ( MeasuredValueForm ) :: &
      Number, &
      Energy, &
      Momentum, &
      AngularMomentum
  end type StressEnergyUnitsForm

end module StressEnergyUnits_Form
