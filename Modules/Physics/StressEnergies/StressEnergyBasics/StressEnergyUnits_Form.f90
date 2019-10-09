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
  contains
    procedure, public, pass :: &
      Show => Show_U
  end type StressEnergyUnitsForm

contains


  subroutine Show_U ( SEU, IgnorabilityOption )

    class ( StressEnergyUnitsForm ), intent ( in ) :: &
      SEU
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call Show ( 'Units -- Phase space', IgnorabilityOption )
    call Show ( SEU % Time,          'Time', IgnorabilityOption )
    call Show ( SEU % Length,        'Length', IgnorabilityOption )
    call Show ( SEU % Coordinate_PS, 'Coordinate_PS', IgnorabilityOption )
    call Show ( SEU % Coordinate_MS, 'Coordinate_MS', IgnorabilityOption )

    call Show ( 'Units -- Local' )
    call Show ( SEU % BaryonMass,        'BaryonMass', IgnorabilityOption )
    call Show ( SEU % NumberDensity,     'NumberDensity', IgnorabilityOption )
    call Show ( SEU % EnergyDensity,     'EnergyDensity', IgnorabilityOption )
    call Show ( SEU % Temperature,       'Temperature', IgnorabilityOption )
    call Show ( SEU % Velocity_U,        'Velocity_U', IgnorabilityOption )
    call Show ( SEU % MomentumDensity_U, 'MomentumDensity_U', &
                IgnorabilityOption )
    call Show ( SEU % MomentumDensity_D, 'MomentumDensity_D', &
                IgnorabilityOption )

    call Show ( 'Units -- Global' )
    call Show ( SEU % Number,          'Number', IgnorabilityOption )
    call Show ( SEU % Energy,          'Energy', IgnorabilityOption )
    call Show ( SEU % Momentum,        'Momentum', IgnorabilityOption )
    call Show ( SEU % AngularMomentum, 'AngularMomentum', IgnorabilityOption )

  end subroutine Show_U


end module StressEnergyUnits_Form
