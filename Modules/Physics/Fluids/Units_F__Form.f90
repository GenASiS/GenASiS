module Units_F__Form
  
  !-- Units_Fluid_Form

  use Basics

  implicit none
  private

  type, public :: Units_F_Form
    !-- Phase space 
    type ( MeasuredValueForm ) :: &
      Time, &
      Length, &
      SqrtDet_M  !-- SquareRoot_Determinant_Metric
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Coordinate_PS, &  !-- Coordinate_PositionSpace
      Coordinate_MS     !-- Coordinate_MomentumSpace
    !-- Local
    type ( MeasuredValueForm ) :: &
      BaryonMass, &
      NumberDensity, &
      MassDensity, &
      EnergyDensity, &
      Temperature
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      Velocity_U, &
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
  end type Units_F_Form

contains


  subroutine Show_U ( U, IgnorabilityOption )

    class ( Units_F_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call Show ( 'Units -- Phase space', IgnorabilityOption )
    call Show ( U % Time,          'Time', IgnorabilityOption )
    call Show ( U % Length,        'Length', IgnorabilityOption )
    call Show ( U % SqrtDet_M,     'SqrtDet_M', IgnorabilityOption )
    call Show ( U % Coordinate_PS, 'Coordinate_PS', IgnorabilityOption )
    call Show ( U % Coordinate_MS, 'Coordinate_MS', IgnorabilityOption )

    call Show ( 'Units -- Local' )
    call Show ( U % BaryonMass,        'BaryonMass', IgnorabilityOption )
    call Show ( U % NumberDensity,     'NumberDensity', IgnorabilityOption )
    call Show ( U % EnergyDensity,     'EnergyDensity', IgnorabilityOption )
    call Show ( U % Temperature,       'Temperature', IgnorabilityOption )
    call Show ( U % Velocity_U,        'Velocity_U', IgnorabilityOption )
    call Show ( U % MomentumDensity_D, 'MomentumDensity_D', &
                IgnorabilityOption )

    call Show ( 'Units -- Global' )
    call Show ( U % Number,          'Number', IgnorabilityOption )
    call Show ( U % Energy,          'Energy', IgnorabilityOption )
    call Show ( U % Momentum,        'Momentum', IgnorabilityOption )
    call Show ( U % AngularMomentum, 'AngularMomentum', IgnorabilityOption )

  end subroutine Show_U


end module Units_F__Form
