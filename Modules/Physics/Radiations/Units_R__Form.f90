module Units_R__Form
  
  !-- Units_Radiation__Form

  use Basics
  use Fluids

  implicit none
  private

  type, public, extends ( Units_F_Form ) :: Units_R_Form
!-- FIXME: GCC 11.3 doesn't like hardwired dimensionality
!    type ( MeasuredValueForm ), dimension ( 3 ) :: &
    !-- Local
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      Coordinate_MS, &     !-- Coordinate_MomentumSpace
      MomentumDensity_U
    !-- Global
    type ( QuantityForm ) :: &
      Luminosity, &
      EnergyAverage
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_U
    final :: &
      Finalize
  end type Units_R_Form

contains


  subroutine Initialize ( U, CoordinateUnit, TypeOption )

    class ( Units_R_Form ), intent ( inout ) :: &
      U
    type ( QuantityForm ), dimension ( : ), intent ( in ) :: &
      CoordinateUnit
    character ( * ), intent ( in ), optional :: &
      TypeOption

    allocate ( U % Coordinate_MS ( 3 ) )
    allocate ( U % MomentumDensity_U ( 3 ) )

    call U % Units_F_Form % Initialize ( CoordinateUnit, TypeOption )

    if ( present ( TypeOption ) ) then
      select case ( trim ( TypeOption ) )
      case ( '' )
      case ( 'MKS' )
        !-- Local
        U % MomentumDensity_U ( 1 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_U ( 2 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_U ( 3 )  =  UNIT % MASS_DENSITY_MKS  &
                                        *  U % Velocity_U ( 3 )
      case ( 'CGS' )
        !-- Local
        U % MomentumDensity_U ( 1 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_U ( 2 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_U ( 3 )  =  UNIT % MASS_DENSITY_CGS  &
                                        *  U % Velocity_U ( 3 )
      case ( 'ASTROPHYSICS' )
        !-- Local
        U % MomentumDensity_U ( 1 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Velocity_U ( 1 )
        U % MomentumDensity_U ( 2 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Velocity_U ( 2 )
        U % MomentumDensity_U ( 3 )  =  UNIT % MASS_DENSITY_NUCLEAR  &
                                        *  U % Velocity_U ( 3 )
        !-- Global
        U % Luminosity     =  UNIT % BETHE  /  UNIT % SECOND
        U % EnergyAverage  =  UNIT % MEGA_ELECTRON_VOLT
      case default
        call Show ( 'Type not recognized', CONSOLE % ERROR )
        call Show ( 'Units_R__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- TypeOption
    end if !-- TypeOption

  end subroutine Initialize


  subroutine Show_U ( U, IgnorabilityOption )

    class ( Units_R_Form ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call Show ( 'Units -- Phase space', IgnorabilityOption )
    call Show ( U % Time,          'Time',          IgnorabilityOption )
    call Show ( U % Length,        'Length',        IgnorabilityOption )
    call Show ( U % SqrtDet_M,     'SqrtDet_M',     IgnorabilityOption )
    call Show ( U % Coordinate_PS, 'Coordinate_PS', IgnorabilityOption )
    call Show ( U % Coordinate_MS, 'Coordinate_MS', IgnorabilityOption )

    call Show ( 'Units -- Local' )
    call Show ( U % NumberDensity,     'NumberDensity', IgnorabilityOption )
    call Show ( U % EnergyDensity,     'EnergyDensity', IgnorabilityOption )
    call Show ( U % Temperature,       'Temperature',   IgnorabilityOption )
    call Show ( U % MomentumDensity_U, 'MomentumDensity_U', &
                                                        IgnorabilityOption )
    call Show ( U % MomentumDensity_D, 'MomentumDensity_D', &
                                                        IgnorabilityOption )

    call Show ( 'Units -- Global' )
    call Show ( U % Number,          'Number', IgnorabilityOption )
    call Show ( U % Energy,          'Energy', IgnorabilityOption )
    call Show ( U % Momentum,        'Momentum', IgnorabilityOption )
    call Show ( U % AngularMomentum, 'AngularMomentum', IgnorabilityOption )

  end subroutine Show_U


  impure elemental subroutine Finalize ( U )

    type ( Units_R_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % MomentumDensity_U ) ) &
      deallocate ( U % MomentumDensity_U )
    if ( allocated ( U % Coordinate_MS ) ) &
      deallocate ( U % Coordinate_MS )

  end subroutine Finalize


end module Units_R__Form
