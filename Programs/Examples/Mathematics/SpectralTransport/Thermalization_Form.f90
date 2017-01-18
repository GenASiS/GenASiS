module Thermalization_Form

  use Basics
  use Mathematics
  use Matter_Form
  use Matter_ASC__Form

  implicit none
  private

  type, public, extends ( IntegratorTemplate ) :: ThermalizationForm
    type ( Matter_ASC_Form ), allocatable :: &
      Matter_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  end type ThermalizationForm

    private :: &
      SetMatter

    real ( KDR ), private, parameter :: &
      T_Scale_MeV = 3.0_KDR, &
      T_Min_MeV = 0.1_KDR, &
      T_Max_MeV = 10.0_KDR

contains


  subroutine Initialize ( T, Name )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      nEnergyCells
    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    if ( T % Type == '' ) &
      T % Type = 'a Thermalization'

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: T % PositionSpace )
    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( T % Geometry_ASC )
    associate ( GA => T % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: T % MomentumSpace )
    select type ( MS => T % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = T_Scale_MeV * UNIT % MEV

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEV

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    nEnergyCells = 8
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Matter

    allocate ( T % Matter_ASC )
    associate ( MA => T % Matter_ASC )
    call MA % Initialize &
           ( PS, TemperatureUnitOption = UNIT % MEV, &
             ChemicalPotentialUnitOption = UNIT % MEV )
    call SetMatter ( T )
    end associate !-- MA

    !-- Cleanup

    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( T )

    type ( ThermalizationForm ), intent ( inout ) :: &
      T

    if ( allocated ( T % Matter_ASC ) ) &
      deallocate ( T % Matter_ASC )

    call T % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetMatter ( T )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T

    real ( KDR ) :: &
      R_Min, R_Max, &
      T_Min, T_Max
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( MatterForm ), pointer :: &
      M

    M => T % Matter_ASC % Matter ( )

    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    associate &
      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
        X  => G % Value ( :, G % CENTER ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER ( 3 ) ), &
        Temperature => M % Value ( :, M % TEMPERATURE ), &
        ChemicalPotential => M % Value ( :, M % CHEMICAL_POTENTIAL ) )
    associate &
      ( R => sqrt ( ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
                    +  ( Z - R0 ( 3 ) ) ** 2 ) )

    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = minval ( R )
    CO % Outgoing % Value ( 2 ) = 1.0_KDR / maxval ( R )
    call CO % Reduce ( REDUCTION % MIN )

    R_Min = CO % Incoming % Value ( 1 )
    R_Max = 1.0_KDR / CO % Incoming % Value ( 2 )

    T_Min = T_Min_MeV * UNIT % MEV
    T_Max = T_Max_MeV * UNIT % MEV

!    Temperature &
!      = T_Max  -  ( T_Max - T_Min ) / ( R_Max - R_Min ) * ( R - R_Min )
    Temperature = T_Max * ( ( R_Min / R ) ** ( log10 ( T_Max / T_Min ) &
                                               / log10 ( R_Max / R_Min ) ) )
    ChemicalPotential = 0.0_KDR

    end associate !-- R
    end associate !-- R0, etc.
    end select !-- C
    end select !-- PS

    nullify ( G, M )

  end subroutine SetMatter


end module Thermalization_Form
