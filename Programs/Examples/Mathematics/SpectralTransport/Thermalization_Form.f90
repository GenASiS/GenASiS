module Thermalization_Form

  use Basics
  use Mathematics

  type, public, extends ( IntegratorTemplate ) :: ThermalizationForm
  contains
    procedure, public, pass :: &
      Initialize
!    final :: &
!      Finalize  
  end type ThermalizationForm

    real ( KDR ), private, parameter :: &
      T_Scale_MeV = 3.0_KDR!, &
!      T_Min_MeV = 0.1_KDR, &
!      T_Max_MeV = 10.0_KDR

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

    !-- Cleanup

    end select !-- MS
    end select !-- PS

  end subroutine Initialize


end module Thermalization_Form
