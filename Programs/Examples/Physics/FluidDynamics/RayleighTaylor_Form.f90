module RayleighTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RayleighTaylorForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RayleighTaylorForm

contains


  subroutine Initialize ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      Acceleration
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace

    if ( RT % Type == '' ) &
      RT % Type = 'a RayleighTaylor'

    call RT % InitializeTemplate ( Name )

    Acceleration = 0.1_KDR
    call PROGRAM_HEADER % GetParameter ( Acceleration, 'Acceleration' )

    !-- Integrator

    allocate ( FluidBoxForm :: RT % Integrator )
    select type ( FB => RT % Integrator )
    type is ( FluidBoxForm )

    associate ( BCF => BoundaryConditionsFace )
    call BCF ( 1 ) % Initialize ( [ 'PERIODIC', 'PERIODIC' ] )
    MinCoordinate ( 1 ) = -0.25_KDR
    MaxCoordinate ( 1 ) = +0.25_KDR
    nCells ( 1 ) = 64
    select case ( trim ( PROGRAM_HEADER % Dimensionality ) )
    case ( '2D' )
      MinCoordinate ( 2 ) = -0.75_KDR
      MaxCoordinate ( 2 ) = +0.75_KDR
      nCells ( 2 ) = 192
      call BCF ( 2 ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )
    case ( '3D' )
      MinCoordinate ( 2 ) = -0.25_KDR
      MaxCoordinate ( 2 ) = +0.25_KDR
      MinCoordinate ( 3 ) = -0.75_KDR
      MaxCoordinate ( 3 ) = +0.75_KDR
      nCells ( 2 ) = 64
      nCells ( 3 ) = 192
      call BCF ( 2 ) % Initialize ( [ 'PERIODIC', 'PERIODIC' ] )
      call BCF ( 3 ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )
    end select

    call FB % Initialize &
           ( Name, FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', &
             BoundaryConditionsFaceOption = BCF, &
             GravitySolverTypeOption = 'UNIFORM', &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             UniformAccelerationOption = Acceleration, &
             nCellsOption = nCells )

    end associate !-- BCF

    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )


    !-- Cleanup

    end select !-- PS
    end select !-- FB

  end subroutine Initialize


  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    call RT % FinalizeTemplate ( )

  end subroutine Finalize


end module RayleighTaylor_Form
