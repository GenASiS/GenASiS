module LaplacianMultipole_Form_Test__Form

  use Basics
  use Manifolds
  use LaplacianMultipole_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'LaplacianMultipole_Form_Test'
    
  type, public :: LaplacianMultipole_Form_Test_Form
    type ( GeometryFlat_ASC_Form ), allocatable :: &
      Geometry
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type LaplacianMultipole_Form_Test_Form

contains


  subroutine Initialize ( LMFT, Name )

    class ( LaplacianMultipole_Form_Test_Form ) :: &
      LMFT
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      RadiusMax
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDL ) :: &
      CoordinateSystem

    !-- Atlas

    allocate ( LMFT % Atlas )
    associate ( A => LMFT % Atlas )
    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )

    CoordinateSystem = 'SPHERICAL'

    nCellsRadial = 64
    call PROGRAM_HEADER % GetParameter ( nCellsRadial, 'nCellsRadial' )

        nCellsPolar = nCellsRadial / 2
    nCellsAzimuthal = nCellsPolar / 2
 
    nCells = [ nCellsRadial, 1, 1 ]
    if ( A % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( A % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    RadiusMax = 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   0.0_KDR, 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    call A % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             nCellsOption = nCells )

    !-- Geometry

    allocate ( LMFT % Geometry )
    associate ( GA => LMFT % Geometry )  !-- GeometryAtlas
    call GA % Initialize ( A )
    call A % SetGeometry ( GA )
    end associate !-- GA

    end associate !-- A

  end subroutine Initialize


  subroutine Finalize ( LMFT )

    type ( LaplacianMultipole_Form_Test_Form ) :: &
      LMFT

    if ( allocated ( LMFT % Geometry ) ) &
      deallocate ( LMFT % Geometry )
    if ( allocated ( LMFT % Atlas ) ) &
      deallocate ( LMFT % Atlas )

  end subroutine Finalize


end module LaplacianMultipole_Form_Test__Form


program LaplacianMultipole_Form_Test

  use Basics
  use LaplacianMultipole_Form_Test__Form

  implicit none

  type ( LaplacianMultipole_Form_Test_Form ), allocatable :: &
    LMFT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( LMFT )
  call LMFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( LMFT )

  deallocate ( PROGRAM_HEADER )

end program LaplacianMultipole_Form_Test

