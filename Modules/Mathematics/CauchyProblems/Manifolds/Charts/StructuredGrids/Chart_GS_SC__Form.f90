module Chart_GS_SC__Form
  
  !-- Chart_GridStructured_SymmetricCurvilinear_Form

  use Basics
  use Chart_GS__Form

  implicit none
  private

  type, public, extends ( Chart_GS_Form ) :: Chart_GS_SC_Form
    integer ( KDI ) :: &
      nCellsRadius = 0
    real ( KDR ) :: &
      RadiusMax
  contains
    procedure, private, pass :: &
      Initialize_GS_SC
    generic, public :: &
      Initialize => Initialize_GS_SC
!    procedure, private, pass :: &
!      Show_C
    final :: &
      Finalize
  end type Chart_GS_SC_Form


contains


  subroutine Initialize_GS_SC &
               ( C, RadiusMax, nDimensions, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, nGhostLayersOption, nCellsRadiusOption )

    class ( Chart_GS_SC_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMax
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS_SC'

    C % RadiusMax  =  RadiusMax

    C % nCellsRadius  =  128
    if ( present ( nCellsRadiusOption ) ) &
      C % nCellsRadius = nCellsRadiusOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsRadius, 'nCellsRadius' )

    ! select case ( Atlas % nDimensions )
    ! case ( 1 )  !-- Spherical coordinates

    !   CoordinateSystem = 'SPHERICAL' 

    !   MinCoordinate = [ 0.0_KDR,   0.0_KDR, 0.0_KDR ]
    !   MaxCoordinate = [ C % RadiusMax, 0.0_KDR, 0.0_KDR ]

    !   nCells = [ C % nCellsRadius, 1, 1 ]

    ! case ( 2 )  !-- Cylindrical coordinates

    !   CoordinateSystem = 'CYLINDRICAL' 

    !   MinCoordinate = [ 0.0_KDR,       - C % RadiusMax, 0.0_KDR ]
    !   MaxCoordinate = [ C % RadiusMax, + C % RadiusMax, 0.0_KDR ]

    !   nCells = [ C % nCellsRadius, 2 * C % nCellsRadius, 1 ]

    ! case ( 3 )  !-- Rectangular coordinates
 
    !   CoordinateSystem = 'RECTANGULAR' 

    !   MinCoordinate  =  - C % RadiusMax
    !   MaxCoordinate  =  + C % RadiusMax

    !   nCells  =  2 * C % nCellsRadius
      
    ! end select !-- nDimensions

  end subroutine Initialize_GS_SC


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_SC_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize


end module Chart_GS_SC__Form
