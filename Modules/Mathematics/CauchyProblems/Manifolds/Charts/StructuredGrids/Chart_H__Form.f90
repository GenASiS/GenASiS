module Chart_H__Form

  !-- Chart_Header__Form

  use Basics
  use Connectivity_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3

  type, public :: Chart_H_Form
    integer ( KDI ) :: &
      IGNORABILITY, &
      nDimensions
    type ( QuantityForm ), dimension ( MAX_DIMENSIONS ) :: &
      CoordinateUnit
    logical ( KDL ), dimension ( MAX_DIMENSIONS ) :: &
      Periodic = .false.
    character ( LDL ) :: &
      Type = '', &
      Name, &
      CoordinateSystem
    character ( LDL ), dimension ( MAX_DIMENSIONS ) :: &
      CoordinateLabel
    type ( ConnectivityForm ), allocatable :: &
      Connectivity
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, private, pass :: &
      Show_C
    generic, public :: &
      Show => Show_C
    final :: &
      Finalize_C
  end type Chart_H_Form

  type, public :: ChartElement
    !-- Chart_Element_Form
    class ( Chart_H_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type ChartElement

    private :: &
      SetDimensionality, &
      SetCoordinateSystem


contains


  subroutine Initialize_H &
               ( C, CoordinateLabelOption, CoordinateSystemOption, NameOption, &
                 CoordinateUnitOption, IgnorabilityOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption, &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nDimensionsOption, &
      iDimensionalityOption

    C % IGNORABILITY  =  CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) &
      C % IGNORABILITY  =  IgnorabilityOption

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart'

    C % Name  =  'Chart'
    if ( present ( NameOption ) ) &
      C % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    call SetDimensionality &
           ( C, nDimensionsOption, iDimensionalityOption )
    call SetCoordinateSystem &
           ( C, CoordinateLabelOption, CoordinateSystemOption, &
             CoordinateUnitOption )

    allocate ( C % Connectivity )
    associate ( Cy  =>  C % Connectivity )
    call Cy % Initialize ( C % nDimensions )
    end associate !-- Cy

  end subroutine Initialize_H


  subroutine Show_C ( C )

    class ( Chart_H_Form ), intent ( in ) :: &
      C

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    associate ( nD  =>  C % nDimensions )

    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )

    call Show ( C % CoordinateSystem, 'CoordinateSystem', C % IGNORABILITY )
    call Show ( C % CoordinateLabel ( : nD ), 'CoordinateLabel', &
                C % IGNORABILITY )
    call Show ( C % CoordinateUnit ( : nD ), 'CoordinateUnit', &
                C % IGNORABILITY )

    call Show ( C % Periodic, 'Periodic', C % IGNORABILITY )

    call C % Connectivity % Show ( C % IGNORABILITY + 1 )

    end associate !-- nD

  end subroutine Show_C


  impure elemental subroutine Finalize_C ( C )

    type ( Chart_H_Form ), intent ( inout ) :: &
      C

    if ( allocated ( C % Connectivity ) ) &
      deallocate ( C % Connectivity )

    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

  end subroutine Finalize_C


  impure elemental subroutine Finalize_E ( CE )
    
    type ( ChartElement ), intent ( inout ) :: &
      CE

    if ( allocated ( CE % Element ) ) &
      deallocate ( CE % Element )

  end subroutine Finalize_E


  subroutine SetDimensionality ( C, nDimensionsOption, iDimensionalityOption )

    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    integer ( KDI ) :: &
      iDimensionality
    character ( LDL ), dimension ( : ), allocatable :: &
      Dimensionality

    if ( present ( nDimensionsOption ) ) then
      C % nDimensions = nDimensionsOption
      call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )
      return
    end if

    !-- Allow for specification of base manifold and bundle 
    !   dimensionalities; take the first element here, the dimensionality
    !   of the base manifold
    call Split ( PROGRAM_HEADER % Dimensionality, '_', Dimensionality )

    iDimensionality = 1
    if ( present ( iDimensionalityOption ) ) &
      iDimensionality = iDimensionalityOption

    if ( iDimensionality > size ( Dimensionality ) ) then
      call Show ( 'Too few dimensionalities specified', CONSOLE % ERROR )
      call Show ( 'Chart_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    select case ( trim ( Dimensionality ( iDimensionality ) ) )
    case ( '1D' )
      C % nDimensions = 1
    case ( '2D' )
      C % nDimensions = 2
    case ( '3D' )
      C % nDimensions = 3
    case default
      call Show ( 'PROGRAM_HEADER % Dimensionality not recognized', &
                  CONSOLE % WARNING )
      call Show ( 'Chart_H__Form', 'module', CONSOLE % WARNING )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % WARNING )
      call Show ( 'Defaulting to 3D', CONSOLE % WARNING )
      C % nDimensions = 3
    end select !-- Dimensionality

  end subroutine SetDimensionality


  subroutine SetCoordinateSystem &
               ( C, CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption )

    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption

    associate ( nD => C % nDimensions )

    C % CoordinateSystem = 'RECTANGULAR'
    if ( present ( CoordinateSystemOption ) ) &
      C % CoordinateSystem = CoordinateSystemOption
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateSystem, 'CoordinateSystem' )

    C % CoordinateLabel  =  [ 'X', 'Y', 'Z' ]
    select case ( trim ( C % CoordinateSystem ) )
    case ( 'CYLINDRICAL' )
      C % CoordinateLabel  =  [ 'R_Perp', 'Z     ', 'Phi   ' ]
    case ( 'SPHERICAL' )
      C % CoordinateLabel  =  [ 'R    ', 'Theta', 'Phi  ' ]
    end select !-- CoordinateSystem
    if ( present ( CoordinateLabelOption ) ) &
      C % CoordinateLabel ( : nD ) = CoordinateLabelOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateLabel ( : nD ), 'CoordinateLabel' )

    C % CoordinateUnit = [ UNIT % IDENTITY, UNIT % IDENTITY, UNIT % IDENTITY ]
    if ( present ( CoordinateUnitOption ) ) &
      C % CoordinateUnit ( : nD ) = CoordinateUnitOption ( : nD )

    end associate !-- nD

  end subroutine SetCoordinateSystem


end module Chart_H__Form
