module Atlas_SCG_SC__Form

  !-- Atlas_SingleChartGrid_SymmetricCurvilinear_Form

  use Basics
  use Charts
  use Atlas_SCG__Form

  implicit none
  private

  type, public, extends ( Atlas_SCG_Form ) :: Atlas_SCG_SC_Form
    class ( Chart_GS_SC_Form ), pointer :: &
      Chart_GS_SC => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG
    procedure, private, pass :: &
      Initialize_SCG_SC
    generic, public :: &
      Initialize => Initialize_SCG_SC
    final :: &
      Finalize
  end type Atlas_SCG_SC_Form


contains


  subroutine Initialize_SCG &
               ( A, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, EvenDecompositionOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, nCellsOption, &
                 nGhostLayersOption, nBricksOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Atlas_SCG_SC_Form ), intent ( inout ), target :: &
      A
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption, &
      NameOption
    logical ( KDL ), dimension ( : ), intent ( in ), optional :: &
      EvenDecompositionOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption, &
      nBricksOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nDimensionsOption, &
      nEqualOption, &
      iDimensionalityOption

      call Show ( 'The method Initialize_SCG is not appropriate for ' &
                  // 'this class.', CONSOLE % ERROR )
      call Show ( 'Please use a different Initialize interface.', &
                  CONSOLE % ERROR )
      call Show ( 'Atlas_SCG_SC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_SCG', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )

  end subroutine Initialize_SCG


  subroutine Initialize_SCG_SC &
               ( A, RadiusMax, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, nGhostLayersOption, nCellsRadiusOption, &
                 nDimensionsOption )

    class ( Atlas_SCG_SC_Form ), intent ( inout ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      RadiusMax
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsRadiusOption, &
      nDimensionsOption

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SCG_SC'

    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( 1 ) )
    
    allocate ( Chart_GS_SC_Form :: A % Chart ( 1 ) % Element )

    call A % Atlas_SCG_Form % Initialize ( NameOption = NameOption )

    select type ( C  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_SC_Form )

      call C % Initialize &
             ( RadiusMax, &
               CommunicatorOption = CommunicatorOption, &
               NameOption = NameOption, &
               CoordinateUnitOption = CoordinateUnitOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsRadiusOption = nCellsRadiusOption, &
               nDimensionsOption = nDimensionsOption )

      A % Chart_GS_SC  =>  C

    end select !--  C

  end subroutine Initialize_SCG_SC


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_SC_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS_SC )

  end subroutine Finalize


end module Atlas_SCG_SC__Form
