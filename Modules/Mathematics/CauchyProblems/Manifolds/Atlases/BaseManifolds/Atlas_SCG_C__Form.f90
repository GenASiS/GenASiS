module Atlas_SCG_C__Form

  !-- Atlas_SingleChartGrid_Central_Form

  use Basics
  use Charts
  use Atlas_SCG__Form

  implicit none
  private

  type, public, extends ( Atlas_SCG_Form ) :: Atlas_SCG_C_Form
    class ( Chart_GS_C_Form ), pointer :: &
      Chart_GS_C => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG
    final :: &
      Finalize
  end type Atlas_SCG_C_Form


contains


  subroutine Initialize_SCG &
               ( A, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Atlas_SCG_C_Form ), intent ( inout ), target :: &
      A
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption, &
      NameOption
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
      nBricksOption, &
      nBricksCompatibleOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption, &
      nDimensionsOption, &
      nEqualOption, &
      iDimensionalityOption

      call Show ( 'The method Initialize_SCG is not appropriate for ' &
                  // 'this class.', CONSOLE % ERROR )
      call Show ( 'Please use a different Initialize interface.', &
                  CONSOLE % ERROR )
      call Show ( 'Atlas_SCG_C_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_SCG', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )

  end subroutine Initialize_SCG


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_C_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS_C )

  end subroutine Finalize


end module Atlas_SCG_C__Form
