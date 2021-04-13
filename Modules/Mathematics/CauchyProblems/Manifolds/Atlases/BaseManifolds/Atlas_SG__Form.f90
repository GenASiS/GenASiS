module Atlas_SG__Form

  !-- Atlas_SingleGrid__Form

  use Basics
  use Charts
  use Atlas_H__Form

  implicit none
  private

  type, public, extends ( Atlas_H_Form ) :: Atlas_SG_Form
    class ( Grid_S_Form ), pointer :: &
      Grid => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SG
    generic, public :: &
      Initialize => Initialize_SG
    final :: &
      Finalize
  end type Atlas_SG_Form


contains


  subroutine Initialize_SG &
               ( A, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, PeriodicOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Atlas_SG_Form ), intent ( inout ), target :: &
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
      PeriodicOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
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

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SG'

    call A % Initialize_H &
           ( NameOption = NameOption, &
             IgnorabilityOption = IgnorabilityOption )

    allocate ( Grid_S_Form :: A % Chart ( 1 ) % Element )
    select type ( G  =>  A % Chart ( 1 ) % Element )
    class is ( Grid_S_Form )

    call G % Initialize &
           ( CommunicatorOption, SpacingOption, CoordinateLabelOption, &
             CoordinateSystemOption, NameOption, PeriodicOption, &
             CoordinateUnitOption, MinCoordinateOption, &
             MaxCoordinateOption, RatioOption, ScaleOption, &
             nCellsOption, nGhostLayersOption, nBricksOption, &
             nBricksCompatibleOption, IgnorabilityOption, &
             nDimensionsOption, nEqualOption, iDimensionalityOption )

    A % Grid  =>  G

    end select !--  G

  end subroutine Initialize_SG


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SG_Form ), intent ( inout ) :: &
      A

    nullify ( A % Grid )

  end subroutine Finalize


end module Atlas_SG__Form
