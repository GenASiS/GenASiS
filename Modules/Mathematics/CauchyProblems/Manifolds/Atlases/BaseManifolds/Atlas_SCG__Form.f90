module Atlas_SCG__Form

  !-- Atlas_SingleChartGrid__Form

  use Basics
  use Charts
  use Atlas_H__Form

  implicit none
  private

  type, public, extends ( Atlas_H_Form ) :: Atlas_SCG_Form
    class ( Chart_GS_Form ), pointer :: &
      Chart_GS => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG
    generic, public :: &
      Initialize => Initialize_SCG
    final :: &
      Finalize
  end type Atlas_SCG_Form


contains


  subroutine Initialize_SCG &
               ( A, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Atlas_SCG_Form ), intent ( inout ), target :: &
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

    logical :: &
      PreviouslyAllocated

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SCG'

    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( 1 ) )
    
    if ( allocated ( A % Chart ( 1 ) % Element ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( Chart_GS_Form :: A % Chart ( 1 ) % Element )
    end if

    call A % Initialize_H &
           ( NameOption = NameOption, &
             IgnorabilityOption = IgnorabilityOption, &
             nChartsOption = 1 )

    if ( .not. PreviouslyAllocated ) then

      select type ( C  =>  A % Chart ( 1 ) % Element )
      class is ( Chart_GS_Form )

      call C % Initialize &
             ( CommunicatorOption, SpacingOption, CoordinateLabelOption, &
               CoordinateSystemOption, NameOption, CoordinateUnitOption, &
               MinCoordinateOption, MaxCoordinateOption, RatioOption, &
               ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
               nBricksCompatibleOption, IgnorabilityOption, &
               nDimensionsOption, nEqualOption, iDimensionalityOption )

      end select !--  C

    end if !-- PreviouslyAllocated  

    select type ( C  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_Form )
      A % Chart_GS  =>  C
    end select !-- C
      
  end subroutine Initialize_SCG


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS )

  end subroutine Finalize


end module Atlas_SCG__Form
