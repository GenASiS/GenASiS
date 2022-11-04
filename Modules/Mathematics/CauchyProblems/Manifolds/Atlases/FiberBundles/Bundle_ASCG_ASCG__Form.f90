module Bundle_ASCG_ASCG__Form

  !-- Bundle_AtlasSingleChartGrid_AtlasSingleChartGrid__Form

  use Basics
  use Charts
  use BaseManifolds
  use Bundle_H__Form

  implicit none
  private

  type, public, extends ( Bundle_H_Form ) :: Bundle_ASCG_ASCG_Form
    integer ( KDI ) :: &
      nBaseValues    = 0, &
      nFibers        = 0, &
      nSections      = 0
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )  !-- distribute Fibers
    class ( Chart_GS_Form ), pointer :: &
      Chart_GS_Base  => null ( ), &
      Chart_GS_Fiber => null ( )
    class ( Atlas_SCG_Form ), pointer :: &
      Atlas_SCG_Base  => null ( ), &
      Atlas_SCG_Fiber => null ( )
  contains
    procedure, private, pass :: &
      Initialize_ASCG_ASCG
    generic, public :: &
      Initialize => Initialize_ASCG_ASCG
    procedure, private, pass :: &
      Show_B
    final :: &
      Finalize
  end type Bundle_ASCG_ASCG_Form


contains


  subroutine Initialize_ASCG_ASCG &
               ( B, Base, CommunicatorOption, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, NameOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, nCellsOption, &
                 nGhostLayersOption, nBricksOption, nDimensionsOption )

    class ( Bundle_ASCG_ASCG_Form ), intent ( inout ), target :: &
      B
    class ( Atlas_SCG_Form ), intent ( in ), target :: &
      Base
    class ( CommunicatorForm ), intent ( in ), target, optional :: &
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
      nBricksOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    if ( B % Type  ==  '' ) &
      B % Type  =  'a Bundle_ASCG_ASCG'

    call B % Initialize_H ( Base, NameOption )

    !-- Base

    B % Base            =>  Base
    B % Atlas_SCG_Base  =>  Base
    B % Chart_GS_Base   =>  Base % Chart_GS

    associate ( CB  =>  B % Chart_GS_Base )
    B % nBaseValues  =  CB % nCellsLocal
    B % nFibers      =  CB % nCellsProper
    end associate !-- CB

    !-- Fibers

    allocate ( Atlas_SCG_Form :: B % Fiber )
    select type ( AF  =>  B % Fiber )
      class is ( Atlas_SCG_Form )

    call AF % Initialize &
           ( SpacingOption = SpacingOption, &
             CoordinateLabelOption = CoordinateLabelOption, &
             CoordinateSystemOption = CoordinateSystemOption, &
             NameOption = NameOption, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             RatioOption = RatioOption, &
             ScaleOption = ScaleOption, &
             nCellsOption = nCellsOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nBricksOption = nBricksOption, &
             nDimensionsOption = nDimensionsOption, &
             iDimensionalityOption = 2 )

    B % Atlas_SCG_Fiber  =>  AF
    B % Chart_GS_Fiber   =>  AF % Chart_GS

    associate ( CF  =>  B % Chart_GS_Fiber )
    B % nSections  =  CF % nCellsLocal
    end associate !-- CF

    if ( present ( CommunicatorOption ) ) then
      B % Communicator  =>  CommunicatorOption
    end if

    end select !-- AF

  end subroutine Initialize_ASCG_ASCG


  subroutine Show_B ( B )

    class ( Bundle_ASCG_ASCG_Form ), intent ( in ) :: &
      B

    call B % Bundle_H_Form % Show ( )

    call Show ( B % nBaseValues, 'nBaseValues', B % IGNORABILITY )
    call Show ( B % nFibers, 'nFibers', B % IGNORABILITY )
    call Show ( B % nSections, 'nSections', B % IGNORABILITY )

    if ( associated ( B % Communicator ) ) &
      call B % Communicator % Show ( B % IGNORABILITY )

  end subroutine Show_B


  subroutine Finalize ( B )

    type ( Bundle_ASCG_ASCG_Form ), intent ( inout ) :: &
      B
    
    nullify ( B % Chart_GS_Fiber )
    nullify ( B % Chart_GS_Base )
    nullify ( B % Atlas_SCG_Fiber )
    nullify ( B % Atlas_SCG_Base )
    nullify ( B % Communicator )

  end subroutine Finalize


end module Bundle_ASCG_ASCG__Form
