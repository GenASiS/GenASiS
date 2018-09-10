!-- Chart_SLD_CC represents a distributed single-level chart with a central
!   core, using spherical coordinates with proportional spacing.

module Chart_SLD_CC__Form

  !-- Chart_SingleLevelDistributed_CentralCore_Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD_C__Template

  implicit none
  private

  type, public, extends ( Chart_SLD_C_Template ) :: Chart_SLD_CC_Form
    integer ( KDI ) :: &
      nCellsCore
    real ( KDR ) :: &
      RadiusCore
  contains
    procedure, private, pass :: &
      Initialize_CC
    generic, public :: &
      Initialize => Initialize_CC
    procedure, private, pass :: &
      SetCore
    procedure, private, pass :: &
      ShowHeader
    procedure, public, pass :: &
      SetCoarsening
    final :: &
      Finalize
  end type Chart_SLD_CC_Form

contains


  subroutine Initialize_CC &
               ( C, Atlas, iChart, CoordinateUnitOption, RadiusMaxOption, &
                 RadiusCoreOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption )

    class ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    C % RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      C % RadiusMax = RadiusMaxOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMax, 'RadiusMax' )

    C % RadiusCore = C % RadiusMax / 8.0_KDR
    if ( present ( RadiusCoreOption ) ) &
      C % RadiusCore = RadiusCoreOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusCore, 'RadiusCore' )

    C % RadiusScale = C % RadiusCore

    call C % InitializeTemplate_C &
           ( Atlas = Atlas, &
             RadiusMin = 0.0_KDR, &
             iChart = iChart, &
             CoordinateUnitOption = CoordinateUnitOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadialRatioOption = RadialRatioOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nEqualOption = C % nCellsCore )

  end subroutine Initialize_CC


  subroutine SetCore ( C )

    class ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    if ( .not. any ( C % nCellsPolar &
                       == [ 32, 64, 128, 256, 512, 1024, 2048, 4096 ] ) ) &
    then 
      call Show ( 'nCellsPolar must be a power of 2 between 32 and 4096', &
                  CONSOLE % ERROR )
      call Show ( 'SetCore', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Chart_SLD_CC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    C % nCellsCore  =  10 * ( C % nCellsPolar / 32 )  

  end subroutine SetCore


  subroutine ShowHeader ( C )

    class ( Chart_SLD_CC_Form ), intent ( in ) :: &
      C

    call C % Chart_SLD_Form % Show ( )
    call C % ShowHeaderTemplate_C ( )

    call Show ( 'Chart_SLD_CC parameters' )
    call Show ( C % RadiusCore, C % CoordinateUnit ( 1 ), 'RadiusCore' )
    call Show ( C % nCellsCore, 'nCellsCore' )

  end subroutine ShowHeader


  subroutine SetCoarsening ( C )

    class ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    call C % SetCoarseningPolar ( )
    call C % SetCoarseningAzimuthal ( )

  end subroutine SetCoarsening


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    call C % FinalizeTemplate_C ( )

  end subroutine Finalize


end module Chart_SLD_CC__Form
