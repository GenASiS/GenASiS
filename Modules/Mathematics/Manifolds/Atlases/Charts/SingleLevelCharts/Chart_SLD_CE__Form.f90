!-- Chart_SLD_CE represents a distributed single-level chart with an inner
!   boundary at finite radius, using spherical coordinates with proportional
!   spacing.

module Chart_SLD_CE__Form

  !-- Chart_SingleLevelDistributed_CentralExcision_Form

  use Basics
  use AtlasBasics
  use Chart_SLD_C__Template

  implicit none
  private

  type, public, extends ( Chart_SLD_C_Template ) :: Chart_SLD_CE_Form
    real ( KDR ) :: &
      RadiusExcision
  contains
    procedure, private, pass :: &
      Initialize_CE
    generic, public :: &
      Initialize => Initialize_CE
    procedure, private, pass :: &
      ShowHeader
    final :: &
      Finalize
  end type Chart_SLD_CE_Form

contains


  subroutine Initialize_CE &
               ( C, Atlas, iChart, CoordinateUnitOption, RadiusMaxOption, &
                 RadiusExcisionOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption )

    class ( Chart_SLD_CE_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    C % RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      C % RadiusMax = RadiusMaxOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMax, 'RadiusMax' )

    C % RadiusExcision = 1.0_KDR
    if ( present ( RadiusExcisionOption ) ) &
      C % RadiusExcision = RadiusExcisionOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusExcision, 'RadiusExcision' )

    C % RadiusScale = C % RadiusExcision

    call C % InitializeTemplate_C &
           ( Atlas = Atlas, &
             RadiusMin = C % RadiusExcision, &
             iChart = iChart, &
             CoordinateUnitOption = CoordinateUnitOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadialRatioOption = RadialRatioOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nCellsPolarOption = nCellsPolarOption )

  end subroutine Initialize_CE


  subroutine ShowHeader ( C )

    class ( Chart_SLD_CE_Form ), intent ( in ) :: &
      C

    call C % Chart_SLD_Form % Show ( )
    call C % ShowHeaderTemplate_C ( )

    call Show ( 'Chart_SLD_CE parameters' )
    call Show ( C % RadiusExcision, C % CoordinateUnit ( 1 ), 'RadiusExcision' )

  end subroutine ShowHeader


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLD_CE_Form ), intent ( inout ) :: &
      C

    call C % FinalizeTemplate_C ( )

  end subroutine Finalize


end module Chart_SLD_CE__Form
