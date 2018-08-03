!-- Chart_SLD_CC represents a distributed single-level chart with a central
!   core, using spherical coordinates with proportional spacing.

module Chart_SLD_CC__Form

  !-- Chart_SingleLevelDistributed_CentralCore_Form

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ) :: Chart_SLD_CC_Form
    integer ( KDI ) :: &
      nCellsPolar, &
      nCellsCore
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax, &
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
  contains
    procedure, private, pass :: &
      Initialize_CC
    generic, public :: &
      Initialize => Initialize_CC
    procedure, private, pass :: &
      ShowHeader
    procedure, public, pass :: &
      SetCoarsening
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

    integer ( KDI ) :: &
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    CoordinateSystem = 'SPHERICAL'

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    IsPeriodic = .false.
    IsPeriodic ( 3 ) = .true.

    C % RadiusMax = 10.0_KDR
    if ( present ( RadiusMaxOption ) ) &
      C % RadiusMax = RadiusMaxOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusMax, 'RadiusMax' )

    associate ( Pi => CONSTANT % PI )
    MinCoordinate = [   0.0_KDR    , 0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]
    end associate !-- Pi

    C % RadiusCore = C % RadiusMax / 8.0_KDR
    if ( present ( RadiusCoreOption ) ) &
      C % RadiusCore = RadiusCoreOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusCore, 'RadiusCore' )

    C % nCellsPolar = 128
    if ( present ( nCellsPolarOption ) ) &
      C % nCellsPolar = nCellsPolarOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsPolar, 'nCellsPolar' )

    if ( .not. any ( C % nCellsPolar &
                       == [ 32, 64, 128, 256, 512, 1024, 2048, 4096 ] ) ) &
    then 
      call Show ( 'nCellsPolar must be a power of 2', CONSOLE % ERROR )
      call Show ( 'Initialize_CC', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Chart_SLD_CC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    C % nCellsCore  =  10 * ( C % nCellsPolar / 32 )  

    C % RadialRatio = 1
    if ( present ( RadialRatioOption ) ) &
      C % RadialRatio = RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    nCellsRadial    = C % RadialRatio  *  C % nCellsPolar !-- Aim for RadiusMax
    nCellsPolar     = C % nCellsPolar
    nCellsAzimuthal = 2 * nCellsPolar
 
    nCells = [ nCellsRadial, 1, 1 ]
    if ( Atlas % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( Atlas % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  CONSTANT % PI / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusCore

    call C % Chart_SLD_Form % Initialize &
           ( Atlas, IsPeriodic, iChart, SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption, &
             nEqualOption = C % nCellsCore )

  end subroutine Initialize_CC


  subroutine ShowHeader ( C )

    class ( Chart_SLD_CC_Form ), intent ( in ) :: &
      C

    call C % Chart_SLD_Form % Show ( )

    call Show ( 'Chart_SLD_CC parameters' )
    call Show ( C % nCellsPolar, 'nCellsPolar' )
    call Show ( C % nCellsCore, 'nCellsCore' )
    call Show ( C % RadiusCore, C % CoordinateUnit ( 1 ), 'RadiusCore' )
    call Show ( C % RadiusCore / C % nCellsCore, C % CoordinateUnit ( 1 ), &
                'CellWidthCore' )
    call Show ( C % RadialRatio, 'RadialRatio' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), &
                'RadiusMax requested' )
    call Show ( C % MaxCoordinate ( 1 ), C % CoordinateUnit ( 1 ), &
                'RadiusMax actual' )

  end subroutine ShowHeader


  subroutine SetCoarsening ( C )

    class ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iV  !-- iValue
    class ( GeometryFlatForm ), pointer :: &
      G

    G => C % Geometry ( )

    if ( C % nDimensions < 2 ) &
      return

    associate &
      (      R_MW => C % RadiusCore, &
           dTheta => C % WidthLeft ( 2 ) % Value ( 1 ) &
                     +  C % WidthRight ( 2 ) % Value ( 1 ), &
                 R => G % Value ( :, G % CENTER_U ( 1 ) ), &
             Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
            Crsn_2 => G % Value ( :, G % COARSENING ( 2 ) ) )
      
    C % MinWidth  =  R_MW * dTheta
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth' )

    do iV = 1, size ( R )
      if ( .not. C % IsProperCell ( iV ) ) &
        cycle
      Crsn_2 ( iV )  =  1.0_KDR
      Coarsen_2: do
        if ( Crsn_2 ( iV )  *  R ( iV )  *  dTheta  >  C % MinWidth ) &
          exit Coarsen_2
        Crsn_2 ( iV )  =  2.0_KDR  *  Crsn_2 ( iV )
      end do Coarsen_2
    end do !-- iV

    if ( C % nDimensions < 3 ) &
      return

    associate &
      (   dPhi => C % WidthLeft ( 3 ) % Value ( 1 ) &
                  +  C % WidthRight ( 3 ) % Value ( 1 ), &
        Crsn_3 => G % Value ( :, G % COARSENING ( 3 ) ) )

    do iV = 1, size ( R )
      if ( .not. C % IsProperCell ( iV ) ) &
        cycle
      Crsn_3 ( iV )  =  1.0_KDR
      Coarsen_3: do
        if ( Crsn_3 ( iV )  *  R ( iV )  *  sin ( Theta ( iV ) )  *  dPhi  &
             >  C % MinWidth ) &
          exit Coarsen_3
        Crsn_3 ( iV )  =  2.0_KDR  *  Crsn_3 ( iV )
      end do Coarsen_3
    end do !-- iV

    nullify ( G )

    end associate !-- dPhi, etc.
    end associate !-- R_MW, etc.

  end subroutine SetCoarsening


end module Chart_SLD_CC__Form
