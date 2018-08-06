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
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      nPillars_2, &
      nPillars_3
    real ( KDR ) :: &
      RadiusCore, &
      RadiusMax, &
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
    type ( CommunicatorForm ), allocatable :: &
      Communicator_2, &
      Communicator_3
  contains
    procedure, private, pass :: &
      Initialize_CC
    generic, public :: &
      Initialize => Initialize_CC
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
      iV, &          !-- iValue
      iB, jB, kB, &  !-- iBrick, etc.
      iC, jC, kC     !-- iCell, etc.
    class ( GeometryFlatForm ), pointer :: &
      G

    G => C % Geometry ( )

    if ( C % nDimensions < 2 ) &
      return

    call Show ( 'SetCoarsening', C % IGNORABILITY )
    call Show ( C % Name, 'Chart', C % IGNORABILITY )

    associate &
      (       nCB => C % nCellsBrick, &
             R_MW => C % RadiusCore, &
              R_G => C % Center ( 1 ) % Value, &  !-- R_Global
           dTheta => C % WidthLeft ( 2 ) % Value ( 1 ) &
                     +  C % WidthRight ( 2 ) % Value ( 1 ), &
                 R => G % Value ( :, G % CENTER_U ( 1 ) ), &
            Crsn_2 => G % Value ( :, G % COARSENING ( 2 ) ) )
      
    C % MinWidth  =  R_MW * dTheta
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth', &
                C % IGNORABILITY )

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

    allocate ( C % nPillars_2 ( C % nBricks ( 1 ), C % nBricks ( 3 ) ) )
    C % nPillars_2 = 0
    do kB = 1, C % nBricks ( 3 )
      do iB = 1, C % nBricks ( 1 )
        do kC = ( kB - 1 ) * nCB ( 3 ) + 1, kB * nCB ( 3 )
          do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
            if ( R_G ( iC ) * dTheta  <  C % MinWidth ) &
              C % nPillars_2 ( iB, kB )  =  C % nPillars_2 ( iB, kB )  +  1
          end do !-- iC
        end do !-- kC
      end do !-- iB
    end do !-- kB
    call Show ( C % nPillars_2, 'nPillars_2', C % IGNORABILITY )

    if ( C % nDimensions < 3 ) &
      return

    associate &
      ( Theta_G => C % Center ( 2 ) % Value, &  !-- Theta_Global
           dPhi => C % WidthLeft ( 3 ) % Value ( 1 ) &
                   +  C % WidthRight ( 3 ) % Value ( 1 ), &
          Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
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

    allocate ( C % nPillars_3 ( C % nBricks ( 1 ), C % nBricks ( 2 ) ) )
    C % nPillars_3 = 0
    do jB = 1, C % nBricks ( 2 )
      do iB = 1, C % nBricks ( 1 )
        do jC = ( jB - 1 ) * nCB ( 2 ) + 1, jB * nCB ( 2 )
          do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
            if ( R_G ( iC ) * sin ( Theta_G ( jC ) ) * dPhi  <  C % MinWidth ) &
              C % nPillars_3 ( iB, jB )  =  C % nPillars_3 ( iB, jB )  +  1
          end do !-- iC
        end do !-- jC
      end do !-- iB
    end do !-- jB
    call Show ( C % nPillars_3, 'nPillars_3', C % IGNORABILITY )

    end associate !-- dPhi, etc.
    end associate !-- nCB, etc.

    nullify ( G )

  end subroutine SetCoarsening


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    if ( allocated ( C % Communicator_3 ) ) &
      deallocate ( C % Communicator_3 )
    if ( allocated ( C % Communicator_2 ) ) &
      deallocate ( C % Communicator_2 )
    if ( allocated ( C % nPillars_3 ) ) &
      deallocate ( C % nPillars_3 )
    if ( allocated ( C % nPillars_2 ) ) &
      deallocate ( C % nPillars_2 )

  end subroutine Finalize


end module Chart_SLD_CC__Form
