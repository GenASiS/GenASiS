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
    type ( CommunicatorForm ), allocatable :: &
      Communicator_2, &
      Communicator_3
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

    C % RadiusCore = C % RadiusMax / 8.0_KDR
    if ( present ( RadiusCoreOption ) ) &
      C % RadiusCore = RadiusCoreOption
    call PROGRAM_HEADER % GetParameter ( C % RadiusCore, 'RadiusCore' )

    call C % InitializeTemplate_C &
           ( Atlas = Atlas, &
             RadiusMin = 0.0_KDR, &
             RadiusScale = C % RadiusCore, &
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
      iC, jC, kC, &  !-- iCell, etc.
      iR, &          !-- iRank
      iRB, iRP, &    !-- iRankBrick, iRankPillar
      iP, &         !-- iPillar
      nRanks_2, &
      nRanks_3, &
      nPillarsTotal, &
      nPillarsSurplus
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank_2, Rank_3, &
      nPillars_2, nPillars_3
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      nPillarsTransverse_2, &
      nPillarsTransverse_3
    class ( GeometryFlatForm ), pointer :: &
      G

    G => C % Geometry ( )

    associate &
      (   nB => C % nBricks, &
         nCB => C % nCellsBrick, &
        R_MW => C % RadiusCore )
             
    !-- Polar coarsening

    if ( C % nDimensions < 2 ) &
      return

    call Show ( 'SetCoarsening_2', C % IGNORABILITY )
    call Show ( C % Name, 'Chart', C % IGNORABILITY )

    associate &
      (    R_G => C % Center ( 1 ) % Value, &  !-- R_Global
        dTheta => C % WidthLeft ( 2 ) % Value ( 1 ) &
                  +  C % WidthRight ( 2 ) % Value ( 1 ), &
             R => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Crsn_2 => G % Value ( :, G % COARSENING ( 2 ) ) )

    !-- Determin MinWidth
    C % MinWidth  =  R_MW * dTheta
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth', &
                C % IGNORABILITY )

    !-- Set coarsening factor
    do iV = 1, size ( R )
      if ( .not. C % IsProperCell ( iV ) ) &
        cycle
      Crsn_2 ( iV )  =  1.0_KDR
      Coarsen_2: do
        if ( Crsn_2 ( iV )  *  R ( iV )  *  dTheta  >  C % MinWidth ) &
          exit Coarsen_2
        Crsn_2 ( iV )  =  2.0_KDR  *  Crsn_2 ( iV )
      end do Coarsen_2
      Crsn_2 ( iV )  =  min ( Crsn_2 ( iV ), real ( C % nCells ( 2 ) ) )
    end do !-- iV

    !-- Count coarsening pillars in transverse brick space 
    allocate ( nPillarsTransverse_2 ( nB ( 1 ), nB ( 3 ) ) )
    nPillarsTransverse_2 = 0
    do kB = 1, nB ( 3 )
      do iB = 1, nB ( 1 )
        do kC = ( kB - 1 ) * nCB ( 3 ) + 1, kB * nCB ( 3 )
          do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
            if ( R_G ( iC ) * dTheta  <  C % MinWidth ) &
              nPillarsTransverse_2 ( iB, kB )  &
                =  nPillarsTransverse_2 ( iB, kB )  +  1
          end do !-- iC
        end do !-- kC
      end do !-- iB
    end do !-- kB
    call Show ( nPillarsTransverse_2, 'nPillarsTransverse_2', &
                C % IGNORABILITY + 1 )

    !-- Determine pillar distribution
    nRanks_2 = count ( nPillarsTransverse_2 > 0 ) * nB ( 2 )
    allocate ( nPillars_2 ( 0 : nRanks_2 - 1 ) )
    nPillarsTotal = sum ( nPillarsTransverse_2 )
    nPillarsSurplus = mod ( nPillarsTotal, nRanks_2 )
    nPillars_2 = nPillarsTotal / nRanks_2
    nPillars_2 ( 0 : nPillarsSurplus - 1 ) &
      = nPillars_2 ( 0 : nPillarsSurplus - 1 ) + 1
    call Show ( nPillars_2, 'nPillars_2', C % IGNORABILITY + 1 )

    !-- Determine subcommunicator ranks
    allocate ( Rank_2 ( 0 : nRanks_2 - 1 ) )
    iR = 0
    do kB = 1, nB ( 3 )
      iLoop: do iB = 1, nB ( 1 )
        if ( nPillarsTransverse_2 ( iB, kB ) == 0 ) &
          cycle iLoop
        do jB = 1, nB ( 2 )
          Rank_2 ( iR ) = ( kB - 1 ) * nB ( 1 ) * nB ( 2 )  &
                          +  ( jB - 1 ) * nB ( 1 )  &
                          +  ( iB - 1 )
          iR = iR + 1
        end do !-- jB
      end do iLoop !-- iB
    end do !-- kB
    call Show ( Rank_2, 'Rank_2', C % IGNORABILITY + 1 )

    !-- Create subcommunicator
    allocate ( C % Communicator_2 )
    associate ( Comm_2 => C % Communicator_2 )
    call Comm_2 % Initialize &
           ( C % Atlas % Communicator, Rank_2, NameOption = 'Coarsening_2' )

    if ( Comm_2 % Initialized ) then

      C % nPillars_2 = nPillars_2 ( Comm_2 % Rank )

      !-- Determine nSegmentsFrom and nSegmentsTo
      allocate ( C % nSegmentsFrom_2 ( 0 : nRanks_2 - 1 ) )
      allocate ( C % nSegmentsTo_2 ( 0 : nRanks_2 - 1 ) )
      iRB = 0
      iRP = 0
      iP = 0
      C % nSegmentsFrom_2 = 0
      C % nSegmentsTo_2   = 0
      do kB = 1, nB ( 3 )
        iLoop: do iB = 1, nB ( 1 )
          if ( nPillarsTransverse_2 ( iB, kB ) == 0 ) &
            cycle iLoop
          do kC = ( kB - 1 ) * nCB ( 3 ) + 1, kB * nCB ( 3 )
            do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
              if ( R_G ( iC ) * dTheta  <  C % MinWidth ) then
                iP = iP + 1
                do jB = 0, nB ( 2 ) - 1
                  if ( Comm_2 % Rank == iRP ) &
                    C % nSegmentsFrom_2 ( iRB + jB )  &
                      =  C % nSegmentsFrom_2 ( iRB + jB )  +  1
                  if ( Comm_2 % Rank == iRB + jB ) &
                    C % nSegmentsTo_2 ( iRP ) &
                      =  C % nSegmentsTo_2 ( iRP )  +  1
                end do !-- jB
                if ( iP == nPillars_2 ( iRP ) ) then
                  iRP = iRP + 1
                  iP = 0
                end if
              end if
            end do !-- iC
          end do !-- kC
          iRB = iRB + nB ( 2 )
        end do iLoop !-- iB
      end do !-- kB
      call Show ( 'Communication counting', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsFrom_2, 'nSegmentsFrom_2', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsTo_2, 'nSegmentsTo_2', C % IGNORABILITY + 1 )

    end if !-- Comm_2 % Initialized
    end associate !-- Comm_2


    !-- Azimuthal coarsening

    if ( C % nDimensions < 3 ) &
      return

    call Show ( 'SetCoarsening_3', C % IGNORABILITY )
    call Show ( C % Name, 'Chart', C % IGNORABILITY )

    associate &
      ( Theta_G => C % Center ( 2 ) % Value, &  !-- Theta_Global
           dPhi => C % WidthLeft ( 3 ) % Value ( 1 ) &
                   +  C % WidthRight ( 3 ) % Value ( 1 ), &
          Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
         Crsn_3 => G % Value ( :, G % COARSENING ( 3 ) ) )

    !-- Set coarsening factor
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
      Crsn_3 ( iV )  =  min ( Crsn_3 ( iV ), real ( C % nCells ( 3 ) ) )
    end do !-- iV

    !-- Count coarsening pillars in transverse brick space 
    allocate ( nPillarsTransverse_3 ( nB ( 1 ), nB ( 2 ) ) )
    nPillarsTransverse_3 = 0
    do jB = 1, nB ( 2 )
      do iB = 1, nB ( 1 )
        do jC = ( jB - 1 ) * nCB ( 2 ) + 1, jB * nCB ( 2 )
          do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
            if ( R_G ( iC ) * sin ( Theta_G ( jC ) ) * dPhi  &
                 <  C % MinWidth ) &
              nPillarsTransverse_3 ( iB, jB )  &
                =  nPillarsTransverse_3 ( iB, jB )  +  1
          end do !-- iC
        end do !-- jC
      end do !-- iB
    end do !-- jB
    call Show ( nPillarsTransverse_3, 'nPillarsTransverse_3', &
                C % IGNORABILITY + 1 )

    !-- Determine pillar distribution
    nRanks_3 = count ( nPillarsTransverse_3 > 0 ) * nB ( 3 )
    allocate ( nPillars_3 ( 0 : nRanks_3 - 1 ) )
    nPillarsTotal = sum ( nPillarsTransverse_3 )
    nPillarsSurplus = mod ( nPillarsTotal, nRanks_3 )
    nPillars_3 = nPillarsTotal / nRanks_3
    nPillars_3 ( 0 : nPillarsSurplus - 1 ) &
      = nPillars_3 ( 0 : nPillarsSurplus - 1 ) + 1
    call Show ( nPillars_3, 'nPillars_3', C % IGNORABILITY + 1 )

    !-- Determine subcommunicator ranks
    allocate ( Rank_3 ( 0 : nRanks_3 - 1 ) )
    iR = 0
    do jB = 1, nB ( 2 )
      iLoop: do iB = 1, nB ( 1 )
        if ( nPillarsTransverse_3 ( iB, jB ) == 0 ) &
          cycle iLoop
        do kB = 1, nB ( 3 )
          Rank_3 ( iR ) = ( kB - 1 ) * nB ( 1 ) * nB ( 2 )  &
                          +  ( jB - 1 ) * nB ( 1 )  &
                          +  ( iB - 1 ) 
          iR = iR + 1
        end do !-- kB
      end do iLoop !-- iB
    end do !-- jB
    call Show ( Rank_3, 'Rank_3', C % IGNORABILITY + 1 )

    !-- Create subcommunicator
    allocate ( C % Communicator_3 )
    associate ( Comm_3 => C % Communicator_3 )
    call Comm_3 % Initialize &
           ( C % Atlas % Communicator, Rank_3, NameOption = 'Coarsening_3' )

    if ( Comm_3 % Initialized ) then

      C % nPillars_3 = nPillars_3 ( Comm_3 % Rank )

      !-- Determine nSegmentsFrom and nSegmentsTo
      allocate ( C % nSegmentsFrom_3 ( 0 : nRanks_3 - 1 ) )
      allocate ( C % nSegmentsTo_3 ( 0 : nRanks_3 - 1 ) )
      iRB = 0
      iRP = 0
      iP = 0
      C % nSegmentsFrom_3 = 0
      C % nSegmentsTo_3   = 0
      do jB = 1, nB ( 2 )
        iLoop: do iB = 1, nB ( 1 )
          if ( nPillarsTransverse_3 ( iB, jB ) == 0 ) &
            cycle iLoop
          do jC = ( jB - 1 ) * nCB ( 2 ) + 1, jB * nCB ( 2 )
            do iC = ( iB - 1 ) * nCB ( 1 ) + 1, iB * nCB ( 1 )
              if ( R_G ( iC ) * sin ( Theta_G ( jC ) ) * dPhi  &
                   <  C % MinWidth ) &
              then
                iP = iP + 1
                do kB = 0, nB ( 3 ) - 1
                  if ( Comm_3 % Rank == iRP ) &
                    C % nSegmentsFrom_3 ( iRB + kB )  &
                      =  C % nSegmentsFrom_3 ( iRB + kB )  +  1
                  if ( Comm_3 % Rank == iRB + kB ) &
                    C % nSegmentsTo_3 ( iRP ) &
                      =  C % nSegmentsTo_3 ( iRP )  +  1
                end do !-- kB
                if ( iP == nPillars_3 ( iRP ) ) then
                  iRP = iRP + 1
                  iP = 0
                end if
              end if
            end do !-- iC
          end do !-- jC
          iRB = iRB + nB ( 3 )
        end do iLoop !-- iB
      end do !-- jB
      call Show ( 'Communication counting', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsFrom_3, 'C % nSegmentsFrom_3', &
                  C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsTo_3, 'C % nSegmentsTo_3', &
                  C % IGNORABILITY + 1 )

    end if !-- Comm_3 % Initialized
    end associate !-- Comm_3


    end associate !-- Theta_G, etc.
    end associate !-- R_G, etc.
    end associate !-- nB, etc.

    nullify ( G )

!call C % Atlas % Communicator % Synchronize ( )
!call PROGRAM_HEADER % Abort ( )
  end subroutine SetCoarsening


  impure elemental subroutine Finalize ( C )

    type ( Chart_SLD_CC_Form ), intent ( inout ) :: &
      C

    if ( allocated ( C % Communicator_3 ) ) &
      deallocate ( C % Communicator_3 )
    if ( allocated ( C % Communicator_2 ) ) &
      deallocate ( C % Communicator_2 )

    if ( allocated ( C % nSegmentsTo_3 ) ) &
      deallocate ( C % nSegmentsTo_3 )
    if ( allocated ( C % nSegmentsFrom_3 ) ) &
      deallocate ( C % nSegmentsFrom_3 )
    if ( allocated ( C % nSegmentsTo_2 ) ) &
      deallocate ( C % nSegmentsTo_2 )
    if ( allocated ( C % nSegmentsFrom_2 ) ) &
      deallocate ( C % nSegmentsFrom_2 )

  end subroutine Finalize


end module Chart_SLD_CC__Form
