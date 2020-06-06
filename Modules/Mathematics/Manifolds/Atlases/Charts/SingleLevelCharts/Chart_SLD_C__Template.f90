!-- Chart_SLD_C represents a distributed single-level chart using spherical 
!   coordinates and proportional radial spacing.

module Chart_SLD_C__Template

  !-- Chart_SingleLevelDistributed_Central_Template

  use Basics
  use AtlasBasics
  use ChartBasics
  use Chart_SLD__Form

  implicit none
  private

  type, public, extends ( Chart_SLD_Form ), abstract :: Chart_SLD_C_Template
    integer ( KDI ) :: &
      nCellsPolar, &
      nPillars_2, &
      nPillars_3
    integer ( KDI ), dimension ( : ), allocatable :: &
      nSegmentsFrom_2, nSegmentsTo_2, &
      nSegmentsFrom_3, nSegmentsTo_3
    real ( KDR ) :: &
      RadiusMax, &
      RadiusScale, &
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
    type ( CommunicatorForm ), allocatable :: &
      Communicator_2, &
      Communicator_3
  contains
    procedure, public, pass :: &
      InitializeTemplate_C
    procedure, public, pass :: &
      ShowHeaderTemplate_C
    procedure, public, pass :: &
      SetCoarseningPolar
    procedure, public, pass :: &
      SetCoarseningAzimuthal
    procedure ( SC ), public, pass, deferred :: &
      SetCoarsening
    procedure, public, pass :: &
      FinalizeTemplate_C
    procedure, private, pass :: &
      SetCore
  end type Chart_SLD_C_Template

  abstract interface
    
    subroutine SC ( C )
      import Chart_SLD_C_Template
      class ( Chart_SLD_C_Template ), intent ( inout ) :: &
        C
    end subroutine SC

  end interface

contains


  subroutine InitializeTemplate_C &
               ( C, Atlas, RadiusMin, iChart, CoordinateUnitOption, &
                 RadiusMaxOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    real ( KDR ), intent ( in ) :: &
      RadiusMin
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption

    integer ( KDI ) :: &
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ) :: &
      Pi
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

    Pi  =  CONSTANT % PI

    CoordinateSystem = 'SPHERICAL'

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    IsPeriodic = .false.
    IsPeriodic ( 3 ) = .true.

    MinCoordinate = [ RadiusMin,     0.0_KDR,      0.0_KDR ]
    MaxCoordinate = [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]

    C % nCellsPolar = 128
    if ( present ( nCellsPolarOption ) ) &
      C % nCellsPolar = nCellsPolarOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsPolar, 'nCellsPolar' )

    call C % SetCore ( )

    C % RadialRatio = 1
    if ( present ( RadialRatioOption ) ) &
      C % RadialRatio = RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    nCellsRadial     =  C % RadialRatio * C % nCellsPolar !-- Aim for RadiusMax
    nCellsPolar      =  C % nCellsPolar
    nCellsAzimuthal  =  2 * nCellsPolar
 
    C % MinWidth  =  C % RadiusScale  *  Pi / nCellsPolar

    nCells = [ nCellsRadial, 1, 1 ]
    if ( Atlas % nDimensions > 1 ) &
      nCells ( 2 ) = nCellsPolar
    if ( Atlas % nDimensions > 2 ) &
      nCells ( 3 ) = nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  Pi / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusScale

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
             nEqualOption = nEqualOption )

    if ( C % nDimensions > 1 .and. C % nCells ( 2 ) /= C % nCellsPolar ) then
      call Show ( 'Choose nBricks such that nCells ( 2 ) need not be', &
                  CONSOLE % ERROR )
      call Show ( 'changed from requested nCellsPolar', &
                  CONSOLE % ERROR )
      call Show ( C % nCellsPolar, 'nCellsPolar', CONSOLE % ERROR )
      call Show ( C % nBricks ( 2 ), 'nBricks ( 2 )', CONSOLE % ERROR )
      call Show ( mod ( C % nCellsPolar, C % nBricks ( 2 ) ), &
                  'mod ( nCellsPolar, nBricks ( 2 ) )', CONSOLE % ERROR )
      call Show ( C % nCells ( 2 ), 'nCells ( 2 )', CONSOLE % ERROR )
      call Show ( 'InitializeTemplate_C', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Chart_SLD_C__Template', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if
      
  end subroutine InitializeTemplate_C


  subroutine ShowHeaderTemplate_C ( C )

    class ( Chart_SLD_C_Template ), intent ( in ) :: &
      C

    call Show ( 'Chart_SLD_C parameters' )
    call Show ( C % nCellsPolar, 'nCellsPolar', C % IGNORABILITY )
    call Show ( C % RadiusScale, C % CoordinateUnit ( 1 ), 'RadiusScale', &
                C % IGNORABILITY )
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth', &
                C % IGNORABILITY )
    call Show ( C % RadialRatio, 'RadialRatio' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), &
                'RadiusMax requested', C % IGNORABILITY )
    call Show ( C % MaxCoordinate ( 1 ), C % CoordinateUnit ( 1 ), &
                'RadiusMax actual', C % IGNORABILITY )

  end subroutine ShowHeaderTemplate_C


  subroutine SetCoarseningPolar ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iV, &          !-- iValue
      iB, jB, kB, &  !-- iBrick, etc.
      iC, jC, kC, &  !-- iCell, etc.
      iR, &          !-- iRank
      iRB, iRP, &    !-- iRankBrick, iRankPillar
      iP, &         !-- iPillar
      nRanks_2, &
      nPillarsTotal, &
      nPillarsSurplus
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank_2, &
      nPillars_2
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      nPillarsTransverse_2
    real ( KDR ) :: &
      Pi
    class ( GeometryFlatForm ), pointer :: &
      G

    Pi  =  CONSTANT % PI

    G => C % Geometry ( )

    associate &
      (   nB => C % nBricks, &
         nCB => C % nCellsBrick )
             
    if ( C % nDimensions < 2 ) &
      return

    call Show ( 'SetCoarseningPolar', C % IGNORABILITY )
    call Show ( C % Name, 'Chart', C % IGNORABILITY )

    associate &
      (    R_G => C % Center ( 1 ) % Value, &  !-- R_Global
        dTheta => Pi / C % nCellsPolar, &
             R => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Crsn_2 => G % Value ( :, G % COARSENING ( 2 ) ) )

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
        iLoop_2: do iB = 1, nB ( 1 )
          if ( nPillarsTransverse_2 ( iB, kB ) == 0 ) &
            cycle iLoop_2
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
        end do iLoop_2 !-- iB
      end do !-- kB
      call Show ( 'Communication counting', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsFrom_2, 'nSegmentsFrom_2', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsTo_2, 'nSegmentsTo_2', C % IGNORABILITY + 1 )

    end if !-- Comm_2 % Initialized
    end associate !-- Comm_2

    !-- Cleanup

    end associate !-- R_G, etc.
    end associate !-- nB, etc.
    nullify ( G )

  end subroutine SetCoarseningPolar


  subroutine SetCoarseningAzimuthal ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

    integer ( KDI ) :: &
      iV, &          !-- iValue
      iB, jB, kB, &  !-- iBrick, etc.
      iC, jC, kC, &  !-- iCell, etc.
      iR, &          !-- iRank
      iRB, iRP, &    !-- iRankBrick, iRankPillar
      iP, &         !-- iPillar
      nRanks_3, &
      nPillarsTotal, &
      nPillarsSurplus
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank_3, &
      nPillars_3
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      nPillarsTransverse_3
    real ( KDR ) :: &
      Pi
    class ( GeometryFlatForm ), pointer :: &
      G

    Pi  =  CONSTANT % PI

    G => C % Geometry ( )

    associate &
      (   nB => C % nBricks, &
         nCB => C % nCellsBrick )
             
    if ( C % nDimensions < 3 ) &
      return

    call Show ( 'SetCoarseningAzimuthal', C % IGNORABILITY )
    call Show ( C % Name, 'Chart', C % IGNORABILITY )

    associate &
      (     R_G => C % Center ( 1 ) % Value, &  !-- R_Global
        Theta_G => C % Center ( 2 ) % Value, &  !-- Theta_Global
           dPhi => 2 * Pi / ( 2 * C % nCellsPolar ), &
              R => G % Value ( :, G % CENTER_U ( 1 ) ), &
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
        iLoop_2: do iB = 1, nB ( 1 )
          if ( nPillarsTransverse_3 ( iB, jB ) == 0 ) &
            cycle iLoop_2
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
        end do iLoop_2 !-- iB
      end do !-- jB
      call Show ( 'Communication counting', C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsFrom_3, 'C % nSegmentsFrom_3', &
                  C % IGNORABILITY + 1 )
      call Show ( C % nSegmentsTo_3, 'C % nSegmentsTo_3', &
                  C % IGNORABILITY + 1 )

    end if !-- Comm_3 % Initialized
    end associate !-- Comm_3

    !-- Cleanup

    end associate !-- R_G, etc.
    end associate !-- nB, etc.
    nullify ( G )

  end subroutine SetCoarseningAzimuthal


  impure elemental subroutine FinalizeTemplate_C ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
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

  end subroutine FinalizeTemplate_C


  subroutine SetCore ( C )

    class ( Chart_SLD_C_Template ), intent ( inout ) :: &
      C

  end subroutine SetCore


end module Chart_SLD_C__Template
