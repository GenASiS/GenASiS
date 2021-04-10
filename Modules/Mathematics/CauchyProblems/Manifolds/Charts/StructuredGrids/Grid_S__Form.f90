module Grid_S__Form

  !-- Grid_Structured__Form

  use Basics
  use Chart_H__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3

  type, public, extends ( Chart_H_Form ) :: Grid_S_Form
    integer ( KDI ) :: &
      nEqual, &
      nCellsProper, &
      nCellsGhost, &
      nCellsLocal
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      iaFirst, &
      iaLast, &
      nCells, &
      nGhostLayers
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      iaBrick, &
      nBricks, &
      nCellsBrick
    real ( KDR ), dimension ( MAX_DIMENSIONS ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    logical ( KDL ) :: &
      Distributed
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      Edge, &
      Width, &
      Center
    logical ( KDL ), dimension ( : ), allocatable :: &
      ProperCell
    character ( LDL ), dimension ( MAX_DIMENSIONS ) :: &
      Spacing
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( PortalHeaderForm ), allocatable :: &
      PortalFace_L_R, &
      PortalFace_R_L, &
      PortalEdge_LL_RR, &
      PortalEdge_RR_LL, &
      PortalEdge_LR_RL, &
      PortalEdge_RL_LR
  contains
    procedure, private, pass :: &
      Initialize_GS
    generic, public :: &
      Initialize => Initialize_GS
    procedure, public, pass :: &
      ComputeCoordinateData
    procedure, private, pass :: &
      Show_C
    procedure, public, pass :: &
      SetFieldPointer_1D_3D
    generic, public :: &
      SetFieldPointer => SetFieldPointer_1D_3D
    final :: &
      Finalize
  end type Grid_S_Form

    private :: &
      SetCoordinateMetadata, &
      SetCells, &
      SetDecomposition

      private :: &
        BrickIndex, &
        SetPortals, &
        SetCellsLocal, &
        SetProperCells, &
        ComputeEdgeEqual, &
        ComputeGeometricRatio, &
        ComputeEdgeGeometric, &
        ComputeEdgeCompactified, &
        ComputeEdgeProportional

        private :: &
          ZeroGeometricRatio


contains


  subroutine Initialize_GS &
               ( G, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, PeriodicOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G
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

    integer ( KDI ) :: &
      iD  !-- iDimension
    logical ( KDL ), dimension ( MAX_DIMENSIONS ) :: &
      Periodic

    if ( G % Type  ==  '' ) &
      G % Type  =  'a Grid_S'

    Periodic  =  .false.
    if ( present ( PeriodicOption ) ) &
      Periodic ( : size ( PeriodicOption ) )  =  PeriodicOption

    call G % Chart_H_Form % Initialize &
           ( Periodic, CoordinateLabelOption, CoordinateSystemOption, &
             NameOption, CoordinateUnitOption, IgnorabilityOption, &
             nDimensionsOption, iDimensionalityOption )

    call SetCoordinateMetadata &
           ( G, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
             RatioOption, ScaleOption, nEqualOption )

    call SetCells &
           ( G, nCellsOption, nGhostLayersOption )

    call SetDecomposition &
           ( G, CommunicatorOption, nBricksOption, nBricksCompatibleOption )

    do iD = 1, G % nDimensions
      call ComputeCoordinateData ( G, iD )
    end do !-- iD

  end subroutine Initialize_GS


  subroutine ComputeCoordinateData ( G, iD, EdgeValueOption )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G
    integer ( KDI ), intent ( in ) :: &
      iD      !-- iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      EdgeValueOption

    integer ( KDI ) :: &
      iC    !-- iCell
    real ( KDL ) :: &
      Width_IG, &
      Width_OG

    if ( .not. allocated ( G % Edge ) ) &
      allocate ( G % Edge ( MAX_DIMENSIONS ) )
    if ( .not. allocated ( G % Width ) ) &
      allocate ( G % Width  ( MAX_DIMENSIONS ) )
    if ( .not. allocated ( G % Center ) ) &
      allocate ( G % Center ( MAX_DIMENSIONS ) )

    associate &
      (  nC => G % nCells ( iD ), &
        nGL => G % nGhostLayers ( iD ) )

    if ( .not. allocated ( G % Edge ( iD ) % Value ) ) &
      call G % Edge ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL + 1, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( G % Width ( iD ) % Value ) ) &
      call G % Width ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( G % Center ( iD ) % Value ) ) &
      call G % Center ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )

    !-- Edge, proper cells
    if ( present ( EdgeValueOption ) ) then
      G % Edge ( iD ) % Value ( 1 : nC + 1 )  =  EdgeValueOption
      G % MinCoordinate ( iD )  =  EdgeValueOption ( 1 )
      G % MaxCoordinate ( iD )  =  EdgeValueOption ( nC + 1 )
    else
      select case ( trim ( G % Spacing ( iD ) ) )
      case ( 'EQUAL' )
        call ComputeEdgeEqual &
               ( G % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 G % MinCoordinate ( iD ), G % MaxCoordinate ( iD ), nC )
      case ( 'GEOMETRIC' )
        if ( G % Scale ( iD ) > 0.0_KDR ) &
          call ComputeGeometricRatio &
                 ( G % CoordinateUnit ( iD ), G % MinCoordinate ( iD ), &
                   G % MaxCoordinate ( iD ), G % Scale ( iD ), nC, &
                   G % Ratio ( iD ) )
        call ComputeEdgeGeometric &
               ( G % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 G % MinCoordinate ( iD ), G % MaxCoordinate ( iD ), &
                 G % Ratio ( iD ), nC )
      case ( 'COMPACTIFIED' )
        call ComputeEdgeCompactified &
               ( G % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 G % Scale ( iD ), nC )
        G % MinCoordinate ( iD )  =  G % Edge ( iD ) % Value ( 1 )
        G % MaxCoordinate ( iD )  =  G % Edge ( iD ) % Value ( nC + 1 )
      case ( 'PROPORTIONAL' )
        call ComputeEdgeProportional &
               ( G % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 G % MinCoordinate ( iD ), G % Ratio ( iD ), &
                 G % Scale ( iD ), nC, G % nEqual )
        G % MaxCoordinate ( iD )  =  G % Edge ( iD ) % Value ( nC + 1 )
      case default
        call Show ( 'Spacing not recognized', CONSOLE % ERROR )
        call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeCoordinateData', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end if

    !-- Edge, ghost cells
    associate ( Edge => G % Edge ( iD ) % Value )
    do iC = 1, nGL
      Width_IG  =  Edge ( iC + 1 )       -  Edge ( iC )
      Width_OG  =  Edge ( nC - iC + 2 )  -  Edge ( nC - iC + 1 )
      Edge ( 1 - iC )       =  Edge ( 2 - iC )   -  Width_IG
      Edge ( nC + 1 + iC )  =  Edge ( nC + iC )  +  Width_OG
    end do !-- iC
    end associate !-- Edge

    !-- Width
    associate &
      ( Edge  => G % Edge ( iD ) % Value, &
        Width => G % Width ( iD ) % Value )
    do iC = lbound ( Width, dim = 1 ), ubound ( Width, dim = 1 )
      Width ( iC )  =  Edge ( iC + 1 )  -  Edge ( iC )
    end do !-- iC
    end associate !-- Edge, etc.

    !-- Center
    associate &
      (   Edge => G % Edge ( iD ) % Value, &
        Center => G % Center ( iD ) % Value )
    do iC = lbound ( Center, dim = 1 ), ubound ( Center, dim = 1 )
      Center ( iC )  =  0.5_KDR * ( Edge ( iC )  +  Edge ( iC + 1 ) )
    end do !-- iC
    end associate !-- Edge, etc.

    end associate !-- nC, etc.

  end subroutine ComputeCoordinateData


  subroutine Show_C ( C )

    class ( Grid_S_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iD  !-- iDimension

    call C % Chart_H_Form % Show ( )

    associate ( nD => C % nDimensions )

    call Show ( C % MinCoordinate ( : nD ), C % CoordinateUnit ( : nD ), &
                'MinCoordinate', C % IGNORABILITY )
    call Show ( C % MaxCoordinate ( : nD ), C % CoordinateUnit ( : nD ), &
                'MaxCoordinate', C % IGNORABILITY )

    call Show ( C % Spacing ( : nD ), 'Spacing', C % IGNORABILITY )
    if ( any ( C % Spacing == 'GEOMETRIC' ) &
         .or. any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % Ratio ( : nD ), 'Ratio', C % IGNORABILITY )
    if ( any ( C % Spacing == 'GEOMETRIC' ) &
         .or. any ( C % Spacing == 'COMPACTIFIED' ) &
         .or. any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % Scale ( : nD ), C % CoordinateUnit, 'Scale', &
                  C % IGNORABILITY )
    if ( any ( C % Spacing == 'PROPORTIONAL' ) ) &
      call Show ( C % nEqual, 'nEqual', C % IGNORABILITY )

    call Show ( C % nCells ( : nD ),       'nCells',       C % IGNORABILITY )
    call Show ( C % nGhostLayers ( : nD ), 'nGhostLayers', C % IGNORABILITY )

    call Show ( C % Distributed, 'Distributed', C % IGNORABILITY )
    if ( C % Distributed ) then
      call Show ( C % Communicator % Name,  'Communicator', &
                  C % IGNORABILITY )
      call Show ( C % nBricks ( : nD ),     'nBricks', &
                  C % IGNORABILITY )
      call Show ( C % iaBrick ( : nD ),     'iaBrick',      &
                  C % IGNORABILITY + 1 )
      call Show ( C % nCellsBrick ( : nD ), 'nCellsBrick', &
                  C % IGNORABILITY )
    end if !-- Distributed

    call Show ( C % iaFirst ( : nD ), 'iaFirst', C % IGNORABILITY + 1 )
    call Show ( C % iaLast  ( : nD ), 'iaLast',  C % IGNORABILITY + 1 )

    call Show ( C % nCellsProper, 'nCellsProper', C % IGNORABILITY )
    call Show ( C % nCellsGhost,  'nCellsGhost',  C % IGNORABILITY )
    call Show ( C % nCellsLocal,  'nCellsLocal',  C % IGNORABILITY )

    if ( C % Distributed ) then
      call C % PortalFace_L_R % Show &
             ( 'PortalFace_L_R', C % IGNORABILITY + 1 )
      call C % PortalFace_R_L % Show &
             ( 'PortalFace_R_L', C % IGNORABILITY + 1 )
      call C % PortalEdge_LL_RR % Show &
             ( 'PortalEdge_LL_RR', C % IGNORABILITY + 1 )
      call C % PortalEdge_RR_LL % Show &
             ( 'PortalEdge_RR_LL', C % IGNORABILITY + 1 )
      call C % PortalEdge_LR_RL % Show &
             ( 'PortalEdge_LR_RL', C % IGNORABILITY + 1 )
      call C % PortalEdge_RL_LR % Show &
             ( 'PortalEdge_RL_LR', C % IGNORABILITY + 1 )
    end if !-- Distributed

    do iD = 1, nD
      call Show ( iD, 'iDimension', C % IGNORABILITY + 1 )
      call Show ( C % Edge ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Edge', C % IGNORABILITY + 1 )
      call Show ( C % Width ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Width', C % IGNORABILITY + 1 )
      call Show ( C % Center ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Center', C % IGNORABILITY + 1 )
    end do !-- iD

    end associate !-- nD

  end subroutine Show_C


  subroutine SetFieldPointer_1D_3D ( G, Field_1D, Field_3D )

    class ( Grid_S_Form ), intent ( in ) :: &
      G
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Field_1D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      Field_3D

    Field_3D &
      ( G % iaFirst ( 1 ) : G % iaLast ( 1 ), &
        G % iaFirst ( 2 ) : G % iaLast ( 2 ), &
        G % iaFirst ( 3 ) : G % iaLast ( 3 ) ) &
          => Field_1D
    
  end subroutine SetFieldPointer_1D_3D

  
  impure elemental subroutine Finalize ( G )

    type ( Grid_S_Form ), intent ( inout ) :: &
      G

    nullify ( G % Communicator )

    if ( allocated ( G % PortalEdge_RL_LR ) ) &
      deallocate ( G % PortalEdge_RL_LR )
    if ( allocated ( G % PortalEdge_LR_RL ) ) &
      deallocate ( G % PortalEdge_LR_RL )
    if ( allocated ( G % PortalEdge_RR_LL ) ) &
      deallocate ( G % PortalEdge_RR_LL )
    if ( allocated ( G % PortalEdge_LL_RR ) ) &
      deallocate ( G % PortalEdge_LL_RR )
    if ( allocated ( G % PortalFace_R_L ) ) &
      deallocate ( G % PortalFace_R_L )
    if ( allocated ( G % PortalFace_L_R ) ) &
      deallocate ( G % PortalFace_L_R )

    if ( allocated ( G % ProperCell ) ) &
      deallocate ( G % ProperCell )

    if ( allocated ( G % Center ) ) &
      deallocate ( G % Center )
    if ( allocated ( G % Width ) ) &
      deallocate ( G % Width )
    if ( allocated ( G % Edge ) ) &
      deallocate ( G % Edge )

  end subroutine Finalize


  subroutine SetCoordinateMetadata &
               ( G, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
                 RatioOption, ScaleOption, nEqualOption )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), intent ( in ), optional :: &
      nEqualOption

    associate ( nD => G % nDimensions )

    G % MinCoordinate = 0.0_KDR
    if ( present ( MinCoordinateOption ) ) &
      G % MinCoordinate ( : nD ) = MinCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( G % MinCoordinate ( : nD ), 'MinCoordinate', &
             InputUnitOption = G % CoordinateUnit ( : nD ) )

    G % MaxCoordinate = 0.0_KDR
    G % MaxCoordinate ( : nD ) = 1.0_KDR
    if ( present ( MaxCoordinateOption ) ) &
      G % MaxCoordinate ( : nD ) = MaxCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( G % MaxCoordinate ( : nD ), 'MaxCoordinate', &
             InputUnitOption = G % CoordinateUnit ( : nD ) )

    G % Spacing = ''
    G % Spacing ( : nD ) = 'EQUAL'
    if ( present ( SpacingOption ) ) &
      G % Spacing ( : nD ) = SpacingOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( G % Spacing ( : nD ), 'Spacing' )

    G % Ratio = 0.0_KDR
    if ( present ( RatioOption ) ) &
      G % Ratio ( : nD ) = RatioOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( G % Ratio ( : nD ), 'Ratio' )

    G % Scale = 0.0_KDR
    if ( present ( ScaleOption ) ) &
      G % Scale ( : nD ) = ScaleOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( G % Scale ( : nD ), 'Scale' )

    G % nEqual = 0
    if ( present ( nEqualOption ) ) &
      G % nEqual = nEqualOption

    end associate !-- nD

  end subroutine SetCoordinateMetadata


  subroutine SetCells ( G, nCellsOption, nGhostLayersOption )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption

    associate ( nD => G % nDimensions )

    G % nCells = 1
    G % nCells ( : nD ) = 32
    if ( present ( nCellsOption ) ) &
      G % nCells ( : nD ) = nCellsOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( G % nCells ( : nD ), 'nCells' )

    G % nGhostLayers = 0
    G % nGhostLayers ( : nD ) = 2
    if ( present ( nGhostLayersOption ) ) &
      G % nGhostLayers ( : nD ) = nGhostLayersOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( G % nGhostLayers ( : nD ), 'nGhostLayers' )

    end associate !-- nD

  end subroutine SetCells


  subroutine SetDecomposition &
                ( G, CommunicatorOption, nBricksOption, &
                  nBricksCompatibleOption )
               
    class ( Grid_S_Form ), intent ( inout ) :: &
      G
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nBricksOption, &
      nBricksCompatibleOption

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      SizeRoot
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      nBricksCompatible

    if ( present ( CommunicatorOption ) ) then
      G % Distributed   =   .true.
      G % Communicator  =>  CommunicatorOption
    else
      G % Distributed  =  .false.
    end if !-- present Communicator 

    if ( G % Distributed ) then

      associate ( nD => G % nDimensions )

      SizeRoot  =  G % Communicator % Size ** ( 1.0_KDR / nD ) + 0.5_KDR

      G % nBricks = 1
      G % nBricks ( : nD ) = SizeRoot
      if ( present ( nBricksOption ) ) &
        G % nBricks  =  nBricksOption 
      call PROGRAM_HEADER % GetParameter ( G % nBricks ( : nD ), 'nBricks' )
    
      nBricksCompatible = G % nBricks
      if ( present ( nBricksCompatibleOption ) ) &
        nBricksCompatible  =  nBricksCompatibleOption 
      call PROGRAM_HEADER % GetParameter &
             ( nBricksCompatible ( : nD ), 'nBricksCompatible' )

      if ( any ( nBricksCompatible /= G % nBricks ) ) then
        call Show ( 'nBricksCompatible /= nBricks', CONSOLE % INFO_1 )
        call Show ( G % nBricks, 'nBricks', CONSOLE % INFO_1 )
        call Show ( nBricksCompatible, 'nBricksCompatible', CONSOLE % INFO_1 )
      end if

      if ( product ( G % nBricks ) /= G % Communicator % Size ) then
        call Show ( 'The total number of bricks must equal ' &
                    // 'the number of MPI processes', CONSOLE % ERROR )
        call Show ( G % Communicator % Size, 'nProcesses', CONSOLE % ERROR )
        call Show ( G % nBricks ( 1 : nD ), 'nBricks', CONSOLE % ERROR )
        call Show ( product ( G % nBricks ), 'product ( nBricks )', &
                    CONSOLE % ERROR )
        call Show ( 'Grid_S__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetDecomposition', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      do iD = 1, nD
        if ( mod ( G % nCells ( iD ), nBricksCompatible ( iD ) ) /= 0 ) then
          call Show ( 'nBricksCompatible in each dimension must divide ' &
                      // 'evenly into nCells in each dimension', &
                      CONSOLE % WARNING )
          call Show ( iD, 'iDimension', CONSOLE % WARNING )
          call Show ( nBricksCompatible ( iD ), 'nBricksCompatible', &
                      CONSOLE % WARNING )
          call Show ( G % nCells ( iD ), 'nCells requested', CONSOLE % WARNING )
          G % nCells ( iD ) &
            =  ( G % nCells ( iD ) / nBricksCompatible ( iD ) ) &
               *  nBricksCompatible ( iD )
          call Show ( G % nCells ( iD ), 'nCells granted', CONSOLE % WARNING )
          call Show ( 'SetDecomposition', 'subroutine', CONSOLE % WARNING )
          call Show ( 'Grid_S__Form', 'module', CONSOLE % WARNING )
        end if
      end do  !-- iD
    
      G % nCellsBrick &
        = G % nCells / G % nBricks
      G % iaBrick &
        = BrickIndex ( G % nBricks, G % nCells, G % Communicator % Rank )

      end associate !-- nD

      call SetPortals ( G )

      call SetCellsLocal ( G, G % nCellsBrick )

    else  !-- not Distributed

      call SetCellsLocal ( G, G % nCells )

    end if  !-- Distributed

    call SetProperCells ( G )

  end subroutine SetDecomposition


  function BrickIndex ( nBricks, nCells, MyRank )  result ( BI ) 

    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nBricks, &
      nCells
    integer ( KDI ), intent ( in ) :: &
      MyRank
    integer ( KDI ) , dimension ( MAX_DIMENSIONS )  :: &
      BI

    associate ( nB => nBricks )

    if ( nCells ( 3 ) > 1 ) then
      BI ( 3 ) = ( MyRank / ( nB ( 1 ) * nB ( 2 ) ) ) + 1
    else
      BI ( 3 ) = 1
    end if

    if ( nCells ( 2 ) > 1 ) then
      BI ( 2 ) = ( mod ( MyRank, nB ( 1 ) * nB ( 2 ) ) / nB ( 1 ) ) + 1
    else
      BI ( 2 ) = 1
    end if

    BI ( 1 ) = mod ( mod ( MyRank, nB ( 1 ) * nB ( 2 ) ), nB ( 1 ) ) + 1

    end associate !-- nB

  end function BrickIndex


  subroutine SetPortals ( G )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iB, jB, kB, &  !-- iBrick, etc.
      iD, jD, kD     !-- iDimension, etc.
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      iaB
    integer ( KDI ), dimension ( : ), allocatable :: &
      !-- Faces
      nCellsFace, &
      Source_L_R, &
      Source_R_L, &
      Target_L_R, &
      Target_R_L, &
      !-- Edges
      nCellsEdge, &
      Source_LL_RR, &
      Source_RR_LL, &
      Source_LR_RL, &
      Source_RL_LR, &
      Target_LL_RR, &
      Target_RR_LL, &
      Target_LR_RL, &
      Target_RL_LR
    integer ( KDI ), dimension ( :, :, : ), allocatable :: &
      Process

    associate &
      ( nB   =>  G % nBricks, &
        nD   =>  G % nDimensions, &
        nGL  =>  G % nGhostLayers, &
        nCB  =>  G % nCellsBrick ) 

    allocate ( Process ( nB ( 1 ), nB ( 2 ), nB ( 3 ) ) )

    iP = 0
    do kB = 1, nB ( 3 )
      do jB = 1, nB ( 2 )
        do iB = 1, nB ( 1 )
          Process ( iB, jB, kB ) = iP
          iP = iP + 1
        end do !-- iB
      end do !-- jB
    end do !-- kB


    !-- Face sibling bricks

    allocate ( nCellsFace ( nD ) )
    allocate ( Source_L_R ( nD ) )
    allocate ( Source_R_L ( nD ) )
    allocate ( Target_L_R ( nD ) )
    allocate ( Target_R_L ( nD ) )
    nCellsFace =  0
    Source_L_R = -1
    Source_R_L = -1
    Target_L_R = -1
    Target_R_L = -1
    
    do iD = 1, nD

      jD  =  mod ( iD, 3 ) + 1
      kD  =  mod ( jD, 3 ) + 1

      nCellsFace ( iD )  =  nGL ( iD )  *  nCB ( jD )  *  nCB ( kD )

      iaB  =  G % iaBrick

      !-- Left brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      Source_L_R ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_R_L ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- Right brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ), nB ( iD ) ) + 1
      Source_R_L ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_L_R ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- iD


    !-- Edge sibling bricks

    allocate ( nCellsEdge ( 3 ) )
    allocate ( Source_LL_RR ( 3 ) )
    allocate ( Source_RR_LL ( 3 ) )
    allocate ( Source_LR_RL ( 3 ) )
    allocate ( Source_RL_LR ( 3 ) )
    allocate ( Target_LL_RR ( 3 ) )
    allocate ( Target_RR_LL ( 3 ) )
    allocate ( Target_LR_RL ( 3 ) )
    allocate ( Target_RL_LR ( 3 ) )
    nCellsEdge =  0
    Source_LL_RR = -1
    Source_RR_LL = -1
    Source_LR_RL = -1
    Source_RL_LR = -1
    Target_LL_RR = -1
    Target_RR_LL = -1
    Target_LR_RL = -1
    Target_RL_LR = -1

    do kD = 3, 1, -1

      iD  =  mod ( kD, 3 ) + 1
      jD  =  mod ( iD, 3 ) + 1

      if ( iD > nD .or. jD > nD ) &
        cycle

      nCellsEdge ( kD )  &
        =  G % nGhostLayers ( iD )  *  G % nGhostLayers ( jD )  &
           *  G % nCellsBrick ( kD )

      iaB  =  G % iaBrick

      !-- LeftLeft brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( G % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_LL_RR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightRight brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( G % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LL_RR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- LeftRight brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( G % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RL_LR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightLeft brick

      iaB ( iD )  &
        =  mod ( G % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( G % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_RL_LR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- kD


    !-- Set face portals
    
    allocate ( G % PortalFace_L_R )
    associate ( PFLR => G % PortalFace_L_R )
    call PFLR % Initialize ( Source_L_R, Target_L_R, nCellsFace, nCellsFace )
    end associate !-- PFLR

    allocate ( G % PortalFace_R_L )
    associate ( PFRL => G % PortalFace_R_L )
    call PFRL % Initialize ( Source_R_L, Target_R_L, nCellsFace, nCellsFace )
    end associate !-- PFRL


    !-- Set edge portals
    
    allocate ( G % PortalEdge_LL_RR )
    associate ( PELLRR => G % PortalEdge_LL_RR )
    call PELLRR % Initialize &
           ( pack ( Source_LL_RR, Source_LL_RR >= 0 ), &
             pack ( Target_LL_RR, Target_LL_RR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELLRR

    allocate ( G % PortalEdge_RR_LL )
    associate ( PERRLL => G % PortalEdge_RR_LL )
    call PERRLL % Initialize &
           ( pack ( Source_RR_LL, Source_RR_LL >= 0 ), &
             pack ( Target_RR_LL, Target_RR_LL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PERRLL

    allocate ( G % PortalEdge_LR_RL )
    associate ( PELRRL => G % PortalEdge_LR_RL )
    call PELRRL % Initialize &
           ( pack ( Source_LR_RL, Source_LR_RL >= 0 ), &
             pack ( Target_LR_RL, Target_LR_RL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELRRL

    allocate ( G % PortalEdge_RL_LR )
    associate ( PERLLR => G % PortalEdge_RL_LR )
    call PERLLR % Initialize &
           ( pack ( Source_RL_LR, Source_RL_LR >= 0 ), &
             pack ( Target_RL_LR, Target_RL_LR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELRRL


    !-- Cleanup

    end associate !-- nB, etc.
    
  end subroutine SetPortals


  subroutine SetCellsLocal ( G, nCellsLocal )

    class ( Grid_S_Form ), intent ( inout ) :: &
      G
    integer, dimension ( : ), intent ( in ) :: &
      nCellsLocal

    G % iaFirst = 1
    G % iaLast = 1
    G % iaFirst ( 1 ) = 1 - G % nGhostLayers ( 1 )
    G % iaLast  ( 1 ) = nCellsLocal ( 1 ) + G % nGhostLayers ( 1 )
    if ( G % nDimensions > 1 ) then
      G % iaFirst ( 2 ) = 1 - G % nGhostLayers ( 2 )
      G % iaLast  ( 2 ) = nCellsLocal ( 2 ) + G % nGhostLayers ( 2 )
    end if
    if ( G % nDimensions > 2 ) then
      G % iaFirst ( 3 ) = 1 - G % nGhostLayers ( 3 )
      G % iaLast  ( 3 ) = nCellsLocal ( 3 ) + G % nGhostLayers ( 3 )
    end if

    G % nCellsLocal   =  product ( nCellsLocal  +  2 * G % nGhostLayers )
    G % nCellsProper  =  product ( nCellsLocal )
    G % nCellsGhost   =  G % nCellsLocal  -  G % nCellsProper

  end subroutine SetCellsLocal


  subroutine SetProperCells ( G )

    class ( Grid_S_Form ), intent ( inout ), target :: &
      G

    integer ( KDI ) :: &
      iC, jC, kC, &
      iV
    logical ( KDL ), dimension ( :, :, : ), pointer :: &
      PC

    allocate ( G % ProperCell ( G % nCellsLocal ) )
    call Clear ( G % ProperCell )

    associate &
      ( iaF  =>  G % iaFirst, &
        iaL  =>  G % iaLast, &
        nGL  =>  G % nGhostLayers )

    PC ( iaF ( 1 ) : iaL ( 1 ), &
         iaF ( 2 ) : iaL ( 2 ), &
         iaF ( 3 ) : iaL ( 3 ) )  &
      =>  G % ProperCell

    associate &
      ( lB  =>  iaF + nGL, &
        uB  =>  iaL - nGL )
    !$OMP parallel do private collapse ( 3 )
    do kC = lB ( 3 ), uB ( 3 )
      do jC = lB ( 2 ), uB ( 2 )
        do iC = lB ( 1 ), uB ( 1 )
          PC ( iC, jC, kC )  =  .true.
        end do
      end do
    end do
    !$OMP end parallel do
    end associate !-- lB, etc.

    end associate !-- iaF, etc.

    nullify ( PC )

  end subroutine SetProperCells


  subroutine ComputeEdgeEqual ( Edge, MinCoordinate, MaxCoordinate, nC )

    !-- Equal cell widths

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) / nC

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
    end do

  end subroutine ComputeEdgeEqual


  subroutine ComputeGeometricRatio &
               ( CoordinateUnit, MinCoordinate, MaxCoordinate, MinWidth, &
                 nCells, Ratio )

    type ( MeasuredValueForm ), intent ( in ) :: &
      CoordinateUnit
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      MinWidth
    integer ( KDI ), intent ( in ) :: &
      nCells
    real ( KDR ), intent ( out ) :: &
      Ratio

    integer ( KDI ) :: &
      i
    real ( KDR ) :: &
      a, b, c, &
      fa, fb, fc

    a  = 1.000001_KDR
    b  = 2.0_KDR
    fa = ZeroGeometricRatio &
           ( a, MinCoordinate, MaxCoordinate, MinWidth, nCells )
    fb = ZeroGeometricRatio &
           ( b, MinCoordinate, MaxCoordinate, MinWidth, nCells )
    if ( fa * fb > 0.0_KDR ) then
      call Show ( 'Solution not bracketed', CONSOLE % ERROR )
      call Show ( a, 'a', CONSOLE % ERROR )
      call Show ( b, 'b', CONSOLE % ERROR )
      call Show ( fa, 'f(a)', CONSOLE % ERROR )
      call Show ( fb, 'f(b)', CONSOLE % ERROR )
      call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'ComputeGeometricRatio', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    !-- Bisection
    do i = 1, 100
      c  = 0.5_KDR * ( a + b )
      if ( ( b - a ) / c  <  1.0e-10_KDR ) then
        Ratio = c
        return
      end if
      fc = ZeroGeometricRatio &
           ( c, MinCoordinate, MaxCoordinate, MinWidth, nCells )
      if ( sign ( 1.0_KDR, fc )  ==  sign ( 1.0_KDR, fa ) ) then
        a  = c
        fa = fc
      else
        b  = c
        fb = fc
      end if
    end do !-- i

    call Show ( 'GeometricRatio failed to converge', CONSOLE % ERROR )
    call Show ( MinCoordinate, CoordinateUnit, 'MinCoordinate', &
                CONSOLE % ERROR )
    call Show ( MaxCoordinate, CoordinateUnit, 'MaxCoordinate', &
                CONSOLE % ERROR )
    call Show ( MinWidth, CoordinateUnit, 'MinWidth', &
                CONSOLE % ERROR )
    call Show ( Ratio, 'Ratio', CONSOLE % ERROR )
    call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
    call Show ( 'ComputeGeometricRatio', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine ComputeGeometricRatio


  subroutine ComputeEdgeGeometric &
               ( Edge, MinCoordinate, MaxCoordinate, Ratio, nC )

    !-- Each successive cell width is larger by Ratio

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    Edge ( 1 )  =  MinCoordinate
    Width       =  ( MaxCoordinate - MinCoordinate ) &
                   * ( Ratio - 1.0_KDR ) / ( Ratio ** nC  -  1.0_KDR )

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio * Width
    end do

  end subroutine ComputeEdgeGeometric


  subroutine ComputeEdgeCompactified ( Edge, Scale, nC )

    !-- Compactify the domain [ 0, Infinity ] to [ 0, 1 ] via the
    !   transformation Coordinate = Scale * S / ( 1 - S )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nC

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      dS, &
      S, &
      Width

    dS = 1.0_KDR / nC
    S  = 0.5_KDR * dS

    Edge ( 1 )  =  0.0_KDR
    Width       =  Scale  *  dS / ( 1.0_KDR - S ) ** 2 

    do iC = 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      S            =  ( 0.5_KDR + ( iC - 1 ) ) * dS
      Width        =  Scale  *  dS / ( 1.0_KDR - S ) ** 2
    end do

  end subroutine ComputeEdgeCompactified


  subroutine ComputeEdgeProportional &
               ( Edge, MinCoordinate, Ratio, Scale, nC, nEqual )

    !-- Width proportional to the inner edge coordinate of the cell

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Edge
    real ( KDR ), intent ( in ) :: &
      MinCoordinate, &
      Ratio, &
      Scale
    integer ( KDI ), intent ( in ) :: &
      nC, &
      nEqual

    integer ( KDI ) :: &
      iC  !-- iCell
    real ( KDR ) :: &
      Width

    if ( nEqual == 0 ) then
      Edge ( 1 )  =  MinCoordinate
    else
      call ComputeEdgeEqual ( Edge, MinCoordinate, Scale, nEqual )
    end if
    
    Width  =  Ratio  *  Edge ( nEqual + 1 )

    do iC = nEqual + 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio  *  Edge ( iC )
    end do

  end subroutine ComputeEdgeProportional


  function ZeroGeometricRatio &
             ( Ratio, MinCoordinate, MaxCoordinate, MinWidth, nCells) &
             result ( ZGR )

    real ( KDR ), intent ( in ) :: &
      Ratio, &
      MinCoordinate, &
      MaxCoordinate, &
      MinWidth
    integer ( KDI ), intent ( in ) :: &
      nCells
    real ( KDR ) :: &
      ZGR

    ZGR  =  ( MaxCoordinate - MinCoordinate ) * ( Ratio - 1.0_KDR ) &
              /  ( Ratio ** nCells  -  1.0_KDR ) &
            -  MinWidth

  end function ZeroGeometricRatio


end module Grid_S__Form
