module Chart_GS__Form

  !-- Chart_GridStructured__Form

  use Basics
  use Chart_H__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3

  type, public, extends ( Chart_H_Form ) :: Chart_GS_Form
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
    procedure, public, pass :: &
      SetFieldPointer_2D_4D
    generic, public :: &
      SetFieldPointer => SetFieldPointer_1D_3D, SetFieldPointer_2D_4D
    final :: &
      Finalize
  end type Chart_GS_Form

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
               ( C, CommunicatorOption, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, NameOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, IgnorabilityOption, &
                 nDimensionsOption, nEqualOption, iDimensionalityOption )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
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

    integer ( KDI ) :: &
      iD  !-- iDimension

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS'

    call C % Initialize_H &
           ( CoordinateLabelOption, CoordinateSystemOption, NameOption, &
             CoordinateUnitOption, IgnorabilityOption, &
             nDimensionsOption, iDimensionalityOption )

    call SetCoordinateMetadata &
           ( C, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
             RatioOption, ScaleOption, nEqualOption )

    call SetCells &
           ( C, nCellsOption, nGhostLayersOption )

    call SetDecomposition &
           ( C, CommunicatorOption, nBricksOption, nBricksCompatibleOption )

    do iD = 1, C % nDimensions
      call ComputeCoordinateData ( C, iD )
    end do !-- iD

  end subroutine Initialize_GS


  subroutine ComputeCoordinateData ( C, iD, EdgeValueOption )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iD      !-- iDimension
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      EdgeValueOption

    integer ( KDI ) :: &
      iC    !-- iCell
    real ( KDL ) :: &
      Width_IG, &
      Width_OG

    if ( .not. allocated ( C % Edge ) ) &
      allocate ( C % Edge ( MAX_DIMENSIONS ) )
    if ( .not. allocated ( C % Width ) ) &
      allocate ( C % Width  ( MAX_DIMENSIONS ) )
    if ( .not. allocated ( C % Center ) ) &
      allocate ( C % Center ( MAX_DIMENSIONS ) )

    associate &
      (  nC => C % nCells ( iD ), &
        nGL => C % nGhostLayers ( iD ) )

    if ( .not. allocated ( C % Edge ( iD ) % Value ) ) &
      call C % Edge ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL + 1, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( C % Width ( iD ) % Value ) ) &
      call C % Width ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )
    if ( .not. allocated ( C % Center ( iD ) % Value ) ) &
      call C % Center ( iD ) % Initialize &
             ( nValues  =  nC  +  2 * nGL, &
               iLowerBoundOption  =  1 - nGL )

    !-- Edge, proper cells
    if ( present ( EdgeValueOption ) ) then
      C % Edge ( iD ) % Value ( 1 : nC + 1 )  =  EdgeValueOption
      C % MinCoordinate ( iD )  =  EdgeValueOption ( 1 )
      C % MaxCoordinate ( iD )  =  EdgeValueOption ( nC + 1 )
    else
      select case ( trim ( C % Spacing ( iD ) ) )
      case ( 'EQUAL' )
        call ComputeEdgeEqual &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), nC )
      case ( 'GEOMETRIC' )
        if ( C % Scale ( iD ) > 0.0_KDR ) &
          call ComputeGeometricRatio &
                 ( C % CoordinateUnit ( iD ), C % MinCoordinate ( iD ), &
                   C % MaxCoordinate ( iD ), C % Scale ( iD ), nC, &
                   C % Ratio ( iD ) )
        call ComputeEdgeGeometric &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), &
                 C % Ratio ( iD ), nC )
      case ( 'COMPACTIFIED' )
        call ComputeEdgeCompactified &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % Scale ( iD ), nC )
        C % MinCoordinate ( iD )  =  C % Edge ( iD ) % Value ( 1 )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case ( 'PROPORTIONAL' )
        call ComputeEdgeProportional &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % Ratio ( iD ), &
                 C % Scale ( iD ), nC, C % nEqual )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case default
        call Show ( 'Spacing not recognized', CONSOLE % ERROR )
        call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeCoordinateData', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end if

    !-- Edge, ghost cells
    associate ( Edge => C % Edge ( iD ) % Value )
    do iC = 1, nGL
      Width_IG  =  Edge ( iC + 1 )       -  Edge ( iC )
      Width_OG  =  Edge ( nC - iC + 2 )  -  Edge ( nC - iC + 1 )
      Edge ( 1 - iC )       =  Edge ( 2 - iC )   -  Width_IG
      Edge ( nC + 1 + iC )  =  Edge ( nC + iC )  +  Width_OG
    end do !-- iC
    end associate !-- Edge

    !-- Width
    associate &
      ( Edge  => C % Edge ( iD ) % Value, &
        Width => C % Width ( iD ) % Value )
    do iC = lbound ( Width, dim = 1 ), ubound ( Width, dim = 1 )
      Width ( iC )  =  Edge ( iC + 1 )  -  Edge ( iC )
    end do !-- iC
    end associate !-- Edge, etc.

    !-- Center
    associate &
      (   Edge => C % Edge ( iD ) % Value, &
        Center => C % Center ( iD ) % Value )
    do iC = lbound ( Center, dim = 1 ), ubound ( Center, dim = 1 )
      Center ( iC )  =  0.5_KDR * ( Edge ( iC )  +  Edge ( iC + 1 ) )
    end do !-- iC
    end associate !-- Edge, etc.

    end associate !-- nC, etc.

  end subroutine ComputeCoordinateData


  subroutine Show_C ( C )

    class ( Chart_GS_Form ), intent ( in ) :: &
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
    if ( C % Distributed ) &
      call Show ( C % Communicator % Name,  'Communicator', &
                  C % IGNORABILITY )
    call Show ( C % nBricks ( : nD ),     'nBricks', &
                C % IGNORABILITY )
    call Show ( C % iaBrick ( : nD ),     'iaBrick',      &
                C % IGNORABILITY + 1 )
    call Show ( C % nCellsBrick ( : nD ), 'nCellsBrick', &
                C % IGNORABILITY )

    call Show ( C % iaFirst ( : nD ), 'iaFirst', C % IGNORABILITY + 1 )
    call Show ( C % iaLast  ( : nD ), 'iaLast',  C % IGNORABILITY + 1 )

    call Show ( C % nCellsProper, 'nCellsProper', C % IGNORABILITY )
    call Show ( C % nCellsGhost,  'nCellsGhost',  C % IGNORABILITY )
    call Show ( C % nCellsLocal,  'nCellsLocal',  C % IGNORABILITY )

    if ( C % Distributed ) then
      call C % PortalFace_L_R % Show &
             ( 'PortalFace_L_R', C % IGNORABILITY + 2 )
      call C % PortalFace_R_L % Show &
             ( 'PortalFace_R_L', C % IGNORABILITY + 2 )
      call C % PortalEdge_LL_RR % Show &
             ( 'PortalEdge_LL_RR', C % IGNORABILITY + 2 )
      call C % PortalEdge_RR_LL % Show &
             ( 'PortalEdge_RR_LL', C % IGNORABILITY + 2 )
      call C % PortalEdge_LR_RL % Show &
             ( 'PortalEdge_LR_RL', C % IGNORABILITY + 2 )
      call C % PortalEdge_RL_LR % Show &
             ( 'PortalEdge_RL_LR', C % IGNORABILITY + 2 )
    end if !-- Distributed

    do iD = 1, nD
      call Show ( iD, 'iDimension', C % IGNORABILITY + 3 )
      call Show ( C % Edge ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Edge', C % IGNORABILITY + 3 )
      call Show ( C % Width ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Width', C % IGNORABILITY + 3 )
      call Show ( C % Center ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Center', C % IGNORABILITY + 3 )
    end do !-- iD

    end associate !-- nD

  end subroutine Show_C


  subroutine SetFieldPointer_1D_3D ( C, Field_1D, Field_3D )

    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Field_1D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      Field_3D

    Field_3D &
      ( C % iaFirst ( 1 ) : C % iaLast ( 1 ), &
        C % iaFirst ( 2 ) : C % iaLast ( 2 ), &
        C % iaFirst ( 3 ) : C % iaLast ( 3 ) ) &
          => Field_1D
    
  end subroutine SetFieldPointer_1D_3D

  
  subroutine SetFieldPointer_2D_4D ( C, Field_2D, Field_4D )

    class ( Chart_GS_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( :, : ), intent ( in ), contiguous, target :: &
      Field_2D
    real ( KDR ), dimension ( :, :, :, : ), intent ( out ), pointer :: &
      Field_4D
      
    integer ( KDI ) :: &
      nFields
      
    nFields = size ( Field_2D, dim = 2 )

    Field_4D &
      ( C % iaFirst ( 1 ) : C % iaLast ( 1 ), &
        C % iaFirst ( 2 ) : C % iaLast ( 2 ), &
        C % iaFirst ( 3 ) : C % iaLast ( 3 ), &
        1 : nFields ) &
          => Field_2D 
    
  end subroutine SetFieldPointer_2D_4D


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_Form ), intent ( inout ) :: &
      C

    nullify ( C % Communicator )

    if ( allocated ( C % PortalEdge_RL_LR ) ) &
      deallocate ( C % PortalEdge_RL_LR )
    if ( allocated ( C % PortalEdge_LR_RL ) ) &
      deallocate ( C % PortalEdge_LR_RL )
    if ( allocated ( C % PortalEdge_RR_LL ) ) &
      deallocate ( C % PortalEdge_RR_LL )
    if ( allocated ( C % PortalEdge_LL_RR ) ) &
      deallocate ( C % PortalEdge_LL_RR )
    if ( allocated ( C % PortalFace_R_L ) ) &
      deallocate ( C % PortalFace_R_L )
    if ( allocated ( C % PortalFace_L_R ) ) &
      deallocate ( C % PortalFace_L_R )

    if ( allocated ( C % ProperCell ) ) &
      deallocate ( C % ProperCell )

    if ( allocated ( C % Center ) ) &
      deallocate ( C % Center )
    if ( allocated ( C % Width ) ) &
      deallocate ( C % Width )
    if ( allocated ( C % Edge ) ) &
      deallocate ( C % Edge )

  end subroutine Finalize


  subroutine SetCoordinateMetadata &
               ( C, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
                 RatioOption, ScaleOption, nEqualOption )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), intent ( in ), optional :: &
      nEqualOption

    associate ( nD => C % nDimensions )

    C % MinCoordinate = 0.0_KDR
    if ( present ( MinCoordinateOption ) ) &
      C % MinCoordinate ( : nD ) = MinCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MinCoordinate ( : nD ), 'MinCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    C % MaxCoordinate = 0.0_KDR
    C % MaxCoordinate ( : nD ) = 1.0_KDR
    if ( present ( MaxCoordinateOption ) ) &
      C % MaxCoordinate ( : nD ) = MaxCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MaxCoordinate ( : nD ), 'MaxCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    C % Spacing = ''
    C % Spacing ( : nD ) = 'EQUAL'
    if ( present ( SpacingOption ) ) &
      C % Spacing ( : nD ) = SpacingOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Spacing ( : nD ), 'Spacing' )

    C % Ratio = 0.0_KDR
    if ( present ( RatioOption ) ) &
      C % Ratio ( : nD ) = RatioOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Ratio ( : nD ), 'Ratio' )

    C % Scale = 0.0_KDR
    if ( present ( ScaleOption ) ) &
      C % Scale ( : nD ) = ScaleOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Scale ( : nD ), 'Scale' )

    C % nEqual = 0
    if ( present ( nEqualOption ) ) &
      C % nEqual = nEqualOption

    end associate !-- nD

  end subroutine SetCoordinateMetadata


  subroutine SetCells ( C, nCellsOption, nGhostLayersOption )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption

    associate ( nD => C % nDimensions )

    C % nCells = 1
    C % nCells ( : nD ) = 32
    if ( present ( nCellsOption ) ) &
      C % nCells ( : nD ) = nCellsOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % nCells ( : nD ), 'nCells' )

    C % nGhostLayers = 0
    C % nGhostLayers ( : nD ) = 2
    if ( present ( nGhostLayersOption ) ) &
      C % nGhostLayers ( : nD ) = nGhostLayersOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % nGhostLayers ( : nD ), 'nGhostLayers' )

    end associate !-- nD

  end subroutine SetCells


  subroutine SetDecomposition &
                ( C, CommunicatorOption, nBricksOption, &
                  nBricksCompatibleOption )
               
    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
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
      C % Distributed   =   .true.
      C % Communicator  =>  CommunicatorOption
    else
      C % Distributed  =  .false.
    end if !-- present Communicator 

    if ( C % Distributed ) then

      associate ( nD => C % nDimensions )

      SizeRoot  =  C % Communicator % Size ** ( 1.0_KDR / nD ) + 0.5_KDR

      C % nBricks = 1
      C % nBricks ( : nD ) = SizeRoot
      if ( present ( nBricksOption ) ) &
        C % nBricks  =  nBricksOption 
      call PROGRAM_HEADER % GetParameter ( C % nBricks ( : nD ), 'nBricks' )
    
      nBricksCompatible = C % nBricks
      if ( present ( nBricksCompatibleOption ) ) &
        nBricksCompatible  =  nBricksCompatibleOption 
      call PROGRAM_HEADER % GetParameter &
             ( nBricksCompatible ( : nD ), 'nBricksCompatible' )

      if ( any ( nBricksCompatible /= C % nBricks ) ) then
        call Show ( 'nBricksCompatible /= nBricks', CONSOLE % INFO_1 )
        call Show ( C % nBricks, 'nBricks', CONSOLE % INFO_1 )
        call Show ( nBricksCompatible, 'nBricksCompatible', CONSOLE % INFO_1 )
      end if

      if ( product ( C % nBricks ) /= C % Communicator % Size ) then
        call Show ( 'The total number of bricks must equal ' &
                    // 'the number of MPI processes', CONSOLE % ERROR )
        call Show ( C % Communicator % Size, 'nProcesses', CONSOLE % ERROR )
        call Show ( C % nBricks ( 1 : nD ), 'nBricks', CONSOLE % ERROR )
        call Show ( product ( C % nBricks ), 'product ( nBricks )', &
                    CONSOLE % ERROR )
        call Show ( 'Chart_GS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetDecomposition', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      do iD = 1, nD
        if ( mod ( C % nCells ( iD ), nBricksCompatible ( iD ) ) /= 0 ) then
          call Show ( 'nBricksCompatible in each dimension must divide ' &
                      // 'evenly into nCells in each dimension', &
                      CONSOLE % WARNING )
          call Show ( iD, 'iDimension', CONSOLE % WARNING )
          call Show ( nBricksCompatible ( iD ), 'nBricksCompatible', &
                      CONSOLE % WARNING )
          call Show ( C % nCells ( iD ), 'nCells requested', CONSOLE % WARNING )
          C % nCells ( iD ) &
            =  ( C % nCells ( iD ) / nBricksCompatible ( iD ) ) &
               *  nBricksCompatible ( iD )
          call Show ( C % nCells ( iD ), 'nCells granted', CONSOLE % WARNING )
          call Show ( 'SetDecomposition', 'subroutine', CONSOLE % WARNING )
          call Show ( 'Chart_GS__Form', 'module', CONSOLE % WARNING )
        end if
      end do  !-- iD
    
      C % nCellsBrick &
        = C % nCells / C % nBricks
      C % iaBrick &
        = BrickIndex ( C % nBricks, C % nCells, C % Communicator % Rank )

      end associate !-- nD

      call SetPortals ( C )

      call SetCellsLocal ( C, C % nCellsBrick )

    else  !-- not Distributed

      C % nBricks      =  [ 1, 1, 1 ]
      C % nCellsBrick  =  C % nCells
      C % iaBrick      =  [ 1, 1, 1 ]

      call SetCellsLocal ( C, C % nCells )

    end if  !-- Distributed

    call SetProperCells ( C )

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


  subroutine SetPortals ( C )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C

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
      ( nB   =>  C % nBricks, &
        nD   =>  C % nDimensions, &
        nGL  =>  C % nGhostLayers, &
        nCB  =>  C % nCellsBrick ) 

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

      iaB  =  C % iaBrick

      !-- Left brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      Source_L_R ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_R_L ( iD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- Right brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
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
        =  C % nGhostLayers ( iD )  *  C % nGhostLayers ( jD )  &
           *  C % nCellsBrick ( kD )

      iaB  =  C % iaBrick

      !-- LeftLeft brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_LL_RR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightRight brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_RR_LL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LL_RR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- LeftRight brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ) - 1 + nB ( iD ) - 1, nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ), nB ( jD ) ) + 1
      Source_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_RL_LR ( KD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

      !-- RightLeft brick

      iaB ( iD )  &
        =  mod ( C % iaBrick ( iD ), nB ( iD ) ) + 1
      iaB ( jD )  &
        =  mod ( C % iaBrick ( jD ) - 1 + nB ( jD ) - 1, nB ( jD ) ) + 1
      Source_RL_LR ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )
      Target_LR_RL ( kD )  =  Process ( iaB ( 1 ), iaB ( 2 ), iaB ( 3 ) )

    end do !-- kD


    !-- Set face portals
    
    allocate ( C % PortalFace_L_R )
    associate ( PFLR => C % PortalFace_L_R )
    call PFLR % Initialize ( Source_L_R, Target_L_R, nCellsFace, nCellsFace )
    end associate !-- PFLR

    allocate ( C % PortalFace_R_L )
    associate ( PFRL => C % PortalFace_R_L )
    call PFRL % Initialize ( Source_R_L, Target_R_L, nCellsFace, nCellsFace )
    end associate !-- PFRL


    !-- Set edge portals
    
    allocate ( C % PortalEdge_LL_RR )
    associate ( PELLRR => C % PortalEdge_LL_RR )
    call PELLRR % Initialize &
           ( pack ( Source_LL_RR, Source_LL_RR >= 0 ), &
             pack ( Target_LL_RR, Target_LL_RR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELLRR

    allocate ( C % PortalEdge_RR_LL )
    associate ( PERRLL => C % PortalEdge_RR_LL )
    call PERRLL % Initialize &
           ( pack ( Source_RR_LL, Source_RR_LL >= 0 ), &
             pack ( Target_RR_LL, Target_RR_LL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PERRLL

    allocate ( C % PortalEdge_LR_RL )
    associate ( PELRRL => C % PortalEdge_LR_RL )
    call PELRRL % Initialize &
           ( pack ( Source_LR_RL, Source_LR_RL >= 0 ), &
             pack ( Target_LR_RL, Target_LR_RL >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELRRL

    allocate ( C % PortalEdge_RL_LR )
    associate ( PERLLR => C % PortalEdge_RL_LR )
    call PERLLR % Initialize &
           ( pack ( Source_RL_LR, Source_RL_LR >= 0 ), &
             pack ( Target_RL_LR, Target_RL_LR >= 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ), &
             pack ( nCellsEdge, nCellsEdge > 0 ) ) 
    end associate !-- PELRRL


    !-- Cleanup

    end associate !-- nB, etc.
    
  end subroutine SetPortals


  subroutine SetCellsLocal ( C, nCellsLocal )

    class ( Chart_GS_Form ), intent ( inout ) :: &
      C
    integer, dimension ( : ), intent ( in ) :: &
      nCellsLocal

    C % iaFirst = 1
    C % iaLast = 1
    C % iaFirst ( 1 ) = 1 - C % nGhostLayers ( 1 )
    C % iaLast  ( 1 ) = nCellsLocal ( 1 ) + C % nGhostLayers ( 1 )
    if ( C % nDimensions > 1 ) then
      C % iaFirst ( 2 ) = 1 - C % nGhostLayers ( 2 )
      C % iaLast  ( 2 ) = nCellsLocal ( 2 ) + C % nGhostLayers ( 2 )
    end if
    if ( C % nDimensions > 2 ) then
      C % iaFirst ( 3 ) = 1 - C % nGhostLayers ( 3 )
      C % iaLast  ( 3 ) = nCellsLocal ( 3 ) + C % nGhostLayers ( 3 )
    end if

    C % nCellsLocal   =  product ( nCellsLocal  +  2 * C % nGhostLayers )
    C % nCellsProper  =  product ( nCellsLocal )
    C % nCellsGhost   =  C % nCellsLocal  -  C % nCellsProper

  end subroutine SetCellsLocal


  subroutine SetProperCells ( C )

    class ( Chart_GS_Form ), intent ( inout ), target :: &
      C

    integer ( KDI ) :: &
      iC, jC, kC, &
      iV
    logical ( KDL ), dimension ( :, :, : ), pointer :: &
      PC

    allocate ( C % ProperCell ( C % nCellsLocal ) )
    call Clear ( C % ProperCell )

    associate &
      ( iaF  =>  C % iaFirst, &
        iaL  =>  C % iaLast, &
        nGL  =>  C % nGhostLayers )

    PC ( iaF ( 1 ) : iaL ( 1 ), &
         iaF ( 2 ) : iaL ( 2 ), &
         iaF ( 3 ) : iaL ( 3 ) )  &
      =>  C % ProperCell

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

    type ( QuantityForm ), intent ( in ) :: &
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


end module Chart_GS__Form
