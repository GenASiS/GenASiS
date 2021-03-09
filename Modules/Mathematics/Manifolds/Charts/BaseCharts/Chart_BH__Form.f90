module Chart_BH__Form

  !-- Chart_BaseHeader_Form

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none
  private

  type, public, extends ( Chart_H_Form ) :: Chart_BH_Form
    integer ( KDI ) :: &
      nValues = 0, &
      nEqual  = 0
    integer ( KDI ), dimension ( : ), pointer :: &
      iaFirst      => null ( ), &
      iaLast       => null ( ), &
      nCells       => null ( ), &
      nGhostLayers => null ( )
    integer ( KDI ), dimension ( : ), pointer :: &
      iaBrick      => null ( ), &
      nBricks      => null ( ), &
      nCellsBrick  => null ( )
    real ( KDR ), dimension ( : ), pointer :: &
      MinCoordinate => null ( ), &
      MaxCoordinate => null ( ), &
      Ratio => null ( ), &
      Scale => null ( )
    class ( Real_1D_Form ), dimension ( : ), pointer :: &
      Edge, &
      Width, &
      Center
    character ( LDL ), dimension ( : ), pointer :: &
      Spacing => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic_BH
    generic, public :: &
      Initialize_BH => InitializeBasic_BH
    procedure, private, pass :: &
      Show_C
    final :: &
      Finalize
    procedure, public, pass :: &
      SetCoordinateData
  end type Chart_BH_Form

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = MANIFOLD % MAX_DIMENSIONS

    private :: &
      SetCoordinateMetadata, &
      SetCells, &
      SetDecomposition

      private :: &
        BrickIndex, &
        SetFirstLast, &
        SetEdgeEqual, &
        ComputeGeometricRatio, &
        SetEdgeGeometric, &
        SetEdgeCompactified, &
        SetEdgeProportional

        private :: &
          ZeroGeometricRatio


contains


  subroutine InitializeBasic_BH &
               ( C, M, IsPeriodic, iChart, CommunicatorOption, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nBricksOption, &
                 nBricksCompatibleOption, nDimensionsOption, nEqualOption )

    class ( Chart_BH_Form ), intent ( inout ) :: &
      C
    class ( Manifold_H_Form ), intent ( in ), target :: &
      M
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
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
      nDimensionsOption, &
      nEqualOption

    integer ( KDI ) :: &
      iD  !-- iDimension

    call C % Chart_H_Form % Initialize_H &
           ( M, IsPeriodic, iChart, CommunicatorOption, &
             CoordinateLabelOption, CoordinateSystemOption, &
             CoordinateUnitOption, nDimensionsOption )

    call SetCoordinateMetadata &
           ( C, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
             RatioOption, ScaleOption, nEqualOption )

    call SetCells ( C, nCellsOption, nGhostLayersOption )

    call SetDecomposition &
           ( C, M, CommunicatorOption, nBricksOption, nBricksCompatibleOption )

    do iD = 1, C % nDimensions
      call SetCoordinateData ( C, iD )
    end do !-- iD

  end subroutine InitializeBasic_BH


  subroutine Show_C ( C )

    class ( Chart_BH_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iD  !-- iDimension

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

    call Show ( C % nCells ( : nD ), 'nCells', C % IGNORABILITY )
    call Show ( C % nGhostLayers ( : nD ), 'nGhostLayers', C % IGNORABILITY )

    call Show ( C % nValues, 'nValues', C % IGNORABILITY )

    call Show ( C % iaFirst ( : nD ), 'iaFirst', C % IGNORABILITY )
    call Show ( C % iaLast  ( : nD ), 'iaLast',  C % IGNORABILITY )

    if ( C % IsDistributed ) then
      call Show ( C % iaBrick ( : nD ), 'iaBrick', C % IGNORABILITY )
      call Show ( C % nBricks ( : nD ), 'nBricks', C % IGNORABILITY )
      call Show ( C % nCellsBrick ( : nD ), 'nCellsBrick', C % IGNORABILITY )
    end if !-- IsDistributed

    do iD = 1, nD
      call Show ( iD, 'iDimension' )
      call Show ( C % Edge ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Edge', C % IGNORABILITY + 1 )
      call Show ( C % Width ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Width', C % IGNORABILITY + 1 )
      call Show ( C % Center ( iD ) % Value, C % CoordinateUnit ( iD ), &
                  'Center', C % IGNORABILITY + 1 )
    end do !-- iD

    end associate !-- nD

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_BH_Form ), intent ( inout ) :: &
      C

    if ( .not. associated ( C % Name ) ) &
      return
    if ( C % Name == '' ) &
      return

    if ( C % AllocatedValues ) then

      deallocate ( C % Spacing )
      deallocate ( C % Center )
      deallocate ( C % Width )
      deallocate ( C % Edge )
      deallocate ( C % Scale )
      deallocate ( C % Ratio )
      deallocate ( C % MaxCoordinate )
      deallocate ( C % MinCoordinate )

      if ( C % IsDistributed ) then
        deallocate ( C % nCellsBrick )
        deallocate ( C % nBricks )
        deallocate ( C % iaBrick )
      end if

      deallocate ( C % nGhostLayers )
      deallocate ( C % nCells )
      deallocate ( C % iaLast )
      deallocate ( C % iaFirst )

    else

      nullify ( C % Spacing )
      nullify ( C % Center )
      nullify ( C % Width )
      nullify ( C % Edge )
      nullify ( C % Scale )
      nullify ( C % Ratio )
      nullify ( C % MaxCoordinate )
      nullify ( C % MinCoordinate )

      nullify ( C % nCellsBrick )
      nullify ( C % nBricks )
      nullify ( C % iaBrick )

      nullify ( C % nGhostLayers )
      nullify ( C % nCells )
      nullify ( C % iaLast )
      nullify ( C % iaFirst )

    end if !-- AllocatedValues

  end subroutine Finalize


  subroutine SetCoordinateData ( C, iD, EdgeValueOption )

    class ( Chart_BH_Form ), intent ( inout ) :: &
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

    if ( .not. C % AllocatedValues ) &
      return

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
        call SetEdgeEqual &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), nC )
      case ( 'GEOMETRIC' )
        if ( C % Scale ( iD ) > 0.0_KDR ) &
          call ComputeGeometricRatio &
                 ( C % CoordinateUnit ( iD ), C % MinCoordinate ( iD ), &
                   C % MaxCoordinate ( iD ), C % Scale ( iD ), nC, &
                   C % Ratio ( iD ) )
        call SetEdgeGeometric &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % MaxCoordinate ( iD ), &
                 C % Ratio ( iD ), nC )
      case ( 'COMPACTIFIED' )
        call SetEdgeCompactified &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % Scale ( iD ), nC )
        C % MinCoordinate ( iD )  =  C % Edge ( iD ) % Value ( 1 )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case ( 'PROPORTIONAL' )
        call SetEdgeProportional &
               ( C % Edge ( iD ) % Value ( 1 : nC + 1 ), &
                 C % MinCoordinate ( iD ), C % Ratio ( iD ), &
                 C % Scale ( iD ), nC, C % nEqual )
        C % MaxCoordinate ( iD )  =  C % Edge ( iD ) % Value ( nC + 1 )
      case default
        call Show ( 'Spacing not recognized', CONSOLE % ERROR )
        call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetGeometryCell', 'subroutine', CONSOLE % ERROR )
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

  end subroutine SetCoordinateData


  subroutine SetCoordinateMetadata &
               ( C, SpacingOption, MinCoordinateOption, MaxCoordinateOption, &
                 RatioOption, ScaleOption, nEqualOption )

    class ( Chart_BH_Form ), intent ( inout ) :: &
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

    allocate ( C % MinCoordinate ( MAX_DIMENSIONS ) )
    C % MinCoordinate = 0.0_KDR
    if ( present ( MinCoordinateOption ) ) &
      C % MinCoordinate ( : nD ) = MinCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MinCoordinate ( : nD ), 'MinCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    allocate ( C % MaxCoordinate ( MAX_DIMENSIONS ) )
    C % MaxCoordinate = 0.0_KDR
    C % MaxCoordinate ( : nD ) = 1.0_KDR
    if ( present ( MaxCoordinateOption ) ) &
      C % MaxCoordinate ( : nD ) = MaxCoordinateOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % MaxCoordinate ( : nD ), 'MaxCoordinate', &
             InputUnitOption = C % CoordinateUnit ( : nD ) )

    allocate ( C % Spacing ( MAX_DIMENSIONS ) )
    C % Spacing = ''
    C % Spacing ( : nD ) = 'EQUAL'
    if ( present ( SpacingOption ) ) &
      C % Spacing ( : nD ) = SpacingOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Spacing ( : nD ), 'Spacing' )

    allocate ( C % Ratio ( MAX_DIMENSIONS ) )
    C % Ratio = 0.0_KDR
    if ( present ( RatioOption ) ) &
      C % Ratio ( : nD ) = RatioOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Ratio ( : nD ), 'Ratio' )

    allocate ( C % Scale ( MAX_DIMENSIONS ) )
    C % Scale = 0.0_KDR
    if ( present ( ScaleOption ) ) &
      C % Scale ( : nD ) = ScaleOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % Scale ( : nD ), 'Scale' )

    C % nEqual = 0
    if ( present ( nEqualOption ) ) &
      C % nEqual = nEqualOption

    allocate ( C % Edge   ( MAX_DIMENSIONS ) )
    allocate ( C % Width  ( MAX_DIMENSIONS ) )
    allocate ( C % Center ( MAX_DIMENSIONS ) )

    end associate !-- nD

  end subroutine SetCoordinateMetadata


  subroutine SetCells ( C, nCellsOption, nGhostLayersOption )

    class ( Chart_BH_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption

    associate ( nD => C % nDimensions )

    allocate ( C % nCells ( MAX_DIMENSIONS ) )
    C % nCells = 1
    C % nCells ( : nD ) = 32
    if ( present ( nCellsOption ) ) &
      C % nCells ( : nD ) = nCellsOption ( : nD )
    call PROGRAM_HEADER % GetParameter ( C % nCells ( : nD ), 'nCells' )

    allocate ( C % nGhostLayers ( MAX_DIMENSIONS ) )
    C % nGhostLayers = 0
    C % nGhostLayers ( : nD ) = 2
    if ( present ( nGhostLayersOption ) ) &
      C % nGhostLayers ( : nD ) = nGhostLayersOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % nGhostLayers ( : nD ), 'nGhostLayers' )

    end associate !-- nD

  end subroutine SetCells


  subroutine SetDecomposition &
                ( C, M, CommunicatorOption, nBricksOption, &
                  nBricksCompatibleOption )
               
    class ( Chart_BH_Form ), intent ( inout ) :: &
      C
    class ( Manifold_H_Form ), intent ( in ), target :: &
      M
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
      C % IsDistributed  =   .true.
      C % Communicator   =>  CommunicatorOption
    else
      C % IsDistributed  =   M % IsDistributed
      C % Communicator   =>  M % Communicator
    end if !-- present Communicator 

    if ( C % IsDistributed ) then

      allocate ( C % iaBrick ( MAX_DIMENSIONS ) )
      allocate ( C % nBricks ( MAX_DIMENSIONS ) )
      allocate ( C % nCellsBrick ( MAX_DIMENSIONS ) )

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
        call Show ( 'ChartHeader_Form', 'module', CONSOLE % ERROR )
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
          C % nCells ( iD ) = ( C % nCells ( iD ) / nBricksCompatible ( iD ) ) &
                          * nBricksCompatible ( iD )
          call Show ( C % nCells ( iD ), 'nCells granted', CONSOLE % WARNING )
          call Show ( 'SetDecomposition', 'subroutine', CONSOLE % WARNING )
          call Show ( 'ChartHeader_Form', 'module', CONSOLE % WARNING )
        end if
      end do
    
      C % nCellsBrick &
        = C % nCells / C % nBricks
      C % iaBrick &
        = BrickIndex ( C % nBricks, C % nCells, C % Communicator % Rank )

      end associate !-- nD

      C % nValues  &
        =  product ( C % nCellsBrick  +  2 * C % nGhostLayers )

      call SetFirstLast ( C, C % nCellsBrick )

    else  !-- not Distributed

      C % nValues  &
        =  product ( C % nCells  +  2 * C % nGhostLayers )

      call SetFirstLast ( C, C % nCells )

    end if  !-- IsDistributed

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


  subroutine SetFirstLast ( C, nCellsLocal )

    class ( Chart_BH_Form ), intent ( inout ) :: &
      C 
    integer, dimension ( : ), intent ( in ) :: &
      nCellsLocal

    allocate ( C % iaFirst ( MAX_DIMENSIONS ) )
    allocate ( C % iaLast ( MAX_DIMENSIONS ) )
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

  end subroutine SetFirstLast


  subroutine SetEdgeEqual ( Edge, MinCoordinate, MaxCoordinate, nC )

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

  end subroutine SetEdgeEqual


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


  subroutine SetEdgeGeometric &
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

  end subroutine SetEdgeGeometric


  subroutine SetEdgeCompactified ( Edge, Scale, nC )

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

  end subroutine SetEdgeCompactified


  subroutine SetEdgeProportional &
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
      call SetEdgeEqual ( Edge, MinCoordinate, Scale, nEqual )
    end if
    
    Width  =  Ratio  *  Edge ( nEqual + 1 )

    do iC = nEqual + 2, nC + 1
      Edge ( iC )  =  Edge ( iC - 1 )  +  Width
      Width        =  Ratio  *  Edge ( iC )
    end do

  end subroutine SetEdgeProportional


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


end module Chart_BH__Form
