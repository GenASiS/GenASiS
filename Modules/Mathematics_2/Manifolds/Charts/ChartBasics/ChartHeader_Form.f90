!-- ChartHeader_Form contains the most basic functionality characterizing a
!   coordinate chart.

module ChartHeader_Form

  use Basics
  use ManifoldBasics

  implicit none
  private

  type, public :: ChartHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iChart       = 0, &
      nDimensions  = 0, &
      nEqual       = 0, &
      nFields      = 0
    real ( KDR ), dimension ( : ), pointer :: &
      MinCoordinate => null ( ), &
      MaxCoordinate => null ( ), &
      Ratio => null ( ), &
      Scale => null ( )
    type ( MeasuredValueForm ), dimension ( : ), pointer :: &
      CoordinateUnit => null ( )
    logical ( KDL ) :: &
      IsDistributed = .false., &
      AllocatedValues = .false.
    logical ( KDL ), dimension ( : ), pointer :: &
      IsPeriodic => null ( )
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( )
    character ( LDL ), pointer :: &
      CoordinateSystem => null ( )
    character ( LDL ), dimension ( : ), pointer :: &
      CoordinateLabel => null ( ), &   
      Spacing => null ( )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    class ( ManifoldHeaderForm ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    generic, public :: &
      Initialize => InitializeBasic
    procedure, private, pass :: &
      Show_CH
    generic, public :: &
      Show => Show_CH
    final :: &
      Finalize
    procedure, public, pass ( C ) :: &
      SetBrick
  end type ChartHeaderForm

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = MANIFOLD % MAX_DIMENSIONS

    private :: &
      BrickIndex


contains


  subroutine InitializeBasic &
               ( C, M, IsPeriodic, iChart, CommunicatorOption, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nDimensionsOption, nEqualOption )

    class ( ChartHeaderForm ), intent ( inout ) :: &
      C
    class ( ManifoldHeaderForm ), intent ( in ), target :: &
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
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    character ( 2 ) :: &
      ChartNumber

    C % IGNORABILITY  =  M % IGNORABILITY
    C % Manifold  =>  M
    C % iChart  =  iChart

    C % AllocatedValues = .true.

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart'
    end if

    allocate ( C % Name )
    write ( ChartNumber, fmt = '(i2.2)' ) iChart
    C % Name = 'Chart_' // ChartNumber // '_' // trim ( M % Name ) 

    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    if ( present ( CommunicatorOption ) ) then
      C % IsDistributed  =   .true.
      C % Communicator   =>  CommunicatorOption
    else
      C % IsDistributed  =   M % IsDistributed
      C % Communicator   =>  M % Communicator
    end if !-- present Communicator 

    if ( present ( nDimensionsOption ) ) then
      C % nDimensions  =  nDimensionsOption
    else
      C % nDimensions  =  M % nDimensions
    end if

    associate ( nD => C % nDimensions )

    allocate ( C % IsPeriodic ( MAX_DIMENSIONS ) )
    C % IsPeriodic = .false.
    C % IsPeriodic ( : nD ) = IsPeriodic ( : nD )

    allocate ( C % CoordinateUnit ( MAX_DIMENSIONS ) )
    C % CoordinateUnit = [ UNIT % IDENTITY, UNIT % IDENTITY, UNIT % IDENTITY ]
    if ( present ( CoordinateUnitOption ) ) &
      C % CoordinateUnit ( : nD ) = CoordinateUnitOption ( : nD )

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

    allocate ( C % CoordinateSystem )
    C % CoordinateSystem = 'RECTANGULAR'
    if ( present ( CoordinateSystemOption ) ) &
      C % CoordinateSystem = CoordinateSystemOption
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateSystem, 'CoordinateSystem' )

    allocate ( C % CoordinateLabel ( MAX_DIMENSIONS ) )
    C % CoordinateLabel  =  [ 'X', 'Y', 'Z' ]
    select case ( trim ( C % CoordinateSystem ) )
    case ( 'CYLINDRICAL' )
      C % CoordinateLabel  =  [ 'R_Perp', 'Z     ', 'Phi   ' ]
    case ( 'SPHERICAL' )
      C % CoordinateLabel  =  [ 'R    ', 'Theta', 'Phi  ' ]
    end select !-- CoordinateSystem
    if ( present ( CoordinateLabelOption ) ) &
      C % CoordinateLabel ( : nD ) = CoordinateLabelOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateLabel ( : nD ), 'CoordinateLabel' )

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

    end associate !-- nD

  end subroutine InitializeBasic


  subroutine Show_CH ( C )

    class ( ChartHeaderForm ), intent ( in ) :: &
      C

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    associate ( nD => C % nDimensions )

    call Show ( C % IsDistributed, 'IsDistributed', C % IGNORABILITY )
    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )
    call Show ( C % nFields, 'nFields', C % IGNORABILITY )

    call Show ( C % IsPeriodic ( : nD ), 'IsPeriodic', C % IGNORABILITY )

    call Show ( C % MinCoordinate ( : nD ), C % CoordinateUnit ( : nD ), &
                'MinCoordinate', C % IGNORABILITY )
    call Show ( C % MaxCoordinate ( : nD ), C % CoordinateUnit ( : nD ), &
                'MaxCoordinate', C % IGNORABILITY )

    call Show ( C % CoordinateSystem, 'CoordinateSystem', C % IGNORABILITY )

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

    end associate !-- nD

  end subroutine Show_CH


  impure elemental subroutine Finalize ( C )

    type ( ChartHeaderForm ), intent ( inout ) :: &
      C

    nullify ( C % Manifold )
    nullify ( C % Communicator )

    if ( .not. associated ( C % Name ) ) &
      return
    if ( C % Name == '' ) &
      return

    if ( C % AllocatedValues ) then
      deallocate ( C % Scale )
      deallocate ( C % Ratio )
      deallocate ( C % Spacing )
      deallocate ( C % CoordinateLabel )
      deallocate ( C % CoordinateSystem )
      deallocate ( C % MaxCoordinate )
      deallocate ( C % MinCoordinate )
      deallocate ( C % CoordinateUnit )
      deallocate ( C % IsPeriodic )
    else
      nullify ( C % Scale )
      nullify ( C % Ratio )
      nullify ( C % Spacing )
      nullify ( C % CoordinateLabel )
      nullify ( C % CoordinateSystem )
      nullify ( C % MaxCoordinate )
      nullify ( C % MinCoordinate )
      nullify ( C % CoordinateUnit )
      nullify ( C % IsPeriodic )
    end if !-- AllocatedValues

    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    if ( C % AllocatedValues ) then
      deallocate ( C % Name )
      deallocate ( C % Type )
    else
      nullify ( C % Name )
      nullify ( C % Type )
    end if !-- AllocatedValues

  end subroutine Finalize


  subroutine SetBrick &
               ( nCells, C, Communicator, nCellsBrick, nBricks, iaBrick, &
                 nBricksOption, nBricksCompatibleOption )

    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      nCells
    class ( ChartHeaderForm ), intent ( in ) :: &
      C
    type ( CommunicatorForm ), intent ( in ) :: &
      Communicator
    integer ( KDI ), dimension ( : ), intent ( out ) :: &
      nCellsBrick, &
      nBricks, &
      iaBrick
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nBricksOption, &
      nBricksCompatibleOption

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      SizeRoot
    integer ( KDI ), dimension ( size ( nBricks ) ) :: &
      nBricksCompatible

    associate ( nD => C % nDimensions )

    SizeRoot  =  Communicator % Size ** ( 1.0_KDR / nD ) + 0.5_KDR

    nBricks = 1
    nBricks ( : nD ) = SizeRoot
    if ( present ( nBricksOption ) ) &
      nBricks  =  nBricksOption 
    call PROGRAM_HEADER % GetParameter ( nBricks ( : nD ), 'nBricks' )
    
    nBricksCompatible = nBricks
    if ( present ( nBricksCompatibleOption ) ) &
      nBricksCompatible  =  nBricksCompatibleOption 
    call PROGRAM_HEADER % GetParameter &
           ( nBricksCompatible ( : nD ), 'nBricksCompatible' )

    if ( any ( nBricksCompatible /= nBricks ) ) then
      call Show ( 'nBricksCompatible /= nBricks', CONSOLE % INFO_1 )
      call Show ( nBricks, 'nBricks', CONSOLE % INFO_1 )
      call Show ( nBricksCompatible, 'nBricksCompatible', CONSOLE % INFO_1 )
    end if

    if ( product ( nBricks ) /= Communicator % Size ) then
      call Show ( 'The total number of bricks must equal ' &
                  // 'the number of MPI processes', CONSOLE % ERROR )
      call Show ( Communicator % Size, 'nProcesses', CONSOLE % ERROR )
      call Show ( nBricks ( 1 : nD ), 'nBricks', CONSOLE % ERROR )
      call Show ( product ( nBricks ), 'product ( nBricks )', CONSOLE % ERROR )
      call Show ( 'Chart_Template', 'module', CONSOLE % ERROR )
      call Show ( 'SetBrick', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if

    do iD = 1, nD
      if ( mod ( nCells ( iD ), nBricksCompatible ( iD ) ) /= 0 ) then
        call Show ( 'nBricksCompatible in each dimension must divide evenly ' &
                    // 'into nCells in each dimension', CONSOLE % WARNING )
        call Show ( iD, 'iDimension', CONSOLE % WARNING )
        call Show ( nBricksCompatible ( iD ), 'nBricksCompatible', &
                    CONSOLE % WARNING )
        call Show ( nCells ( iD ), 'nCells requested', CONSOLE % WARNING )
        nCells ( iD ) = ( nCells ( iD ) / nBricksCompatible ( iD ) ) &
                        * nBricksCompatible ( iD )
        call Show ( nCells ( iD ), 'nCells granted', CONSOLE % WARNING )
        call Show ( 'SetBrick', 'subroutine', CONSOLE % WARNING )
        call Show ( 'Chart_Template', 'module', CONSOLE % WARNING )
      end if
    end do
    
    nCellsBrick = nCells / nBricks
    iaBrick = BrickIndex ( nBricks, nCells, Communicator % Rank )

    end associate !-- nD

  end subroutine SetBrick


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


end module ChartHeader_Form
