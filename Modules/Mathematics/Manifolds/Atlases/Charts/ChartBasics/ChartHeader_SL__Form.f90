!-- ChartHeader_SL handles metadata of a single-level Chart.

module ChartHeader_SL__Form

  !-- ChartHeaderSingleLevelForm

  use Basics
  use AtlasBasics
  use Chart_Template
  use FieldChart_Template
  use Field_CSL__Template

  implicit none
  private

  type, public, extends ( ChartTemplate ) :: ChartHeader_SL_Form
    integer ( KDI ) :: &
      iFieldGeometry = 0, &
      nProperCells   = 0, &
      nGhostCells    = 0
    integer ( KDI ), dimension ( : ), pointer :: &
      iaFirst      => null ( ), &
      iaLast       => null ( ), &
      iaBrick      => null ( ), &
      nCells       => null ( ), &
      nCellsBrick  => null ( ), &
      nBricks      => null ( ), &
      nGhostLayers => null ( )
    logical ( KDL ), dimension ( : ), pointer :: &
      IsProperCell => null ( )
    type ( FieldChartPointer ), dimension ( : ), allocatable :: &
      Field
    procedure ( CS ), pointer, nopass :: &
      CoarsenSingularities => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeBasic, InitializeClone
    procedure, public, pass :: &
      AddField
    procedure, public, pass :: &
      SetVariablePointer
    procedure, private, pass :: &
      ShowHeader
    generic, public :: &
      Show => ShowHeader
    final :: &
      Finalize
  end type ChartHeader_SL_Form

  abstract interface
    subroutine CS ( S, iAngular )
      use Basics
      class ( StorageForm ), intent ( inout ) :: &
        S
      integer ( KDI ), intent ( in ) :: &
        iAngular
    end subroutine CS
  end interface

    private :: &
      SetChartCells, &
      SetChartFirstLast

    !-- Avoid conflict with Atlas dummy variable in Initialize
    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = ATLAS % MAX_DIMENSIONS, &
      MAX_FIELDS     = ATLAS % MAX_FIELDS

contains


  subroutine InitializeBasic &
               ( C, Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 IsDistributedOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

    class ( ChartHeader_SL_Form ), intent ( inout ) :: &
      C
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      Atlas
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    integer ( KDI ), intent ( in ) :: &
      iChart
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    logical ( KDL ), intent ( in ), optional :: &
      IsDistributedOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    integer ( KDI ) :: &
      iC, jC, kC, &
      iV
    real ( KDR ), dimension ( : ), allocatable :: &
      IsProperCell
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      IPC

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart_SL'
    end if

    call C % InitializeTemplateBasic &
               ( Atlas, IsPeriodic, iChart, SpacingOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 IsDistributedOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nDimensionsOption, nEqualOption )

    allocate ( C % nCells ( MAX_DIMENSIONS ) )
    allocate ( C % nGhostLayers ( MAX_DIMENSIONS ) )
    call SetChartCells ( C, nCellsOption, nGhostLayersOption )

    if ( C % IsDistributed ) then
      allocate ( C % iaBrick ( MAX_DIMENSIONS ) )
      allocate ( C % nCellsBrick ( MAX_DIMENSIONS ) )
      allocate ( C % nBricks ( MAX_DIMENSIONS ) )
      call C % SetBrick &
             ( C % nCells, Atlas % Communicator, C % nCellsBrick, &
               C % nBricks, C % iaBrick )
      call SetChartFirstLast ( C, C % nCellsBrick )
    else
      call SetChartFirstLast ( C, C % nCells )
    end if !-- IsDistributed

    allocate ( C % IsProperCell ( C % nProperCells + C % nGhostCells ) )
    allocate ( IsProperCell ( C % nProperCells + C % nGhostCells ) )
    call Clear ( C % IsProperCell )
    call Clear ( IsProperCell )
    call C % SetVariablePointer ( IsProperCell, IPC )

    associate &
      ( lB => C % iaFirst  +  C % nGhostLayers, &
        uB => C % iaLast   -  C % nGhostLayers )  
    !$OMP parallel do private ( iC, jC, kC )
    do kC = lB ( 3 ), uB ( 3 )
      do jC = lB ( 2 ), uB ( 2 )
        do iC = lB ( 1 ), uB ( 1 )
          IPC ( iC, jC, kC ) = 1.0_KDR
        end do
      end do
    end do
    !$OMP end parallel do
    end associate !-- lB, etc.

    !$OMP parallel do private ( iV )
    do iV = 1, size ( IsProperCell )
      if ( IsProperCell ( iV ) > 0.0_KDR ) &
        C % IsProperCell ( iV ) = .true.
    end do
    !$OMP end parallel do
    nullify ( IPC )

    allocate ( C % Field ( MAX_FIELDS ) )

  end subroutine InitializeBasic


  subroutine InitializeClone ( C, C_Source )

    class ( ChartHeader_SL_Form ), intent ( inout ) :: &
      C
    class ( ChartHeader_SL_Form ), intent ( in ), target :: &
      C_Source

    call C % InitializeTemplateClone ( C_Source )

    C % iFieldGeometry =  0  !-- NOT COPIED!
    C % nProperCells   =  C_Source % nProperCells
    C % nGhostCells    =  C_Source % nGhostCells
    C % iaFirst        => C_Source % iaFirst
    C % iaLast         => C_Source % iaLast
    C % iaBrick        => C_Source % iaBrick
    C % nCells         => C_Source % nCells
    C % nCellsBrick    => C_Source % nCellsBrick
    C % nBricks        => C_Source % nBricks
    C % nGhostLayers   => C_Source % nGhostLayers
    C % IsProperCell   => C_Source % IsProperCell

    allocate ( C % Field ( MAX_FIELDS ) )

  end subroutine InitializeClone


  subroutine AddField ( C, FC )

    class ( ChartHeader_SL_Form ), intent ( inout ) :: &
      C
    class ( Field_CSL_Template ), intent ( in ), target :: &
      FC

    C % nFields = C % nFields + 1
    C % Field ( C % nFields ) % Pointer => FC

    call Show ( 'Adding a field', C % IGNORABILITY + 2 )
    call Show ( C % Name, 'Chart', C % IGNORABILITY + 2 )
    call Show ( FC % Name, 'Field', C % IGNORABILITY + 2 )
    call Show ( C % nFields, 'nFields', C % IGNORABILITY + 2 )

  end subroutine AddField


  subroutine SetVariablePointer ( C, Variable_1D, Variable_3D )

    class ( ChartHeader_SL_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      Variable_1D
    real ( KDR ), dimension ( :, :, : ), intent ( out ), pointer :: &
      Variable_3D

    Variable_3D &
      ( C % iaFirst ( 1 ) : C % iaLast ( 1 ), &
        C % iaFirst ( 2 ) : C % iaLast ( 2 ), &
        C % iaFirst ( 3 ) : C % iaLast ( 3 ) ) &
          => Variable_1D
    
  end subroutine SetVariablePointer


  subroutine ShowHeader ( C )

    class ( ChartHeader_SL_Form ), intent ( in ) :: &
      C

    call C % ShowTemplate ( )

    call Show ( C % nCells, 'nCells', C % IGNORABILITY )
    call Show ( C % nGhostLayers, 'nGhostLayers', C % IGNORABILITY )
    call Show ( C % nProperCells, 'nProperCells', C % IGNORABILITY )
    call Show ( C % nGhostCells, 'nGhostCells', C % IGNORABILITY )

    if ( C % IsDistributed ) then
      call Show ( C % iaBrick, 'iaBrick', C % IGNORABILITY )
      call Show ( C % nBricks, 'nBricks', C % IGNORABILITY )
      call Show ( C % nCellsBrick, 'nCellsBrick', C % IGNORABILITY )
    end if !-- IsDistributed

  end subroutine ShowHeader


  impure elemental subroutine Finalize ( C )

    type ( ChartHeader_SL_Form ), intent ( inout ) :: &
      C

    if ( allocated ( C % Field ) ) &
      deallocate ( C % Field )

    if ( C % AllocatedValues ) then
      deallocate ( C % iaFirst )
      deallocate ( C % iaLast )
      deallocate ( C % nCells )
      deallocate ( C % nGhostLayers )
      deallocate ( C % IsProperCell )
      if ( C % IsDistributed ) then
        allocate ( C % iaBrick ( MAX_DIMENSIONS ) )
        allocate ( C % nCellsBrick ( MAX_DIMENSIONS ) )
        allocate ( C % nBricks ( MAX_DIMENSIONS ) )
      end if
    else
      nullify ( C % iaFirst )
      nullify ( C % iaLast )
      nullify ( C % iaBrick )
      nullify ( C % nCells )
      nullify ( C % nCellsBrick )
      nullify ( C % nBricks )
      nullify ( C % nGhostLayers )
      nullify ( C % IsProperCell )
    end if !-- AllocatedValues

    call C % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetChartCells ( C, nCellsOption, nGhostLayersOption )

    class ( ChartHeader_SL_Form ), intent ( inout ) :: &
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

  end subroutine SetChartCells


  subroutine SetChartFirstLast ( C, nCellsLocal )

    class ( ChartHeader_SL_Form ), intent ( inout ) :: &
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

    C % nProperCells = product ( nCellsLocal )
    C % nGhostCells  &
      = product ( C % iaLast - ( C % iaFirst - 1 ) ) - C % nProperCells

  end subroutine SetChartFirstLast


end module ChartHeader_SL__Form
