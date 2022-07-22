module Stream_Form

  use Basics
  use Manifolds
  use FieldSet_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3, &
      MAX_FIELD_SETS = 96

  type, public :: StreamForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFieldSets   = 0
    integer ( KDI ) :: &
      iTimerWrite = 0, &
      iTimerRead  = 0
    logical ( KDL ) :: &
      Verbose = .false.
    character ( LDL ) :: &
      Name
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    type ( CurveImageForm ), dimension ( : ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), dimension ( : ), allocatable :: &
      GridImage
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( FieldSetElement ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Show => Show_S
    procedure, public, pass :: &
      TimerWrite
    procedure, public, pass :: &
      TimerRead
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    final :: &
      Finalize
  end type StreamForm

   private :: &
     SetEdgeValues

    
contains


  subroutine Initialize ( S, A, GIS, NameOption, VerboseOption )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    integer ( KDI ) :: &
      iC  !-- iChart

    S % IGNORABILITY  =  A % IGNORABILITY

    S % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      S % Name  =  NameOption

    call Show ( 'Initializing a Stream', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    S % Verbose  =  .false.
    if ( present ( VerboseOption ) ) &
      S % Verbose  =  VerboseOption
    call PROGRAM_HEADER % GetParameter ( S % Verbose, 'VerboseStream' )
    
    S % GridImageStream  =>  GIS
    S % Atlas            =>  A

    associate &
      ( nC    =>  S % Atlas % nCharts, &
         GIS  =>  S % GridImageStream )

    allocate ( S % CurveImage ( nC ) )
    allocate ( S % GridImage  ( nC ) )

    do iC  =  1,  nC
      associate ( C  =>  S % Atlas % Chart ( iC ) % Element )
      select case ( C % nDimensions )
      case ( 1 ) 
        associate ( CI => S % CurveImage ( iC ) )
        call CI % Initialize ( GIS )
        end associate !-- CI
      case default
        associate ( GI => S % GridImage ( iC ) )
        call GI % Initialize ( GIS ) 
        end associate !-- GI
      end select !-- nDimensions
      end associate !-- C
    end do !-- iC

    allocate ( S % FieldSet ( MAX_FIELD_SETS ) )

    end associate !-- nC, etc.

  end subroutine Initialize


  subroutine AddFieldSet ( S, FS, NameOption, iaSelectedOption )

    class ( StreamForm ), intent ( inout ) :: &
      S
    class ( FieldSetForm ), intent ( in ) :: &
      FS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iC  !-- iChart
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( S % Verbose ) then
      allocate ( iaSelected, source = FS % iaSelected )
    else
      if ( present ( iaSelectedOption ) ) then
        allocate ( iaSelected, source = iaSelectedOption )
      else
        allocate ( iaSelected, source = FS % iaSelected )
      end if !-- iaSelectedOption
    end if !-- Verbose
      
    associate ( nFS  =>  S % nFieldSets )

    nFS  =  nFS + 1

    allocate ( S % FieldSet ( nFS ) % Element )
    associate ( FS_S  =>  S % FieldSet ( nFS ) % Element )

    call Show ( 'Adding a FieldSet to a Stream', &
                S % IGNORABILITY  +  1 )
    call Show (  S % Name, 'Stream',   S % IGNORABILITY  +  1 )
    call Show ( FS % Name, 'FieldSet', S % IGNORABILITY  +  1 )

    call FS_S % Initialize &
           ( FS, iaSelected, &
             NameOption = NameOption, &
             IgnorabilityOption = S % IGNORABILITY + 2 )
    
    do iC  =  1, S % Atlas % nCharts
      associate &
        (  C      =>  S % Atlas % Chart ( iC ) % Element, &
          FS_S_S  =>  FS_S % Storage ( iC ) )
      select case ( C % nDimensions )
      case ( 1 ) 
        associate ( CI => S % CurveImage ( iC ) )
        call CI % AddStorage ( FS_S_S )
        end associate !-- CI
      case default
        associate ( GI => S % GridImage ( iC ) )
        call GI % AddStorage ( FS_S_S )
        end associate !-- GI
      end select !-- nDimensions
      end associate !-- C
    end do !-- iC

    end associate !-- FS_S
    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Show_S ( S )

    class ( StreamForm ), intent ( in ) :: &
      S

    integer ( KDI ) :: &
      iFS

    call Show ( 'Stream Parameters', S % IGNORABILITY )

    associate &
      (   A  =>  S % Atlas, &
        GIS  =>  S % GridImageStream )
    call Show (   S % Name,    'Name',            S % IGNORABILITY )
    call Show (   A % Name,    'Atlas',           S % IGNORABILITY )
    call Show ( GIS % Name,    'GridImageStream', S % IGNORABILITY )
    call Show (   S % Verbose, 'Verbose',         S % IGNORABILITY )
    end associate !-- C, etc.

    call Show ( S % nFieldSets, 'nFieldSets', S % IGNORABILITY )
    do iFS  =  1, S % nFieldSets
      associate ( FS  =>  S % FieldSet ( iFS ) % Element )
      call Show ( FS % Name, 'FieldSet', S % IGNORABILITY )
      call FS % Show ( )
      end associate !-- FS
    end do !-- iFS

  end subroutine Show_S


  function TimerWrite ( S, Level ) result ( T )

    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimerWrite, &
               Name = trim ( S % Name ) // '_Wrt', &
               Level = Level )

  end function TimerWrite


  function TimerRead ( S, Level ) result ( T )

    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimerRead, &
               Name = trim ( S % Name ) // '_Rd', &
               Level = Level )

  end function TimerRead


  subroutine Write ( S, DirectoryOption, TimeOption, CycleNumberOption )

    class ( StreamForm ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( QuantityForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      nCellsProper, &
      nCellsGhost
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter, &
      nCellsWrite
    type ( Real_1D_Form ), dimension ( MAX_DIMENSIONS ) :: &
      Edge
    character ( LDF ) :: &
      Directory

    call Show ( 'Writing a Stream', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    associate ( A  =>  S % Atlas )

    do iC  =  1,  A % nCharts 
      select type ( C  =>  A % Chart ( iC ) % Element )
      class is ( Chart_GS_Form )

      if ( C % Distributed ) then
        nCellsProper    =  C % nCellsProper
        nCellsGhost     =  C % nCellsGhost
        nGhostInner     =  C % nGhostLayers
        nGhostOuter     =  C % nGhostLayers
        nExteriorInner  =  0
        nExteriorOuter  =  0
        where ( C % iaBrick  ==  1 )
          nGhostInner     =  0
          nExteriorInner  =  C % nGhostLayers
        end where
        where ( C % iaBrick  ==  C % nBricks )
          nGhostOuter     =  0
          nExteriorOuter  =  C % nGhostLayers
        end where
      else  ! .not. Distributed
        nCellsProper    =  C % nCellsProper
        nCellsGhost     =  0
        nGhostInner     =  0
        nGhostOuter     =  0
        nExteriorInner  =  C % nGhostLayers
        nExteriorOuter  =  C % nGhostLayers
      end if !-- Distributed
      nCellsWrite  =  C % nCellsBrick  +  nGhostInner  +  nGhostOuter

      call SetEdgeValues ( Edge, C )

      if ( A % nCharts  ==  1 ) then
        Directory  =  trim ( C % Name )  //  '/'
      else
        Directory  =  trim ( A % Name ) // '/' // trim ( C % Name )  //  '/'
      end if
      if ( present ( DirectoryOption ) ) &
        Directory  =  DirectoryOption

      select case ( C % nDimensions )
      case ( 1 ) 
        associate ( CI => S % CurveImage ( iC ) )
        call CI % SetGridWrite &
               ( Directory, Edge ( 1 ), nCellsProper, &
                 oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ), &
                 CoordinateLabelOption = C % CoordinateLabel ( 1 ), &
                 CoordinateUnitOption = C % CoordinateUnit ( 1 ) )
        call CI % Write &
               ( TimeOption = TimeOption, &
                 CycleNumberOption = CycleNumberOption )
        call CI % ClearGrid ( )
        end associate !-- CI
      case default
        associate ( GI => S % GridImage ( iC ) )
        call GI % SetGridWrite &
               ( Directory, Edge, nCellsWrite, nGhostInner, nGhostOuter, &
                 nExteriorInner, nExteriorOuter, C % nDimensions, &
                 nCellsProper, nCellsGhost, &
                 CoordinateLabelOption = C % CoordinateLabel, &
                 CoordinateUnitOption = C % CoordinateUnit )
        call GI % Write &
               ( TimeOption = TimeOption, &
                 CycleNumberOption = CycleNumberOption )
        call GI % ClearGrid ( )
        end associate !-- GI
      end select !-- nDimensions

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Stream_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Write', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- A

  end subroutine Write


  subroutine Read ( S, DirectoryOption, TimeOption, CycleNumberOption )

    class ( StreamForm ), intent ( inout ) :: &
      S
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( QuantityForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption

    integer ( KDI ) :: &
      iC, &  !-- iChart
      nCellsProper, &
      nCellsGhost
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter, &
      nCellsRead
    character ( LDF ) :: &
      Directory

    call Show ( 'Reading a Stream', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    associate ( A  =>  S % Atlas )

    do iC  =  1,  A % nCharts 
      select type ( C  =>  A % Chart ( iC ) % Element )
      class is ( Chart_GS_Form )

      if ( C % Distributed ) then
        nCellsProper    =  C % nCellsProper
        nCellsGhost     =  C % nCellsGhost
        nGhostInner     =  C % nGhostLayers
        nGhostOuter     =  C % nGhostLayers
        nExteriorInner  =  0
        nExteriorOuter  =  0
        where ( C % iaBrick  ==  1 )
          nGhostInner     =  0
          nExteriorInner  =  C % nGhostLayers
        end where
        where ( C % iaBrick  ==  C % nBricks )
          nGhostOuter     =  0
          nExteriorOuter  =  C % nGhostLayers
        end where
      else  ! .not. Distributed
        nCellsProper    =  C % nCellsProper
        nCellsGhost     =  0
        nGhostInner     =  0
        nGhostOuter     =  0
        nExteriorInner  =  C % nGhostLayers
        nExteriorOuter  =  C % nGhostLayers
      end if !-- Distributed
      nCellsRead  =  C % nCellsBrick  +  nGhostInner  +  nGhostOuter

      if ( A % nCharts  ==  1 ) then
        Directory  =  trim ( C % Name )  //  '/'
      else
        Directory  =  trim ( A % Name ) // '/' // trim ( C % Name )  //  '/'
      end if
      if ( present ( DirectoryOption ) ) &
        Directory  =  DirectoryOption

      select case ( C % nDimensions )
      case ( 1 ) 
        associate ( CI => S % CurveImage ( iC ) )
        call CI % SetGridRead &
               ( Directory, nCellsProper, &
                 oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ) )
        call CI % Read &
               ( StorageOnlyOption = .true., &
                 TimeOption = TimeOption, &
                 CycleNumberOption = CycleNumberOption )
        call CI % ClearGrid ( )
        end associate !-- CI
      case default
        associate ( GI => S % GridImage ( iC ) )
        call GI % SetGridRead &
               ( Directory, nCellsRead, nGhostInner, nGhostOuter, &
                 nExteriorInner, nExteriorOuter, C % nDimensions, &
                 nCellsProper, nCellsGhost )
        call GI % Read &
               ( StorageOnlyOption = .true., &
                 TimeOption = TimeOption, &
                 CycleNumberOption = CycleNumberOption )
        call GI % ClearGrid ( )
        end associate !-- GI
      end select !-- nDimensions

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Stream_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Read', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- A

  end subroutine Read


  impure elemental subroutine Finalize ( S )

    type ( StreamForm ), intent ( inout ) :: &
      S

    if ( allocated ( S % FieldSet ) ) &
      deallocate ( S % FieldSet )

    nullify ( S % Atlas )

    if ( allocated ( S % GridImage ) ) &
      deallocate ( S % GridImage )
    if ( allocated ( S % CurveImage ) ) &
      deallocate ( S % CurveImage )
    
    nullify ( S % GridImageStream )

    call Show ( 'Finalizing a Stream', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

  end subroutine Finalize


  subroutine SetEdgeValues ( Edge, C )

    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge
    class ( Chart_H_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      oE, &  !-- oEdge
      nE     !-- nEdges

    select type ( C )
    class is ( Chart_GS_Form )

    do iD = 1, C % nDimensions

      associate &
        ( nCB  =>  C % nCellsBrick ( iD ), &
          nGL  =>  C % nGhostLayers ( iD ), &
          iaB  =>  C % iaBrick ( iD ) )

      oE  =  ( iaB - 1 ) * nCB  -  nGL
      nE  =  nCB  +  2 * nGL  +  1

      call Edge ( iD ) % Initialize ( nE )
      Edge ( iD ) % Value  =  C % Edge ( iD ) % Value ( oE + 1 : oE + nE )

      end associate !-- nCB, etc.

    end do !-- iD

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Stream_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetEdgeValues', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine SetEdgeValues


end module Stream_Form
