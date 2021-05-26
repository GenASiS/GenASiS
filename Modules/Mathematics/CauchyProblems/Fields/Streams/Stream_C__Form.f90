module Stream_C__Form

  !-- Stream_Chart__Form

  use Basics
  use Manifolds
  use FieldSets

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = 3, &
      MAX_FIELD_SETS = 96

  type, public :: Stream_C_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iTimerWrite  = 0, &
      iTimerRead   = 0, &
      nFieldSets   = 0
    logical ( KDL ) :: &
      Verbose = .false.
    character ( LDL ) :: &
      Type = '', &
      Name
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    type ( FieldSet_C_Element ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    procedure, public, pass :: &
      Show => Show_SC
    final :: &
      Finalize
  end type Stream_C_Form

  type, public :: Stream_C_Element
    !-- Stream_Chart_Element
    class ( Stream_C_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Stream_C_Element

    private :: &
      AddStorage, &
      SetEdgeValues

    
contains


  subroutine Initialize ( SC, C, GIS, NameOption, VerboseOption )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    SC % IGNORABILITY  =  C % IGNORABILITY

    if ( SC % Type == '' ) &
      SC % Type = 'a Stream_C' 
    
    SC % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      SC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    SC % Verbose  =  .false.
    if ( present ( VerboseOption ) ) &
      SC % Verbose  =  VerboseOption
    call PROGRAM_HEADER % GetParameter ( SC % Verbose, 'VerboseStream' )
    
    SC % GridImageStream  =>  GIS
    SC % Chart            =>  C

    associate ( GIS  =>  SC % GridImageStream )
    select case ( C % nDimensions )
    case ( 1 ) 
      allocate ( SC % CurveImage )
      associate ( CI => SC % CurveImage )
      call CI % Initialize ( GIS )
      end associate !-- CI
    case default
      allocate ( SC % GridImage )
      associate ( GI => SC % GridImage )
      call GI % Initialize ( GIS ) 
      end associate !-- GI
    end select !-- nDimensions
    end associate !-- GIS

    allocate ( SC % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize


  subroutine AddFieldSet ( SC, FSC, NameOption, iaSelectedOption )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_C_Form ), intent ( in ) :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( SC % Verbose ) then
      allocate ( iaSelected, source = FSC % iaSelected )
    else
      if ( present ( iaSelectedOption ) ) then
        allocate ( iaSelected, source = iaSelectedOption )
      else
        allocate ( iaSelected, source = FSC % iaSelected )
      end if !-- iaSelectedOption
    end if !-- Verbose
      
    associate ( nFS  =>  SC % nFieldSets )

    nFS  =  nFS + 1

    allocate ( SC % FieldSet ( nFS ) % Element )
    associate ( FSC_SC  =>  SC % FieldSet ( nFS ) % Element )

    call Show ( 'Adding a FieldSet to ' // trim ( SC % Type ), &
                SC % IGNORABILITY  +  1 )
    call Show (  SC % Name, 'Stream',   SC % IGNORABILITY  +  1 )
    call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY  +  1 )

    call FSC_SC % Initialize &
           ( FSC, iaSelected, &
             NameOption = NameOption, &
             IgnorabilityOption = SC % IGNORABILITY + 1 )
    
    call AddStorage ( SC, FSC_SC % Storage_FSC )

    end associate !-- FSC_SC
    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Write ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                     TimerLevelOption )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    integer ( KDI ) :: &
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
      Directory, &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    call Show ( 'Writing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    associate ( iT  =>  SC % iTimerWrite )
    if ( iT == 0 ) then
      TimerName  =  'Write ' // trim ( SC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( SC % iTimerWrite )
    call T % Start ( )

    select type ( C  =>  SC % Chart )
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

    call SetEdgeValues ( SC, Edge )

    Directory  =  trim ( C % Name )  //  '/'
    if ( present ( DirectoryOption ) ) &
      Directory  =  DirectoryOption

    select case ( C % nDimensions )
    case ( 1 ) 
      associate ( CI => SC % CurveImage )
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
      associate ( GI => SC % GridImage )
      call GI % SetGridWrite &
             ( Directory, Edge, nCellsWrite, nGhostInner, nGhostOuter, &
               nExteriorInner, nExteriorOuter, C % nDimensions, nCellsProper, &
               nCellsGhost, CoordinateLabelOption = C % CoordinateLabel, &
               CoordinateUnitOption = C % CoordinateUnit )
      call GI % Write &
             ( TimeOption = TimeOption, &
               CycleNumberOption = CycleNumberOption )
      call GI % ClearGrid ( )
      end associate !-- GI
    end select !-- nDimensions

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Stream_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Read', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    call T % Stop ( )

  end subroutine Write


  subroutine Read ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                    TimerLevelOption )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    integer ( KDI ) :: &
      nCellsProper, &
      nCellsGhost
    integer ( KDI ), dimension ( MAX_DIMENSIONS ) :: &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter, &
      nCellsRead
    character ( LDF ) :: &
      Directory, &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    call Show ( 'Reading ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    associate ( iT  =>  SC % iTimerRead )
    if ( iT == 0 ) then
      TimerName  =  'Read ' // trim ( SC % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( SC % iTimerRead )
    call T % Start ( )

    select type ( C  =>  SC % Chart )
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

    Directory  =  trim ( C % Name )  //  '/'
    if ( present ( DirectoryOption ) ) &
      Directory  =  DirectoryOption

    select case ( C % nDimensions )
    case ( 1 ) 
      associate ( CI => SC % CurveImage )
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
      associate ( GI => SC % GridImage )
      call GI % SetGridRead &
             ( Directory, nCellsRead, nGhostInner, nGhostOuter, &
               nExteriorInner, nExteriorOuter, C % nDimensions, nCellsProper, &
               nCellsGhost )
      call GI % Read &
             ( StorageOnlyOption = .true., &
               TimeOption = TimeOption, &
               CycleNumberOption = CycleNumberOption )
      call GI % ClearGrid ( )
      end associate !-- GI
    end select !-- nDimensions

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Stream_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Read', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

    call T % Stop ( )

  end subroutine Read


  subroutine Show_SC ( SC )

    class ( Stream_C_Form ), intent ( in ) :: &
      SC

    integer ( KDI ) :: &
      iFS
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( SC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SC % IGNORABILITY )

    associate &
      (   C  =>  SC % Chart, &
        GIS  =>  SC % GridImageStream )
    call Show (  SC % Name,       'Name',            SC % IGNORABILITY )
    call Show (   C % Name,       'Chart',           SC % IGNORABILITY )
    call Show ( GIS % Name,       'GridImageStream', SC % IGNORABILITY )
    call Show (  SC % Verbose,    'Verbose',         SC % IGNORABILITY )
    end associate !-- C, etc.

    call Show ( SC % nFieldSets, 'nFieldSets', SC % IGNORABILITY )
    do iFS  =  1, SC % nFieldSets
      associate ( FSC  =>  SC % FieldSet ( iFS ) % Element )
      call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY )
      call FSC % Show ( )
      end associate !-- FS
    end do !-- iFS

  end subroutine Show_SC


  impure elemental subroutine Finalize ( SC )

    type ( Stream_C_Form ), intent ( inout ) :: &
      SC

    if ( allocated ( SC % FieldSet ) ) &
      deallocate ( SC % FieldSet )

    nullify ( SC % Chart )

    if ( allocated ( SC % GridImage ) ) &
      deallocate ( SC % GridImage )
    if ( allocated ( SC % CurveImage ) ) &
      deallocate ( SC % CurveImage )
    
    nullify ( SC % GridImageStream )

    call Show ( 'Finalizing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

  end subroutine Finalize


  impure elemental subroutine Finalize_E ( SE )
    
    type ( Stream_C_Element ), intent ( inout ) :: &
      SE

    if ( allocated ( SE % Element ) ) &
      deallocate ( SE % Element )

  end subroutine Finalize_E


  subroutine AddStorage ( SC, S_FSC )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    class ( Storage_FSC_Form ), intent ( in ) :: &
      S_FSC

    if ( allocated ( SC % CurveImage ) ) then
      call SC % CurveImage % AddStorage ( S_FSC % Storage )
    else if ( allocated ( SC % GridImage ) ) then
      call SC % GridImage % AddStorage ( S_FSC % Storage )
    end if

  end subroutine AddStorage


  subroutine SetEdgeValues ( SC, Edge )

    class ( Stream_C_Form ), intent ( inout ) :: &
      SC
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      oE, &  !-- oEdge
      nE     !-- nEdges

    select type ( C  =>  SC % Chart )
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
      call Show ( 'Stream_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetEdgeValues', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine SetEdgeValues


end module Stream_C__Form
