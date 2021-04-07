module Stream_CB__Form

  !-- Stream_ChartBase_Form

  use Basics
  use ManifoldBasics
  use ChartBasics
  use Chart_BH__Form
  use FieldSet_CB__Form

  implicit none
  private

  type, public, extends ( Stream_CH_Form ) :: Stream_CB_Form
    integer ( KDI ) :: &
      iTimerWrite = 0, &
      iTimerRead  = 0
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    final :: &
      Finalize
  end type Stream_CB_Form

    private :: &
      AddStorage, &
      SetEdgeValues


contains


  subroutine Initialize ( SC, C, SM )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    class ( Stream_MH_Form ), intent ( in ), target :: &
      SM

    if ( SC % Type == '' ) &
      SC % Type = 'a Stream_CB' 

    call SC % Stream_CH_Form % Initialize ( C, SM )

    associate ( GIS  =>  SC % Stream_M % GridImageStream )
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

  end subroutine Initialize


  subroutine AddFieldSet ( SC, FSC )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC

    select type ( FSC )
    class is ( FieldSet_CB_Form )
      if ( SC % Stream_M % Verbose ) then
        call AddStorage ( SC, FSC % FieldSet )
      else
        call AddStorage ( SC, FSC % FieldSetStream )
      end if
    end select !-- FSC

  end subroutine AddFieldSet


  subroutine Write ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                     TimerLevelOption )

    class ( Stream_CB_Form ), intent ( inout ) :: &
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
    integer ( KDI ), dimension ( MANIFOLD % MAX_DIMENSIONS ) :: &
      nGhostInner, &
      nGhostOuter, &
      nExteriorInner, &
      nExteriorOuter, &
      nCellsWrite
    type ( Real_1D_Form ), dimension ( MANIFOLD % MAX_DIMENSIONS ) :: &
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
    class is ( Chart_BH_Form )

    if ( C % Manifold % Distributed ) then
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

    end select !-- C

    call T % Stop ( )

  end subroutine Write


  subroutine Read ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                    TimerLevelOption )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    ! integer ( KDI ) :: &
    !   nCellsProper, &
    !   nCellsGhost
    ! integer ( KDI ), dimension ( ATLAS % MAX_DIMENSIONS ) :: &
    !   nCells, &
    !   nGhostInner, &
    !   nGhostOuter, &
    !   nExteriorInner, &
    !   nExteriorOuter
    ! character ( 2 ) :: &
    !   ChartNumber
    character ( LDF ) :: &
    !   Directory, &
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
    class is ( Chart_BH_Form )

    ! if ( C % IsDistributed ) then
    !   nCellsProper = C % nCellsProper
    !   nCellsGhost  = C % nCellsGhost
    !   nGhostInner = C % nGhostLayers
    !   nGhostOuter = C % nGhostLayers
    !   nExteriorInner = 0
    !   nExteriorOuter = 0
    !   where ( C % iaBrick == 1 )
    !     nGhostInner = 0
    !     nExteriorInner = C % nGhostLayers
    !   end where
    !   where ( C % iaBrick == C % nBricks )
    !     nGhostOuter = 0
    !     nExteriorOuter = C % nGhostLayers
    !   end where
    ! else  ! .not. IsDistributed
    !   nCellsProper = C % nCellsProper
    !   nCellsGhost  = 0
    !   nGhostInner = 0
    !   nGhostOuter = 0
    !   nExteriorInner = C % nGhostLayers
    !   nExteriorOuter = C % nGhostLayers
    ! end if !-- IsDistributed

    ! write ( ChartNumber, fmt = '(i2.2)' ) C % iChart
    ! Directory = 'Chart_' // ChartNumber // '/'
    ! if ( present ( DirectoryOption ) ) &
    !   Directory = DirectoryOption

    ! select case ( C % nDimensions )
    ! case ( 1 ) 
    !   associate ( CI => SC % CurveImage )
    !   call CI % SetGridRead &
    !          ( Directory, nCellsProper, &
    !            oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ) )
    !   call CI % Read &
    !          ( StorageOnlyOption = .true., &
    !            TimeOption = TimeOption, &
    !            CycleNumberOption = CycleNumberOption )
    !   call CI % ClearGrid ( )
    !   end associate !-- CI
    ! case default
    !   associate ( GI => SC % GridImage )
    !   call GI % SetGridRead &
    !          ( Directory, nCells, nGhostInner, nGhostOuter, &
    !            nExteriorInner, nExteriorOuter, C % nDimensions, nCellsProper, &
    !            nCellsGhost )
    !   call GI % Read &
    !          ( StorageOnlyOption = .true., &
    !            TimeOption = TimeOption, &
    !            CycleNumberOption = CycleNumberOption )
    !   call GI % ClearGrid ( )
    !   end associate !-- GI
    ! end select !-- nDimensions

    end select !-- C

    call T % Stop ( )

  end subroutine Read


  impure elemental subroutine Finalize ( SC )

    type ( Stream_CB_Form ), intent ( inout ) :: &
      SC

    if ( allocated ( SC % GridImage ) ) &
      deallocate ( SC % GridImage )
    if ( allocated ( SC % CurveImage ) ) &
      deallocate ( SC % CurveImage )

  end subroutine Finalize


  subroutine AddStorage ( SC, S )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( StorageForm ), intent ( in ) :: &
      S

    if ( allocated ( SC % CurveImage ) ) then
      call SC % CurveImage % AddStorage ( S )
    else if ( allocated ( SC % GridImage ) ) then
      call SC % GridImage % AddStorage ( S )
    end if

  end subroutine AddStorage


  subroutine SetEdgeValues ( SC, Edge )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      oE, &  !-- oEdge
      nE     !-- nEdges

    select type ( C  =>  SC % Chart )
    class is ( Chart_BH_Form )

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

    end select !-- C

  end subroutine SetEdgeValues


end module Stream_CB__Form
