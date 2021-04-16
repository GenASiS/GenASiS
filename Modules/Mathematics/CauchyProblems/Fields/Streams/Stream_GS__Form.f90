module Stream_GS__Form

  !-- Stream_GridStream__Form

  use Basics
  use Manifolds
  use FieldSets
  use Stream_CH__Form

  implicit none
  private

  type, public, extends ( Stream_CH_Form ) :: Stream_GS_Form
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
  end type Stream_GS_Form

    private :: &
      AddStorage, &
      SetEdgeValues


contains


  subroutine Initialize ( SG, G, GIS, NameOption, VerboseOption )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    class ( Grid_S_Form ), intent ( inout ), target :: &
      G
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    if ( SG % Type == '' ) &
      SG % Type = 'a Stream_GS' 

    call SG % Initialize_H ( G, GIS, NameOption, VerboseOption )

    associate ( GIS  =>  SG % GridImageStream )
    select case ( G % nDimensions )
    case ( 1 ) 
      allocate ( SG % CurveImage )
      associate ( CI => SG % CurveImage )
      call CI % Initialize ( GIS )
      end associate !-- CI
    case default
      allocate ( SG % GridImage )
      associate ( GI => SG % GridImage )
      call GI % Initialize ( GIS ) 
      end associate !-- GI
    end select !-- nDimensions
    end associate !-- GIS

  end subroutine Initialize


  subroutine AddFieldSet ( SC, FSC, NameOption, iaSelectedOption )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    call SC % Stream_CH_Form % AddFieldSet ( FSC, NameOption, iaSelectedOption )

    associate ( nFS  =>  SC % nFieldSets )

    select type ( FSC )
    class is ( FieldSet_GS_Form )
      if ( SC % Verbose ) then
        select type ( FSC )
        class is ( FieldSet_GS_Form )
        call AddStorage ( SC, FSC % FieldSet )
        end select !-- FSC
      else
        select type ( FSC_S  =>  SC % FieldSet ( nFS ) % Element )
        class is ( FieldSet_GS_Form )
        call AddStorage ( SC, FSC_S % FieldSet )
        end select !-- FSC_S
      end if
    end select !-- FSC

    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Write ( SG, DirectoryOption, TimeOption, CycleNumberOption, &
                     TimerLevelOption )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    ! integer ( KDI ) :: &
    !   nCellsProper, &
    !   nCellsGhost
    ! integer ( KDI ), dimension ( MANIFOLD % MAX_DIMENSIONS ) :: &
    !   nGhostInner, &
    !   nGhostOuter, &
    !   nExteriorInner, &
    !   nExteriorOuter, &
    !   nCellsWrite
    ! type ( Real_1D_Form ), dimension ( MANIFOLD % MAX_DIMENSIONS ) :: &
    !   Edge
    ! character ( LDF ) :: &
    !   Directory, &
    !   TimerName
    ! type ( TimerForm ), pointer :: &
    !   T 

    ! call Show ( 'Writing ' // trim ( SG % Type ), SG % IGNORABILITY )
    ! call Show ( SG % Name, 'Name', SG % IGNORABILITY )

    ! associate ( iT  =>  SG % iTimerWrite )
    ! if ( iT == 0 ) then
    !   TimerName  =  'Write ' // trim ( SG % Name )
    !   if ( present ( TimerLevelOption ) ) then
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
    !   else
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
    !   end if
    ! end if
    ! end associate !-- iT

    ! T  =>  PROGRAM_HEADER % TimerPointer ( SG % iTimerWrite )
    ! call T % Start ( )

    ! select type ( C  =>  SG % Chart )
    ! class is ( Chart_BH_Form )

    ! if ( C % Manifold % Distributed ) then
    !   nCellsProper    =  C % nCellsProper
    !   nCellsGhost     =  C % nCellsGhost
    !   nGhostInner     =  C % nGhostLayers
    !   nGhostOuter     =  C % nGhostLayers
    !   nExteriorInner  =  0
    !   nExteriorOuter  =  0
    !   where ( C % iaBrick  ==  1 )
    !     nGhostInner     =  0
    !     nExteriorInner  =  C % nGhostLayers
    !   end where
    !   where ( C % iaBrick  ==  C % nBricks )
    !     nGhostOuter     =  0
    !     nExteriorOuter  =  C % nGhostLayers
    !   end where
    ! else  ! .not. Distributed
    !   nCellsProper    =  C % nCellsProper
    !   nCellsGhost     =  0
    !   nGhostInner     =  0
    !   nGhostOuter     =  0
    !   nExteriorInner  =  C % nGhostLayers
    !   nExteriorOuter  =  C % nGhostLayers
    ! end if !-- Distributed
    ! nCellsWrite  =  C % nCellsBrick  +  nGhostInner  +  nGhostOuter

    ! call SetEdgeValues ( SG, Edge )

    ! Directory  =  trim ( C % Name )  //  '/'
    ! if ( present ( DirectoryOption ) ) &
    !   Directory  =  DirectoryOption

    ! select case ( C % nDimensions )
    ! case ( 1 ) 
    !   associate ( CI => SG % CurveImage )
    !   call CI % SetGridWrite &
    !          ( Directory, Edge ( 1 ), nCellsProper, &
    !            oValue = nGhostInner ( 1 ) + nExteriorInner ( 1 ), &
    !            CoordinateLabelOption = C % CoordinateLabel ( 1 ), &
    !            CoordinateUnitOption = C % CoordinateUnit ( 1 ) )
    !   call CI % Write &
    !          ( TimeOption = TimeOption, &
    !            CycleNumberOption = CycleNumberOption )
    !   call CI % ClearGrid ( )
    !   end associate !-- CI
    ! case default
    !   associate ( GI => SG % GridImage )
    !   call GI % SetGridWrite &
    !          ( Directory, Edge, nCellsWrite, nGhostInner, nGhostOuter, &
    !            nExteriorInner, nExteriorOuter, C % nDimensions, nCellsProper, &
    !            nCellsGhost, CoordinateLabelOption = C % CoordinateLabel, &
    !            CoordinateUnitOption = C % CoordinateUnit )
    !   call GI % Write &
    !          ( TimeOption = TimeOption, &
    !            CycleNumberOption = CycleNumberOption )
    !   call GI % ClearGrid ( )
    !   end associate !-- GI
    ! end select !-- nDimensions

    ! end select !-- C

    ! call T % Stop ( )

  end subroutine Write


  subroutine Read ( SG, DirectoryOption, TimeOption, CycleNumberOption, &
                    TimerLevelOption )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
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
    ! integer ( KDI ), dimension ( MANIFOLD % MAX_DIMENSIONS ) :: &
    !   nGhostInner, &
    !   nGhostOuter, &
    !   nExteriorInner, &
    !   nExteriorOuter, &
    !   nCellsRead
    ! character ( LDF ) :: &
    !   Directory, &
    !   TimerName
    ! type ( TimerForm ), pointer :: &
    !   T 

    ! call Show ( 'Reading ' // trim ( SG % Type ), SG % IGNORABILITY )
    ! call Show ( SG % Name, 'Name', SG % IGNORABILITY )

    ! associate ( iT  =>  SG % iTimerRead )
    ! if ( iT == 0 ) then
    !   TimerName  =  'Read ' // trim ( SG % Name )
    !   if ( present ( TimerLevelOption ) ) then
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
    !   else
    !     call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
    !   end if
    ! end if
    ! end associate !-- iT

    ! T  =>  PROGRAM_HEADER % TimerPointer ( SG % iTimerRead )
    ! call T % Start ( )

    ! select type ( C  =>  SG % Chart )
    ! class is ( Chart_BH_Form )

    ! if ( C % Manifold % Distributed ) then
    !   nCellsProper    =  C % nCellsProper
    !   nCellsGhost     =  C % nCellsGhost
    !   nGhostInner     =  C % nGhostLayers
    !   nGhostOuter     =  C % nGhostLayers
    !   nExteriorInner  =  0
    !   nExteriorOuter  =  0
    !   where ( C % iaBrick  ==  1 )
    !     nGhostInner     =  0
    !     nExteriorInner  =  C % nGhostLayers
    !   end where
    !   where ( C % iaBrick  ==  C % nBricks )
    !     nGhostOuter     =  0
    !     nExteriorOuter  =  C % nGhostLayers
    !   end where
    ! else  ! .not. Distributed
    !   nCellsProper    =  C % nCellsProper
    !   nCellsGhost     =  0
    !   nGhostInner     =  0
    !   nGhostOuter     =  0
    !   nExteriorInner  =  C % nGhostLayers
    !   nExteriorOuter  =  C % nGhostLayers
    ! end if !-- Distributed
    ! nCellsRead  =  C % nCellsBrick  +  nGhostInner  +  nGhostOuter

    ! Directory  =  trim ( C % Name )  //  '/'
    ! if ( present ( DirectoryOption ) ) &
    !   Directory  =  DirectoryOption

    ! select case ( C % nDimensions )
    ! case ( 1 ) 
    !   associate ( CI => SG % CurveImage )
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
    !   associate ( GI => SG % GridImage )
    !   call GI % SetGridRead &
    !          ( Directory, nCellsRead, nGhostInner, nGhostOuter, &
    !            nExteriorInner, nExteriorOuter, C % nDimensions, nCellsProper, &
    !            nCellsGhost )
    !   call GI % Read &
    !          ( StorageOnlyOption = .true., &
    !            TimeOption = TimeOption, &
    !            CycleNumberOption = CycleNumberOption )
    !   call GI % ClearGrid ( )
    !   end associate !-- GI
    ! end select !-- nDimensions

    ! end select !-- C

    ! call T % Stop ( )

  end subroutine Read


  impure elemental subroutine Finalize ( SG )

    type ( Stream_GS_Form ), intent ( inout ) :: &
      SG

    if ( allocated ( SG % GridImage ) ) &
      deallocate ( SG % GridImage )
    if ( allocated ( SG % CurveImage ) ) &
      deallocate ( SG % CurveImage )

  end subroutine Finalize


  subroutine AddStorage ( SG, S )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    class ( StorageForm ), intent ( in ) :: &
      S

    if ( allocated ( SG % CurveImage ) ) then
      call SG % CurveImage % AddStorage ( S )
    else if ( allocated ( SG % GridImage ) ) then
      call SG % GridImage % AddStorage ( S )
    end if

  end subroutine AddStorage


  subroutine SetEdgeValues ( SG, Edge )

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SG
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      oE, &  !-- oEdge
      nE     !-- nEdges

    select type ( G  =>  SG % Chart )
    class is ( Grid_S_Form )

    do iD = 1, G % nDimensions

      associate &
        ( nCB  =>  G % nCellsBrick ( iD ), &
          nGL  =>  G % nGhostLayers ( iD ), &
          iaB  =>  G % iaBrick ( iD ) )

      oE  =  ( iaB - 1 ) * nCB  -  nGL
      nE  =  nCB  +  2 * nGL  +  1

      call Edge ( iD ) % Initialize ( nE )
      Edge ( iD ) % Value  =  G % Edge ( iD ) % Value ( oE + 1 : oE + nE )

      end associate !-- nCB, etc.

    end do !-- iD

    end select !-- C

  end subroutine SetEdgeValues


end module Stream_GS__Form
