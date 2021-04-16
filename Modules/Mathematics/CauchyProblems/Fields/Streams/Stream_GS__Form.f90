module Stream_GS__Form

  !-- Stream_GridStream__Form

  use Basics
  use Manifolds
  use FieldSets
  use Stream_CH__Form

  implicit none
  private

  type, public, extends ( Stream_CH_Form ) :: Stream_GS_Form
  !   integer ( KDI ) :: &
  !     iTimerWrite = 0, &
  !     iTimerRead  = 0
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
  !   procedure, public, pass :: &
  !     Write
  !   procedure, public, pass :: &
  !     Read
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

    ! select type ( FSC )
    ! class is ( FieldSet_GS_Form )
    !   if ( SC % Verbose ) then
    !     call AddStorage ( SC, FSC % FieldSet )
    !   else
    !     call AddStorage ( SC, FSC % FieldSetStream )
    !   end if
    ! end select !-- FSC

  end subroutine AddFieldSet


  impure elemental subroutine Finalize ( SG )

    type ( Stream_GS_Form ), intent ( inout ) :: &
      SG

    if ( allocated ( SG % GridImage ) ) &
      deallocate ( SG % GridImage )
    if ( allocated ( SG % CurveImage ) ) &
      deallocate ( SG % CurveImage )

  end subroutine Finalize


  subroutine AddStorage ( SC, S )

    class ( Stream_GS_Form ), intent ( inout ) :: &
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

    class ( Stream_GS_Form ), intent ( inout ) :: &
      SC
    type ( Real_1D_Form ), dimension ( : ), intent ( inout ) :: &
      Edge

    integer ( KDI ) :: &
      iD, &  !-- iDimension
      oE, &  !-- oEdge
      nE     !-- nEdges

    select type ( C  =>  SC % Chart )
    class is ( Grid_S_Form )

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


end module Stream_GS__Form
