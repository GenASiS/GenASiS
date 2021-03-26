module Stream_CB__Form

  !-- Stream_ChartBase_Form

  use Basics
  use ManifoldBasics
  use ChartBasics
  use FieldSet_CB__Form

  implicit none
  private

  type, public, extends ( Stream_CH_Form ) :: Stream_CB_Form
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    final :: &
      Finalize
  end type Stream_CB_Form

    private :: &
      AddStorage


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

  end subroutine Initialize


  subroutine AddFieldSet ( SC, FSC )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      FSC

    call SC % Stream_CH_Form % AddFieldSet ( FSC )

    select type ( FSC )
    class is ( FieldSet_CB_Form )
      if ( SC % Verbose ) then
        call AddStorage ( SC, FSC % FieldSet )
      else
        call AddStorage ( SC, FSC % FieldSetStream )
      end if
    end select !-- FSC

  end subroutine AddFieldSet


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


end module Stream_CB__Form
