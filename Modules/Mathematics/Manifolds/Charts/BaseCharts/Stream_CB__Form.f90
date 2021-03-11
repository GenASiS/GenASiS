module Stream_CB__Form

  !-- Stream_ChartBase_Form

  use Basics
  use Chart_BH__Form

  implicit none
  private

  type, public :: Stream_CB_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = ''
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    type ( CurveImageForm ), allocatable :: &
      CurveImage
    type ( StructuredGridImageForm ), allocatable :: &
      GridImage
    class ( Chart_BH_Form ), pointer :: &
      Chart => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    final :: &
      Finalize
  end type Stream_CB_Form


contains


  subroutine Initialize ( SC, C, GIS, Name )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( Chart_BH_Form ), intent ( in ), target :: &
      C
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ) :: &
      Name

    SC % IGNORABILITY = CONSOLE % INFO_4
    SC % Name = Name

    call Show ( 'Initializing a Stream_CB', SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    SC % GridImageStream  =>  GIS
    SC % Chart            =>    C 


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

  end subroutine Initialize


  subroutine AddFieldSet ( SC, S )

    class ( Stream_CB_Form ), intent ( inout ) :: &
      SC
    class ( StorageForm ), intent ( in ) :: &
      S

    if ( allocated ( SC % CurveImage ) ) then
      call SC % CurveImage % AddStorage ( S )
    else if ( allocated ( SC % GridImage ) ) then
      call SC % GridImage % AddStorage ( S )
    end if

  end subroutine AddFieldSet


  impure elemental subroutine Finalize ( SC )

    type ( Stream_CB_Form ), intent ( inout ) :: &
      SC

    if ( allocated ( SC % GridImage ) ) &
      deallocate ( SC % GridImage )
    if ( allocated ( SC % CurveImage ) ) &
      deallocate ( SC % CurveImage )

    nullify ( SC % GridImageStream )
    nullify ( SC % Chart )

    if ( SC % Name == '' ) return

    call Show ( 'Finalizing a Stream_CB', SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

  end subroutine Finalize


end module Stream_CB__Form
