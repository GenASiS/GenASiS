module Stream_CH__Form

  use Basics
  use ManifoldBasics
  use Chart_H__Form
  use FieldSet_CH__Form

  implicit none
  private

  type, public :: Stream_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFieldSets   = 0
    character ( LDF ) :: &
      Name = '', &
      Type = ''
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    type ( FieldSet_CH_Pointer ), dimension ( : ), allocatable :: &
      FieldSet
    class ( Stream_MH_Form ), pointer :: &
      Stream_M => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    final :: &
      Finalize
  end type Stream_CH_Form

  type, public :: Stream_CH_Pointer
    class ( Stream_CH_Form ), pointer :: &
      Pointer => null ( )
  end type Stream_CH_Pointer

    integer ( KDI ), private, parameter :: &
      MAX_FIELD_SETS = MANIFOLD % MAX_FIELD_SETS


contains


  subroutine Initialize ( SC, SM, C, GIS, Name )

    class ( Stream_CH_Form ), intent ( inout ), target :: &
      SC
    class ( Stream_MH_Form ), intent ( in ), target :: &
      SM
    class ( Chart_H_Form ), intent ( in ), target :: &
      C
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ) :: &
      Name

    SC % IGNORABILITY  =  C % IGNORABILITY  +  1
    SC % Name          =  Name

    if ( SC % Type == '' ) &
      SC % Type = 'a Stream_C' 
    
    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    SC % GridImageStream  =>  GIS
    SC % Chart            =>    C
    SC % Stream_M         =>   SM

    allocate ( SC % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize


  subroutine AddFieldSet ( SC, FSC )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ), target :: &
      FSC
    
    integer ( KDI ) :: &
      iFS

    associate ( nFS  =>  SC % nFieldSets )

    do iFS  =  1, nFS
      if ( associated ( SC % FieldSet ( iFS ) % Pointer, FSC ) ) then
        call Show ( 'FieldSet already added to ' // SC % Type, &
                    CONSOLE % WARNING )
        call Show (  SC % Name, 'Stream',   CONSOLE % WARNING )
        call Show ( FSC % Name, 'FieldSet', CONSOLE % WARNING )
        return
      end if
    end do !-- iFS

    nFS  =  nFS + 1
    SC % FieldSet ( iFS ) % Pointer  =>  FSC
    call Show ( 'Adding a FieldSet to ' // trim ( SC % Type ), &
                SC % IGNORABILITY )
    call Show (  SC % Name, 'Stream',   SC % IGNORABILITY )
    call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY )

    end associate !-- nFS

  end subroutine AddFieldSet


  impure elemental subroutine Finalize ( SC )

    type ( Stream_CH_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % Stream_M )

    if ( allocated ( SC % FieldSet ) ) &
      deallocate ( SC % FieldSet )

    nullify ( SC % Chart )
    nullify ( SC % GridImageStream )

    if ( SC % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

  end subroutine Finalize


end module Stream_CH__Form
