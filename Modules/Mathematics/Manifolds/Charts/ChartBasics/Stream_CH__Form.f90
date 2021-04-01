module Stream_CH__Form

  !-- Stream_ChartHeader_Form

  use Basics
  use ManifoldBasics
  use Chart_H__Form
  use FieldSet_CH__Form

  implicit none
  private

  type, public :: Stream_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDL ) :: &
      Name = '', &
      Type = ''
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    class ( Stream_MH_Form ), pointer :: &
      Stream_M => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_SC
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


  subroutine Initialize ( SC, C, SM )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    class ( Stream_MH_Form ), intent ( in ), target :: &
      SM

    SC % IGNORABILITY  =  C % IGNORABILITY

    if ( SC % Type == '' ) &
      SC % Type = 'a Stream_C' 
    
    SC % Name  =  SM % Name

    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    SC % Chart            =>    C
    SC % Stream_M         =>   SM

  end subroutine Initialize


  subroutine Show_SC ( SC )

    class ( Stream_CH_Form ), intent ( in ) :: &
      SC

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( SC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SC % IGNORABILITY )

    associate &
      ( M  =>  SC % Chart % Manifold, &
        C  =>  SC % Chart )
    call Show (  SC % Name,    'Name',            SC % IGNORABILITY )
    call Show (   M % Name,    'Manifold',        SC % IGNORABILITY )   
    call Show (   C % Name,    'Chart',           SC % IGNORABILITY )   
    end associate !-- GIS

  end subroutine Show_SC


  impure elemental subroutine Finalize ( SC )

    type ( Stream_CH_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % Stream_M )
    nullify ( SC % Chart )

    if ( SC % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

  end subroutine Finalize


end module Stream_CH__Form
