module Stream_MH__Form

  use Basics
  use MANIFOLD_Singleton
  use Manifold_H__Form
  use FieldSet_MH__Form

  implicit none
  private

  type, public :: Stream_MH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFieldSets   = 0
    character ( LDF ) :: &
      Name = '', &
      Type = ''
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    class ( Manifold_H_Form ), pointer :: &
      Manifold => null ( )
    type ( FieldSet_MH_Pointer ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    final :: &
      Finalize
  end type Stream_MH_Form

    integer ( KDI ), private, parameter :: &
      MAX_FIELD_SETS = MANIFOLD % MAX_FIELD_SETS


contains


  subroutine Initialize ( SM, M, GIS, Name )

    class ( Stream_MH_Form ), intent ( inout ) :: &
      SM
    class ( Manifold_H_Form ), intent ( in ), target :: &
      M
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ) :: &
      Name

    SM % IGNORABILITY  =  M % IGNORABILITY  +  1
    SM % Name          =  Name

    if ( SM % Type == '' ) &
      SM % Type = 'a Stream_M' 
    
    call Show ( 'Initializing ' // trim ( SM % Type ), SM % IGNORABILITY )
    call Show ( SM % Name, 'Name', SM % IGNORABILITY )

    SM % GridImageStream  =>  GIS
    SM % Manifold         =>    M 

    allocate ( SM % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize


  subroutine AddFieldSet ( SM, FSM )

    class ( Stream_MH_Form ), intent ( inout ) :: &
      SM
    class ( FieldSet_MH_Form ), intent ( in ), target :: &
      FSM
    
    integer ( KDI ) :: &
      iFS
    associate ( nFS  =>  SM % nFieldSets )

    do iFS  =  1, nFS
      if ( associated ( SM % FieldSet ( iFS ) % Pointer, FSM ) ) then
        call Show ( 'FieldSet already added to a ' // SM % Type, &
                    CONSOLE % WARNING )
        call Show (  SM % Name, 'Stream',   CONSOLE % WARNING )
        call Show ( FSM % Name, 'FieldSet', CONSOLE % WARNING )
        return
      end if
    end do !-- iFS

    nFS = nFS + 1
    SM % FieldSet ( iFS ) % Pointer  =>  FSM
    call Show ( 'Adding a FieldSet to a ' // trim ( SM % Type ), &
                SM % IGNORABILITY )
    call Show (  SM % Name, 'Stream',   SM % IGNORABILITY )
    call Show ( FSM % Name, 'FieldSet', SM % IGNORABILITY )

    end associate !-- nFS

  end subroutine AddFieldSet


  impure elemental subroutine Finalize ( SM )

    type ( Stream_MH_Form ), intent ( inout ) :: &
      SM

    if ( allocated ( SM % FieldSet ) ) &
      deallocate ( SM % FieldSet )

    nullify ( SM % Manifold )
    nullify ( SM % GridImageStream )

    if ( SM % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( SM % Type ), SM % IGNORABILITY )
    call Show ( SM % Name, 'Name', SM % IGNORABILITY )

  end subroutine Finalize


end module Stream_MH__Form
