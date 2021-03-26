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
      IGNORABILITY = 0, &
      iStream      = 0, &
      nFieldSets   = 0
    logical ( KDL ) :: &
      Verbose = .false.
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameShort = ''
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
    
    SC % Name  =  trim ( SM % NameShort ) // '_' // trim ( C % Name )

    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

     C % nStreams  =   C % Manifold % nStreams
    SC % iStream   =  SM % iStream

    SC % NameShort  =  SM % NameShort
    SC % Verbose    =  SM % Verbose

    SC % GridImageStream  =>  SM % GridImageStream
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

    nFS  =  SC % Stream_M % nFieldSets
    SC % FieldSet ( iFS ) % Pointer  =>  FSC
    call Show ( 'Adding a FieldSet to ' // trim ( SC % Type ), &
                SC % IGNORABILITY  +  1 )
    call Show (  SC % Name, 'Stream',   SC % IGNORABILITY  +  1 )
    call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY  +  1 )

    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Show_SC ( SC )

    class ( Stream_CH_Form ), intent ( in ) :: &
      SC

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( SC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SC % IGNORABILITY )

    associate ( GIS  =>  SC % GridImageStream )
    call Show (  SC % Name,             'Name', SC % IGNORABILITY )
    call Show ( GIS % Name,  'GridImageStream', SC % IGNORABILITY )
    call Show (  SC % iStream,       'iStream', SC % IGNORABILITY )
    call Show (  SC % Verbose,       'Verbose', SC % IGNORABILITY )
    end associate !-- GIS

  end subroutine Show_SC


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
