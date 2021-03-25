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
    class ( Manifold_H_Form ), pointer :: &
      Manifold => null ( )
    type ( FieldSet_MH_Pointer ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Show => Show_SM
    final :: &
      Finalize
  end type Stream_MH_Form

    integer ( KDI ), private, parameter :: &
      MAX_FIELD_SETS = MANIFOLD % MAX_FIELD_SETS


contains


  subroutine Initialize ( SM, M, GIS, NameShort, VerboseOption )

    class ( Stream_MH_Form ), intent ( inout ) :: &
      SM
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    logical ( KDL ) :: &
      Verbose

    Verbose = .false.
    if ( present ( VerboseOption ) ) &
      Verbose = VerboseOption
    
    SM % IGNORABILITY  =  M % IGNORABILITY

    if ( SM % Type == '' ) &
      SM % Type = 'a Stream_M' 
    
    SM % Name  =  trim ( NameShort ) // '_' // trim ( M % Name )

    call Show ( 'Initializing ' // trim ( SM % Type ), SM % IGNORABILITY )
    call Show ( SM % Name, 'Name', SM % IGNORABILITY )

     M % nStreams  =  M % nStreams  +  1
    SM % iStream   =  M % nStreams

    SM % NameShort  =  NameShort
    SM % Verbose    =  Verbose

    SM % GridImageStream  =>  GIS
    SM % Manifold         =>    M 

    allocate ( SM % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize


  subroutine AddFieldSet ( SM, FSM )

    class ( Stream_MH_Form ), intent ( inout ) :: &
      SM
    class ( FieldSet_MH_Form ), intent ( inout ), target :: &
      FSM
    
    integer ( KDI ) :: &
      iFS

    associate ( nFS  =>  SM % nFieldSets )

    do iFS  =  1, nFS
      if ( associated ( SM % FieldSet ( iFS ) % Pointer, FSM ) ) then
        call Show ( 'FieldSet already added to ' // SM % Type, &
                    CONSOLE % WARNING )
        call Show (  SM % Name, 'Stream',   CONSOLE % WARNING )
        call Show ( FSM % Name, 'FieldSet', CONSOLE % WARNING )
        return
      end if
    end do !-- iFS

    nFS  =  nFS + 1
    SM % FieldSet ( iFS ) % Pointer  =>  FSM
    call Show ( 'Adding a FieldSet to ' // trim ( SM % Type ), &
                SM % IGNORABILITY  +  1 )
    call Show (  SM % Name, 'Stream',   SM % IGNORABILITY  +  1 )
    call Show ( FSM % Name, 'FieldSet', SM % IGNORABILITY  +  1 )

    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Show_SM ( SM )

    class ( Stream_MH_Form ), intent ( in ) :: &
      SM

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( SM % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SM % IGNORABILITY )

    associate ( GIS  =>  SM % GridImageStream )
    call Show (  SM % Name,             'Name', SM % IGNORABILITY )
    call Show ( GIS % Name,  'GridImageStream', SM % IGNORABILITY )
    call Show (  SM % iStream,       'iStream', SM % IGNORABILITY )
    call Show (  SM % Verbose,       'Verbose', SM % IGNORABILITY )
    end associate !-- GIS

  end subroutine Show_SM


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
