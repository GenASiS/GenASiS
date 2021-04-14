module Stream_CH__Form

  !-- Stream_ChartHeader__Form

  use Basics
  use Manifolds
  use FieldSets

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      MAX_FIELD_SETS = 96

  type, public :: Stream_CH_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFieldSets   = 0
    logical ( KDL ) :: &
      Verbose = .false.
    character ( LDL ) :: &
      Type = '', &
      Name
    type ( GridImageStreamForm ), pointer :: &
      GridImageStream => null ( )
    class ( Chart_H_Form ), pointer :: &
      Chart => null ( )
    type ( FieldSet_C_Pointer ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Show => Show_SC
    final :: &
      Finalize
  end type Stream_CH_Form


contains


  subroutine Initialize_H ( SC, C, GIS, NameOption, VerboseOption )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( Chart_H_Form ), intent ( inout ), target :: &
      C
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    SC % IGNORABILITY  =  C % IGNORABILITY

    if ( SC % Type == '' ) &
      SC % Type = 'a Stream_C' 
    
    SC % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      SC % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

    SC % Verbose  =  .false.
    if ( present ( VerboseOption ) ) &
      SC % Verbose  =  VerboseOption
    
    SC % GridImageStream  =>  GIS
    SC % Chart            =>  C

    allocate ( SC % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize_H


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
                SC % IGNORABILITY  +  1 )
    call Show (  SC % Name, 'Stream',   SC % IGNORABILITY  +  1 )
    call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY  +  1 )

    end associate !-- nFS

  end subroutine AddFieldSet


  subroutine Show_SC ( SC )

    class ( Stream_CH_Form ), intent ( in ) :: &
      SC

    integer ( KDI ) :: &
      iFS
    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( SC % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SC % IGNORABILITY )

    associate &
      (   C  =>  SC % Chart, &
        GIS  =>  SC % GridImageStream )
    call Show (  SC % Name,       'Name',            SC % IGNORABILITY )
    call Show (   C % Name,       'Chart',           SC % IGNORABILITY )
    call Show ( GIS % Name,       'GridImageStream', SC % IGNORABILITY )
    call Show (  SC % Verbose,    'Verbose',         SC % IGNORABILITY )
    end associate !-- C, etc.

    call Show ( SC % nFieldSets, 'nFieldSets', SC % IGNORABILITY )
    do iFS  =  1, SC % nFieldSets
      associate ( FSC  =>  SC % FieldSet ( iFS ) % Pointer )
      call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY )
      end associate !-- FS
    end do !-- iFS

  end subroutine Show_SC


  impure elemental subroutine Finalize ( SC )

    type ( Stream_CH_Form ), intent ( inout ) :: &
      SC

    if ( allocated ( SC % FieldSet ) ) &
      deallocate ( SC % FieldSet )

    nullify ( SC % Chart )
    nullify ( SC % GridImageStream )

    call Show ( 'Finalizing ' // trim ( SC % Type ), SC % IGNORABILITY )
    call Show ( SC % Name, 'Name', SC % IGNORABILITY )

  end subroutine Finalize


end module Stream_CH__Form
