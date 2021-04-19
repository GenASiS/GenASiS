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
    type ( FieldSet_C_Element ), dimension ( : ), allocatable :: &
      FieldSet
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, public, pass :: &
      AddFieldSet_H
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Show => Show_SC
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    final :: &
      Finalize
    procedure, private, nopass :: &
      AllocateFieldSetElement
  end type Stream_CH_Form

  type, public :: Stream_C_Element
    !-- Stream_Chart_Element
    class ( Stream_CH_Form ), allocatable :: &
      Element
  contains
    final :: &
      Finalize_E
  end type Stream_C_Element


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
    call PROGRAM_HEADER % GetParameter ( SC % Verbose, 'VerboseStream' )
    
    SC % GridImageStream  =>  GIS
    SC % Chart            =>  C

    allocate ( SC % FieldSet ( MAX_FIELD_SETS ) )

  end subroutine Initialize_H


  subroutine AddFieldSet_H ( SC, FSC, NameOption, iaSelectedOption )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ), dimension ( : ), allocatable :: &
      iaSelected

    if ( present ( iaSelectedOption ) ) then
      allocate ( iaSelected, source = iaSelectedOption )
    else
      allocate ( iaSelected, source = FSC % iaSelected )
    end if

    associate ( nFS  =>  SC % nFieldSets )

    nFS  =  nFS + 1

    call SC % AllocateFieldSetElement ( SC % FieldSet ( nFS ) % Element )
    associate ( FSC_SC  =>  SC % FieldSet ( nFS ) % Element )
    call Show ( 'Adding a FieldSet to ' // trim ( SC % Type ), &
                SC % IGNORABILITY  +  1 )
    call Show (  SC % Name, 'Stream',   SC % IGNORABILITY  +  1 )
    call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY  +  1 )
    call FSC_SC % Clone &
           ( FSC, &
             NameOption = NameOption, &
             iaSelectedOption = iaSelected, &
             IgnorabilityOption = SC % IGNORABILITY + 1 )
    end associate !-- FSC_SC

    end associate !-- nFS

  end subroutine AddFieldSet_H


  subroutine AddFieldSet ( SC, FSC, NameOption, iaSelectedOption )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    class ( FieldSet_CH_Form ), intent ( in ) :: &
      FSC
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption
    
    call SC % AddFieldSet_H ( FSC, NameOption, iaSelectedOption )

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
      associate ( FSC  =>  SC % FieldSet ( iFS ) % Element )
      call Show ( FSC % Name, 'FieldSet', SC % IGNORABILITY )
      end associate !-- FS
    end do !-- iFS

  end subroutine Show_SC


  subroutine Write ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                     TimerLevelOption )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

  end subroutine Write


  subroutine Read ( SC, DirectoryOption, TimeOption, CycleNumberOption, &
                    TimerLevelOption )

    class ( Stream_CH_Form ), intent ( inout ) :: &
      SC
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

  end subroutine Read


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


  impure elemental subroutine Finalize_E ( SE )
    
    type ( Stream_C_Element ), intent ( inout ) :: &
      SE

    if ( allocated ( SE % Element ) ) &
      deallocate ( SE % Element )

  end subroutine Finalize_E


  subroutine AllocateFieldSetElement ( FSC )

    class ( FieldSet_CH_Form ), intent ( out ), allocatable :: &
      FSC

    allocate ( FSC )

  end subroutine AllocateFieldSetElement


end module Stream_CH__Form
