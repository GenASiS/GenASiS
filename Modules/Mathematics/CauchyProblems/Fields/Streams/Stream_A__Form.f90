module Stream_A__Form

  !-- Stream_Atlas__Form

  use Basics
  use Manifolds
  use FieldSets
  use Stream_C__Form

  implicit none
  private

  type, public :: Stream_A_Form
    integer ( KDI ) :: &
      IGNORABILITY
    character ( LDL ) :: &
      Type = '', &
      Name
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    type ( Stream_C_Element ), dimension ( : ), allocatable :: &
      Stream_C
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      AddFieldSet
    procedure, public, pass :: &
      Write
    procedure, public, pass :: &
      Read
    procedure, public, pass :: &
      Show => Show_SA
    final :: &
      Finalize
  end type Stream_A_Form


contains


  subroutine Initialize ( SA, A, GIS, NameOption, VerboseOption )

    class ( Stream_A_Form ), intent ( inout ), target :: &
      SA
    class ( Atlas_H_Form ), intent ( inout ), target :: &
      A
    type ( GridImageStreamForm ), intent ( in ), target :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      VerboseOption

    integer ( KDI ) :: &
      iC  !-- iChart

    SA % IGNORABILITY  =  A % IGNORABILITY

    if ( SA % Type  ==  '' ) &
      SA % Type  =  'a Stream_A'

    SA % Name  =  'Stream'
    if ( present ( NameOption ) ) &
      SA % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( SA % Type ), A % IGNORABILITY )
    call Show ( SA % Name, 'Name', A % IGNORABILITY )

    SA % Atlas  =>  A

    associate ( nC  =>  A % nCharts )
    allocate ( SA % Stream_C ( nC ) )
    do iC  =  1, nC
      allocate ( SA % Stream_C ( iC ) % Element )
      associate ( SC  =>  SA % Stream_C ( iC ) % Element )
      associate (  C  =>   A %    Chart ( iC ) % Element )
      call SC % Initialize ( C, GIS, NameOption, VerboseOption )
      end associate !-- C
      end associate !-- SC
    end do !-- iC
    end associate !-- nC

  end subroutine Initialize


  subroutine AddFieldSet ( SA, FSA, NameOption, iaSelectedOption )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    class ( FieldSet_A_Form ), intent ( in ) :: &
      FSA
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, SA % Atlas % nCharts
      associate &
        (  SC  =>   SA %   Stream_C ( iC ) % Element, &
          FSC  =>  FSA % FieldSet_C ( iC ) % Element )
      call SC % AddFieldSet ( FSC, NameOption, iaSelectedOption )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine AddFieldSet


  subroutine Write ( SA, DirectoryOption, TimeOption, CycleNumberOption, &
                     TimerLevelOption )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, SA % Atlas % nCharts
      associate ( SC  =>  SA % Stream_C ( iC ) % Element )
      call SC % Write &
             ( DirectoryOption, TimeOption, CycleNumberOption, &
               TimerLevelOption )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine Write


  subroutine Read ( SA, DirectoryOption, TimeOption, CycleNumberOption, &
                    TimerLevelOption )

    class ( Stream_A_Form ), intent ( inout ) :: &
      SA
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption, &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, SA % Atlas % nCharts
      associate ( SC  =>  SA % Stream_C ( iC ) % Element )
      call SC % Write &
             ( DirectoryOption, TimeOption, CycleNumberOption, &
               TimerLevelOption )
      end associate !-- SC, etc.
    end do !-- iC

  end subroutine Read


  subroutine Show_SA ( SA )

    class ( Stream_A_Form ), intent ( in ) :: &
      SA

   integer ( KDI ) :: &
     iC  !-- iChart
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( SA % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', SA % IGNORABILITY )

    associate ( A  =>  SA % Atlas )

    call Show ( SA % Name, 'Name',  SA % IGNORABILITY )
    call Show (  A % Name, 'Atlas', SA % IGNORABILITY )

    do iC  =  1, A % nCharts
      if ( allocated ( SA % Stream_C ( iC ) % Element ) ) then
        associate ( SC  =>  SA % Stream_C ( iC ) % Element )
        call SC % Show ( )
        end associate !-- SC
      end if  
    end do !-- iC

    end associate  !-- A

  end subroutine Show_SA


  impure elemental subroutine Finalize ( SA )

    type ( Stream_A_Form ), intent ( inout ) :: &
      SA

    if ( allocated ( SA % Stream_C ) ) &
      deallocate ( SA % Stream_C )  

    nullify ( SA % Atlas )

    call Show ( 'Finalizing ' // trim ( SA % Type ), SA % IGNORABILITY )
    call Show ( SA % Name, 'Name', SA % IGNORABILITY )

  end subroutine Finalize


end module Stream_A__Form
