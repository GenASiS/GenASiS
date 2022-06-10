module Atlas_H__Form

  !-- Atlas_Header__Form

  use Basics
  use Charts

  implicit none
  private

  type, public :: Atlas_H_Form
    integer ( KDI ) :: &
      IGNORABILITY, &
      nCharts
    character ( LDL ) :: &
      Type = '', &
      Name
    type ( ChartElement ), dimension ( : ), allocatable :: &
      Chart
  contains
    procedure, public, pass :: &
      Initialize_H
    procedure, private, pass :: &
      Show_A
    generic, public :: &
      Show => Show_A
    final :: &
      Finalize
  end type Atlas_H_Form


contains


  subroutine Initialize_H &
               ( A, NameOption, IgnorabilityOption, nChartsOption )

    class ( Atlas_H_Form ), intent ( inout ) :: &
      A
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
    integer ( KDI ), intent ( in ), optional :: &
      nChartsOption

    A % IGNORABILITY  =  CONSOLE % INFO_1
    if ( present ( IgnorabilityOption ) ) &
      A % IGNORABILITY  =  IgnorabilityOption

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas'

    A % Name  =  'Atlas'
    if ( present ( NameOption ) ) &
      A % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    A % nCharts  =  1
    if ( present ( nChartsOption ) ) &
      A % nCharts  =  nChartsOption
    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( A % nCharts ) )

  end subroutine Initialize_H


  subroutine Show_A ( A )

    class ( Atlas_H_Form ), intent ( in ) :: &
      A

   integer ( KDI ) :: &
     iC  !-- iC
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( A % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    call Show ( A % nCharts, 'nCharts', A % IGNORABILITY )
    do iC  =  1, A % nCharts
      if ( allocated ( A % Chart ( iC ) % Element ) ) then
        associate ( C  =>  A % Chart ( iC ) % Element )
        call C % Show ( )
        end associate !-- C
      end if  
    end do !-- iC

  end subroutine Show_A


  impure elemental subroutine Finalize ( A )

    type ( Atlas_H_Form ), intent ( inout ) :: &
      A

    if ( allocated ( A % Chart ) ) &
      deallocate ( A % Chart )  

    call Show ( 'Finalizing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

  end subroutine Finalize


end module Atlas_H__Form
