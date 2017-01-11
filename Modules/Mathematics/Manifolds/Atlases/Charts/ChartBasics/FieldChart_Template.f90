!-- FieldChart is a template for a set of fields on a Chart.

module FieldChart_Template

  use Basics
  use Chart_Template

  implicit none
  private

  type, public, abstract :: FieldChartTemplate
    integer ( KDI ) :: &
      IGNORABILITY = 0
    character ( LDF ) :: &
      Name = '', &
      Type = '', &
      NameOutput = ''
    class ( ChartTemplate ), pointer :: &
      Chart => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
  end type FieldChartTemplate

contains


  subroutine InitializeTemplate ( FC, C, NameOutputOption )

    class ( FieldChartTemplate ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    FC % IGNORABILITY = C % IGNORABILITY

    if ( FC % Type == '' ) &
      FC % Type = 'a FieldChart' 

    call Split ( FC % Type, ' ', TypeWord )
    FC % Name = trim ( TypeWord ( 2 ) ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
    if ( present ( NameOutputOption ) ) then
      FC % NameOutput = NameOutputOption
      call Show ( FC % NameOutput, 'NameOutput', FC % IGNORABILITY )
    end if

    FC % Chart => C

  end subroutine InitializeTemplate


  impure elemental subroutine FinalizeTemplate ( FC )

    class ( FieldChartTemplate ), intent ( inout ) :: &
      FC

    nullify ( FC % Chart )

    if ( FC % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
  end subroutine FinalizeTemplate


end module FieldChart_Template
