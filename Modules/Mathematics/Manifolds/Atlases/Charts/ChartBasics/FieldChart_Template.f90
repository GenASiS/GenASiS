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
      NameShort = ''
    class ( ChartTemplate ), pointer :: &
      Chart => null ( )
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      FinalizeTemplate
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type FieldChartTemplate

    abstract interface 
      subroutine SF ( FC )
        import FieldChartTemplate
        class ( FieldChartTemplate ), intent ( inout ) :: &
          FC
      end subroutine
    end interface

  type, public :: FieldChartPointer
    class ( FieldChartTemplate ), pointer :: &
      Pointer => null ( )
  end type FieldChartPointer

contains


  subroutine InitializeTemplate ( FC, C, NameShort, IgnorabilityOption )

    class ( FieldChartTemplate ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ), optional :: &
      NameShort
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    FC % IGNORABILITY = C % IGNORABILITY
    if ( present ( IgnorabilityOption ) ) &
      FC % IGNORABILITY = IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a FieldChart' 

    FC % Name = trim ( NameShort ) // '_' // trim ( C % Name ) 

    call Show ( 'Initializing ' // trim ( FC % Type ), FC % IGNORABILITY )
    call Show ( FC % Name, 'Name', FC % IGNORABILITY )
   
    FC % NameShort = NameShort
    call Show ( FC % NameShort, 'NameShort', FC % IGNORABILITY )

    FC % Chart => C

    call FC % SetField ( )

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
