!-- Field_CSL is a template for a set of fields on a Chart_SL.

module Field_CSL__Template

  !-- Field_ChartSingleLevel_Template

  use Basics
  use Chart_Template
  use FieldChart_Template

  implicit none
  private

  type, public, extends ( FieldChartTemplate ), abstract :: Field_CSL_Template
    integer ( KDI ) :: &
      nValues
    class ( VariableGroupForm ), allocatable :: &
      Field, &
      FieldOutput
  contains
    procedure, public, pass :: &
      InitializeTemplate_CSL
    procedure, public, pass :: &
      FinalizeTemplate_CSL
    procedure ( SF ), private, pass, deferred :: &
      SetField
  end type Field_CSL_Template

    abstract interface 
      subroutine SF ( FC, NameOption )
        import Field_CSL_Template
        class ( Field_CSL_Template ), intent ( inout ) :: &
          FC
        character ( * ), intent ( in ), optional :: &
          NameOption
      end subroutine
    end interface

  type, public :: Field_CSL_Pointer
    class ( Field_CSL_Template ), pointer :: &
      Pointer => null ( )
  end type Field_CSL_Pointer

contains


  subroutine InitializeTemplate_CSL ( FC, C, nValues, NameOutputOption )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    call FC % InitializeTemplate ( C, NameOutputOption = NameOutputOption )

    FC % nValues = nValues

    if ( FC % NameOutput == '' ) then
      call FC % SetField ( )
    else
      call FC % SetField ( NameOption = FC % NameOutput )
    end if

  end subroutine InitializeTemplate_CSL


  impure elemental subroutine FinalizeTemplate_CSL ( FC )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % FieldOutput ) ) &
      deallocate ( FC % FieldOutput )
    if ( allocated ( FC % Field ) ) &
      deallocate ( FC % Field )

    call FC % FinalizeTemplate ( )

  end subroutine FinalizeTemplate_CSL


end module Field_CSL__Template
