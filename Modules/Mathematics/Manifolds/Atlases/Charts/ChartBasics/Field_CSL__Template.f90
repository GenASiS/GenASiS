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
    class ( StorageForm ), allocatable :: &
      Field, &
      FieldOutput
  contains
    procedure, public, pass :: &
      InitializeTemplate_CSL
    procedure, public, pass :: &
      FinalizeTemplate_CSL
  end type Field_CSL_Template

contains


  subroutine InitializeTemplate_CSL &
               ( FC, C, NameShort, nValues, IgnorabilityOption )

    class ( Field_CSL_Template ), intent ( inout ) :: &
      FC
    class ( ChartTemplate ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nValues
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( FC % Type == '' ) &
      FC % Type = 'a Field_CSL' 

    FC % nValues = nValues

    call FC % InitializeTemplate ( C, NameShort, IgnorabilityOption )

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
