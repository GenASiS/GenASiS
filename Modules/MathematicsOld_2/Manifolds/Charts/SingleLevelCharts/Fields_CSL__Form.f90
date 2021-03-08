module Fields_CSL__Form

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none
  private

  type, public, extends ( FieldsHeader_C_Form ) :: Fields_CSL_Form
    class ( StorageForm ), allocatable :: &
      Fields, &
      FieldsOutput
  contains
    procedure, private, pass :: &
      InitializeAllocate
    generic, public :: &
      Initialize => InitializeAllocate
    final :: &
      Finalize
  end type Fields_CSL_Form


contains


  subroutine InitializeAllocate &
               ( FC, FM, C, NameShort, nFields, FieldOption, WriteOption, &
                 PinnedOption, IgnorabilityOption )

    class ( Fields_CSL_Form ), intent ( inout ) :: &
      FC
    class ( FieldsHeader_M_Form ), intent ( in ), target :: &
      FM
    class ( ChartHeaderForm ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    logical ( KDL ), intent ( in ), optional :: &
      WriteOption, &
      PinnedOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    call FC % Initialize ( FM, C, NameShort, PinnedOption, IgnorabilityOption )

  end subroutine InitializeAllocate


  impure elemental subroutine Finalize ( FC )

    type ( Fields_CSL_Form ), intent ( inout ) :: &
      FC

    if ( allocated ( FC % FieldsOutput ) ) &
      deallocate ( FC % FieldsOutput )
    if ( allocated ( FC % Fields ) ) &
      deallocate ( FC % Fields )

  end subroutine Finalize


end module Fields_CSL__Form
