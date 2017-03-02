!-- Storage_CSL is generic storage for a set of fields on a Chart_SL.

module Storage_CSL__Form

  !-- Storage_ChartSingleLevel_Form

  use Basics
  use ChartBasics
  use Chart_SL__Template

  implicit none
  private

  type, public, extends ( Field_CSL_Template ) :: Storage_CSL_Form
    integer ( KDI ) :: &
      nFields
    class ( Chart_SL_Template ), pointer :: &
      Chart_SL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Storage_CSL_Form

contains


  subroutine Initialize ( SC, C, nFields, nValues, NameOutputOption )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      SC
    class ( Chart_SL_Template ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      nFields, &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOutputOption

    if ( SC % Type == '' ) &
      SC % Type = 'a Storage_CSL' 

    SC % nFields  =  nFields
    SC % Chart_SL => C

    call SC % InitializeTemplate_CSL ( C, nValues, NameOutputOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SC )

    type ( Storage_CSL_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % Chart_SL )

    call SC % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetField ( FC, NameOption )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      FC
    character ( * ), intent ( in ), optional :: &
      NameOption

    allocate ( FC % Field )

    call FC % Field % Initialize &
           ( [ FC % nValues, FC % nFields ], NameOption = NameOption )

  end subroutine SetField


end module Storage_CSL__Form
