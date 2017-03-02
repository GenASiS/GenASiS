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


  subroutine Initialize &
               ( SC, C, NameShort, nFields, nValues, IgnorabilityOption )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      SC
    class ( Chart_SL_Template ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields, &
      nValues
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( SC % Type == '' ) &
      SC % Type = 'a Storage_CSL' 

    SC % nFields  =  nFields
    SC % Chart_SL => C

    call SC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  impure elemental subroutine Finalize ( SC )

    type ( Storage_CSL_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % Chart_SL )

    call SC % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % Field )

    call FC % Field % Initialize &
           ( [ FC % nValues, FC % nFields ], NameOption = FC % NameShort )

  end subroutine SetField


end module Storage_CSL__Form
