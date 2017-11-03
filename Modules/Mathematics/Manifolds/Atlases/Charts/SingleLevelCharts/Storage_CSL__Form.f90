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
    logical ( KDL ) :: &
      Write
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    class ( Chart_SL_Template ), pointer :: &
      Chart_SL => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Storage
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Storage_CSL_Form

contains


  subroutine Initialize &
               ( SC, C, NameShort, nFields, nValues, VariableOption, &
                 WriteOption, IgnorabilityOption )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      SC
    class ( Chart_SL_Template ), intent ( in ), target :: &
      C
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields, &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    logical ( KDL ), intent ( in ), optional :: &
      WriteOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( SC % Type == '' ) &
      SC % Type = 'a Storage_CSL' 

    SC % nFields = nFields

    SC % Write = .false.
    if ( present ( WriteOption ) ) &
      SC % Write = WriteOption

    allocate ( SC % Variable ( nFields ) )
    SC % Variable = ''
    if ( present ( VariableOption ) ) &
      SC % Variable = VariableOption   

    SC % Chart_SL => C

    call SC % InitializeTemplate_CSL &
           ( C, NameShort, nValues, IgnorabilityOption )

  end subroutine Initialize


  function Storage ( SC ) result ( S )

    class ( Storage_CSL_Form ), intent ( in ), target :: &
      SC
    class ( VariableGroupForm ), pointer :: &
      S
      
    class ( VariableGroupForm ), pointer :: &
      Field

    S => null ( )

    Field => SC % Field
    select type ( Field )
    class is ( VariableGroupForm )
    S => Field
    end select !-- Field

  end function Storage


  impure elemental subroutine Finalize ( SC )

    type ( Storage_CSL_Form ), intent ( inout ) :: &
      SC

    nullify ( SC % Chart_SL )

    if ( allocated ( SC % Variable ) ) &
      deallocate ( SC % Variable )

    call SC % FinalizeTemplate_CSL ( )

  end subroutine Finalize


  subroutine SetField ( FC )

    class ( Storage_CSL_Form ), intent ( inout ) :: &
      FC

    allocate ( FC % Field )

    call FC % Field % Initialize &
           ( [ FC % nValues, FC % nFields ], &
             VariableOption = FC % Variable, NameOption = FC % NameShort )

    if ( FC % Write ) then
      allocate ( FC % FieldOutput )
      call FC % FieldOutput % Initialize ( FC % Field )
    end if

  end subroutine SetField


end module Storage_CSL__Form
