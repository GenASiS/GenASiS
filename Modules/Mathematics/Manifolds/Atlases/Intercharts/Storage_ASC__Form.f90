!-- Storage_ASC is generic storage for a set of fields on a single-chart 
!   Atlas.

module Storage_ASC__Form

  !-- Storage_AtlasSingleChart_Form

  use Basics
  use Charts
  use Field_ASC__Template
  use Atlas_SC__Form

  implicit none
  private

  type, public, extends ( Field_ASC_Template ) :: Storage_ASC_Form
    integer ( KDI ) :: &
      nFields
    logical ( KDL ) :: &
      Write
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    class ( Atlas_SC_Form ), pointer :: &
      Atlas_SC => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      Storage_CSL
    generic, public :: &
      Storage => Storage_CSL
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Storage_ASC_Form

contains


  subroutine Initialize &
               ( SA, A, NameShort, nFields, VariableOption, WriteOption, &
                 IgnorabilityOption )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      SA
    class ( Atlas_SC_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      NameShort
    integer ( KDI ), intent ( in ) :: &
      nFields
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    logical ( KDL ), intent ( in ), optional :: &
      WriteOption
    integer ( KDL ), intent ( in ), optional :: &
      IgnorabilityOption

    if ( SA % Type == '' ) &
      SA % Type = 'a Storage_ASC' 

    SA % nFields = nFields

    SA % Write = .false.
    if ( present ( WriteOption ) ) &
      SA % Write = WriteOption

    allocate ( SA % Variable ( nFields ) )
    SA % Variable = ''
    if ( present ( VariableOption ) ) &
      SA % Variable = VariableOption   

    SA % Atlas_SC => A

    call SA % InitializeTemplate_ASC ( A, NameShort, IgnorabilityOption )

  end subroutine Initialize


  function Storage_CSL ( SA ) result ( S )

    class ( Storage_ASC_Form ), intent ( in ) :: &
      SA
    class ( VariableGroupForm ), pointer :: &
      S

    select type ( SC => SA % Chart )
    class is ( Storage_CSL_Form )
      S => SC % Storage ( )
    class default
      call Show ( 'Storage Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Storage_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Storage_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SC

  end function Storage_CSL


  impure elemental subroutine Finalize ( SA )

    type ( Storage_ASC_Form ), intent ( inout ) :: &
      SA

    nullify ( SA % Atlas_SC )

    call SA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      FA

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Form )

    select type ( C => FA % Atlas_SC % Chart )
    class is ( Chart_SL_Template )

    allocate ( Storage_CSL_Form :: FA % Chart )

    select type ( SC => FA % Chart )
    class is ( Storage_CSL_Form )
      associate ( nValues => C % nProperCells + C % nGhostCells )
      call SC % Initialize &
                  ( C, FA % NameShort, FA % nFields, nValues, &
                    VariableOption = FA % Variable, &
                    WriteOption = FA % Write, &
                    IgnorabilityOption = FA % IGNORABILITY + 1 )
      end associate !-- nValues
    end select !-- GC

    if ( FA % Write ) &
      call A % AddField ( FA )

    end select !-- C
    end select !-- A

  end subroutine SetField


end module Storage_ASC__Form
