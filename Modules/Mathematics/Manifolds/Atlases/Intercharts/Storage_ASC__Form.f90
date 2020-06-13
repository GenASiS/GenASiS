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
      nFields = 0
    logical ( KDL ) :: &
      Write = .false.
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable
    class ( Atlas_SC_Form ), pointer :: &
      Atlas_SC => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeAllocate, InitializeClone
    procedure, public, pass :: &
      Storage_CSL
    procedure, public, pass :: &
      Storage
    final :: &
      Finalize
    procedure, private, pass :: &
      SetField
  end type Storage_ASC_Form

contains


  subroutine InitializeAllocate &
               ( SA, A, NameShort, nFields, VariableOption, WriteOption, &
                 UsePinnedMemoryOption, IgnorabilityOption )

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
      WriteOption, &
      UsePinnedMemoryOption
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

    call SA % InitializeTemplate_ASC &
           ( A, NameShort, UsePinnedMemoryOption, IgnorabilityOption )

  end subroutine InitializeAllocate


  subroutine InitializeClone &
               ( SA_Target, FA_Source, NameShort, UsePinnedMemoryOption, &
                 iaSelectedOption )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      SA_Target
    class ( Field_ASC_Template ), intent ( in ) :: &
      FA_Source
    character ( * ), intent ( in ) :: &
      NameShort
    logical ( KDL ), intent ( in ), optional :: &
      UsePinnedMemoryOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      iaSelectedOption
    
    logical ( KDL ) :: &
      UsePinnedMemory
      
    UsePinnedMemory = .false.
    if ( present ( UsePinnedMemoryOption ) ) &
      UsePinnedMemory = UsePinnedMemoryOption

    select type ( FC => FA_Source % Chart )
    class is ( Field_CSL_Template )

      if ( SA_Target % Type == '' ) &
        SA_Target % Type = 'a Storage_ASC' 

      allocate ( Storage_CSL_Form :: SA_Target % Chart )
      select type ( SC => SA_Target % Chart )
      type is ( Storage_CSL_Form )
      call SC % Initialize &
             ( FC, NameShort, UsePinnedMemory, &
               iaSelectedOption = iaSelectedOption )
      end select !-- SC

      call SA_Target % InitializeTemplate_ASC &
             ( FA_Source % Atlas, NameShort, &
               UsePinnedMemoryOption = UsePinnedMemory, &
               IgnorabilityOption = FA_Source % IGNORABILITY )

    class default
      call Show ( 'Field Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Storage_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeClone', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SC
    
  end subroutine InitializeClone


  function Storage_CSL ( SA ) result ( S )

    class ( Storage_ASC_Form ), intent ( in ) :: &
      SA
    class ( Storage_CSL_Form ), pointer :: &
      S

    select type ( SC => SA % Chart )
    class is ( Storage_CSL_Form )
      S => SC
    class default
      call Show ( 'Storage Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Storage_ASC__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Storage_CSL', 'function', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- SC

  end function Storage_CSL
  
  
  function Storage ( SA ) result ( S )

    class ( Storage_ASC_Form ), intent ( in ) :: &
      SA
    class ( StorageForm ), pointer :: &
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

  end function Storage


  impure elemental subroutine Finalize ( SA )

    type ( Storage_ASC_Form ), intent ( inout ) :: &
      SA

    nullify ( SA % Atlas_SC )

    call SA % FinalizeTemplate_ASC ( )

  end subroutine Finalize


  subroutine SetField ( FA )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      FA

    if ( allocated ( FA % Chart ) ) &
      return

    select type ( A => FA % Atlas )
    class is ( Atlas_SC_Form )

    select type ( C => FA % Atlas_SC % Chart )
    class is ( Chart_SL_Template )

    allocate ( Storage_CSL_Form :: FA % Chart )

    select type ( SC => FA % Chart )
    class is ( Storage_CSL_Form )
      associate ( nValues => C % nProperCells + C % nGhostCells )
      call SC % Initialize &
                  ( C, FA % NameShort, FA % UsePinnedMemory, FA % nFields, &
                    nValues, VariableOption = FA % Variable, &
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
