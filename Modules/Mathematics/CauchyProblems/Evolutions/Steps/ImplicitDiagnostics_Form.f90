module ImplicitDiagnostics_Form

  use Basics
  use Manifolds
  use Fields

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_ID = 2

  type, public, extends ( FieldSet_BM_Form ) :: ImplicitDiagnosticsForm
    integer ( KDI ) :: &
      N_FIELDS_ID = N_FIELDS_ID
    integer ( KDI ) :: &
      N_ITERATIONS = 0, &
      RESIDUAL_MAX = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_ID
    generic, public :: &
      Initialize => InitializeAllocate_ID
    final :: &
      Finalize
  end type ImplicitDiagnosticsForm


contains


  subroutine InitializeAllocate_ID &
               ( FS, A, iStage, FieldOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, nFieldsOption, IgnorabilityOption )

    class ( ImplicitDiagnosticsForm ), intent ( inout ), target :: &
      FS
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iStage
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption, &
      DevicesCommunicateOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    integer ( KDI ) :: &
      nFields
    character ( 1 ) :: &
      StageNumber
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'an ImplicitDiagnostics' 

    write ( StageNumber, fmt = '(i1.1)' ) iStage    
    Name  =  'ImplicitDiagnostics_Stage_' // StageNumber
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else

      FS % N_ITERATIONS  =  1
      FS % RESIDUAL_MAX  =  2

      nFields  =  FS % N_FIELDS_ID

    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
      Field ( 1 )  =  'nIterations'
      Field ( 2 )  =  'ResidualMax'
    end if !-- FieldOption

    !-- FieldSet

    call FS % FieldSet_BM_Form % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_ID


  impure elemental subroutine Finalize ( ID )

    type ( ImplicitDiagnosticsForm ), intent ( inout ) :: &
      ID

  end subroutine Finalize


end module ImplicitDiagnostics_Form
