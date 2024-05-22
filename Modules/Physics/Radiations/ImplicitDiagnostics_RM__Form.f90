module ImplicitDiagnostics_RM__Form
  
  !-- ImplicitDiagnostics_RadiationMoments_Form

  use Basics
  use Mathematics

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_RM = 2

  type, public, extends ( ImplicitDiagnosticsForm ) &
    :: ImplicitDiagnostics_RM_Form
      integer ( KDI ) :: &
        N_FIELDS_RM = N_FIELDS_RM
      integer ( KDI ) :: &
        RESIDUAL_RADIATION_ENERGY = 0, &
        RESIDUAL_FLUID_ENERGY     = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_ID
    final :: &
      Finalize
  end type ImplicitDiagnostics_RM_Form


contains


  subroutine InitializeAllocate_ID &
               ( ID, A, FieldSetName, iStage, FieldOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, nFieldsOption, IgnorabilityOption )

    class ( ImplicitDiagnostics_RM_Form ), intent ( inout ), target :: &
      ID
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A
    character ( * ), intent ( in ) :: &
      FieldSetName
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
      oF, &
      nFields
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( ID % Type  ==  '' ) &
      ID % Type  =  'an ImplicitDiagnostics_RM' 

    !-- Field indices

    oF  =  ID % N_FIELDS_ID

    ID % RESIDUAL_RADIATION_ENERGY  =  oF  +  1
    ID % RESIDUAL_FLUID_ENERGY      =  oF  +  2

    nFields  =  oF  +  ID % N_FIELDS_RM
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + ID % N_FIELDS_RM ) &
      = [ 'ResidualRadiationEnergy', &
          'ResidualFluidEnergy    ' ]
          
    !-- FieldSet

    call ID % ImplicitDiagnosticsForm % Initialize &
           ( A, FieldSetName, iStage, &
             FieldOption = Field, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_ID


  impure elemental subroutine Finalize ( ID )

    type ( ImplicitDiagnostics_RM_Form ), intent ( inout ) :: &
      ID

  end subroutine Finalize


end module ImplicitDiagnostics_RM__Form
