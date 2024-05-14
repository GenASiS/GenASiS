module ImplicitDiagnostics_NM_G__Form
  
  !-- ImplicitDiagnostics_NeutrinoMoments_Grey_Form

  use Basics
  use Mathematics
  use ImplicitDiagnostics_RM__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_NM_G = 2

  type, public, extends ( ImplicitDiagnostics_RM_Form ) &
    :: ImplicitDiagnostics_NM_G_Form
      integer ( KDI ) :: &
        N_FIELDS_NM_G = N_FIELDS_NM_G
      integer ( KDI ) :: &
        RESIDUAL_RADIATION_NUMBER = 0, &
        RESIDUAL_FLUID_NUMBER     = 0
  contains
    procedure, private, pass :: &
      InitializeAllocate_ID
    final :: &
      Finalize
  end type ImplicitDiagnostics_NM_G_Form


contains


  subroutine InitializeAllocate_ID &
               ( ID, A, iStage, FieldOption, NameOption, &
                 DeviceMemoryOption, PinnedMemoryOption, &
                 DevicesCommunicateOption, nFieldsOption, IgnorabilityOption )

    class ( ImplicitDiagnostics_NM_G_Form ), intent ( inout ), target :: &
      ID
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
      oF, &
      nFields
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( ID % Type  ==  '' ) &
      ID % Type  =  'an ImplicitDiagnostics_NM_G' 

    !-- Field indices

    oF  =  ID % N_FIELDS_ID  +  ID % N_FIELDS_RM

    ID % RESIDUAL_RADIATION_NUMBER  =  oF  +  1
    ID % RESIDUAL_FLUID_NUMBER      =  oF  +  2

    nFields  =  oF  +  ID % N_FIELDS_NM_G
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + ID % N_FIELDS_NM_G ) &
      = [ 'ResidualRadiationNumber', &
          'ResidualFluidNumber    ' ]
          
    !-- FieldSet

    call ID % ImplicitDiagnostics_RM_Form % Initialize &
           ( A, iStage, &
             FieldOption = Field, &
             DeviceMemoryOption = DeviceMemoryOption, &
             PinnedMemoryOption = PinnedMemoryOption, &
             DevicesCommunicateOption = DevicesCommunicateOption, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_ID


  impure elemental subroutine Finalize ( ID )

    type ( ImplicitDiagnostics_NM_G_Form ), intent ( inout ) :: &
      ID

  end subroutine Finalize


end module ImplicitDiagnostics_NM_G__Form
