module DiffusionFactor_CS__Form

  !-- DiffusionFactor_CurrentSet__Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_CS  = 1, &
      N_VECTORS_CS = 0

  type, public, extends ( FieldSet_BM_Form ) :: DiffusionFactor_CS_Form
    integer ( KDI ) :: &
      N_FIELDS_CS  = N_FIELDS_CS, &
      N_VECTORS_CS = N_VECTORS_CS
    integer ( KDI ) :: &
      DIFFUSION_FACTOR = 0
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_DF
    generic, public :: &
      Initialize => InitializeAllocate_DF
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type DiffusionFactor_CS_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel ( DF, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          DF
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_DF &
               ( DF, CS, FieldOption, nFieldsOption )

    class ( DiffusionFactor_CS_Form ), intent ( inout ) :: &
      DF
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
 
    integer ( KDI ) :: &
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( DF % Type  ==  '' ) &
      DF % Type  =  'an DiffusionFactor_CS' 
    
    Name  =  trim ( CS % Name ) // '_DffsnFctr'

    DF % CurrentSet   =>  CS

    !-- Field indices

    DF % DIFFUSION_FACTOR  =  1

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      nFields  =  DF % N_FIELDS_CS
    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : DF % N_FIELDS_CS ) &
      =  [ 'DiffusionFactor' ]
          
    !-- FieldSet

    call DF % FieldSet_BM_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_DF


  subroutine Compute ( DF, iC, iD )

    class ( DiffusionFactor_CS_Form ), intent ( inout ) :: &
      DF
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions

    call Show ( 'Computing ' // trim ( DF % Type ), DF % IGNORABILITY + 3 )
    call Show ( DF % Name, 'Name', DF % IGNORABILITY + 3 )

    associate &
      ( DFV  =>  DF % Storage ( iC ) % Value )

      call ComputeKernel &
             ( DFV ( : , DF % DIFFUSION_FACTOR ), &
               UseDeviceOption = DF % DeviceMemory )

    end associate !-- DFV

  end subroutine Compute


  ! subroutine ComputeDiffusionFactor ( DF, DP, iC, iD )

  !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  !     DF
  !   class ( DivergencePart_CS_Form ), intent ( in ) :: &
  !     DP
  !   integer ( KDI ), intent ( in ) :: &
  !     iC, &  !-- iChart
  !     iD     !-- iDimension

  !   associate ( CS  =>  DP % CurrentSet )

  !     call Compute_DF_Kernel ( DF, UseDeviceOption = CS % DeviceMemory )

  !   end associate !-- CS

  ! end subroutine ComputeDiffusionFactor


  impure elemental subroutine Finalize ( DF )

    type ( DiffusionFactor_CS_Form ), intent ( inout ) :: &
      DF

    nullify ( DF % CurrentSet )

  end subroutine Finalize


end module DiffusionFactor_CS__Form
