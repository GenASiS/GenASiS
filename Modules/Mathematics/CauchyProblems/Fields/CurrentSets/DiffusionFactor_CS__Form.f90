module DiffusionFactor_CS__Form

  !-- DiffusionFactor_CurrentSet__Form

  use Basics
  use Manifolds
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
    procedure, public, pass ( DF ) :: &
      ComputeReconstruction
    final :: &
      Finalize
  end type DiffusionFactor_CS_Form

    private :: &
      ComputeKernel, &
      Compute_I_CGS_Kernel

    interface

      module subroutine ComputeKernel ( DF, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          DF
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

      module subroutine Compute_I_CGS_Kernel &
               ( DF, iD, oV, DF_I, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
          DF
        integer ( KDI ), intent ( in ) :: &
          iD, &
          oV   
        real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
          DF_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_I_CGS_Kernel

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
      DF % Type  =  'a DiffusionFactor_CS' 
    
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


  subroutine ComputeReconstruction ( FS, DF, iDF, iC, iD )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    integer ( KDI ), intent ( in ) :: &
      iDF, &  !-- iDiffusionFactor
      iC, &   !-- iChart
      iD      !-- iDimensions

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      dX_3D, &
      DF_3D, &
      DF_I_3D

    call Show ( 'Reconstructing ' // trim ( DF % Type ), DF % IGNORABILITY + 3 )
    call Show ( DF % Name, 'Name', DF % IGNORABILITY + 3 )

    associate &
      ( DFV  =>  DF % Storage ( iC ) % Value, &
        FSV  =>  FS % Storage ( iC ) % Value )

    select type ( C  =>  FS % Atlas % Chart ( iC ) % Element )
    class is ( Chart_GS_Form )

      call C % SetFieldPointer ( DFV ( :, DF % DIFFUSION_FACTOR ), DF_3D )
      call C % SetFieldPointer ( FSV ( :, iDF ), DF_I_3D )

      call Compute_I_CGS_Kernel &
             ( DF_3D, iD, C % nGhostLayers ( iD ), DF_I_3D, &
               UseDeviceOption = DF % DeviceMemory )

    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Reconstruction_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C
    
    end associate !-- DFV, etc.

  end subroutine ComputeReconstruction


  impure elemental subroutine Finalize ( DF )

    type ( DiffusionFactor_CS_Form ), intent ( inout ) :: &
      DF

    nullify ( DF % CurrentSet )

  end subroutine Finalize


end module DiffusionFactor_CS__Form
