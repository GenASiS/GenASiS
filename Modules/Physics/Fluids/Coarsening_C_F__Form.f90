module Coarsening_C_F__Form

  !-- Coarsening_Central_Fluid_Form

  use Basics
  use Mathematics
  use Fluid_D__Form

  implicit none
  private

  type, public, extends ( Coarsening_C_Form ) :: Coarsening_C_F_Form
    integer ( KDI ) :: &
      nRadiusZero, &
      nPolarZero!, &
!      nBlocksThreshold
    real ( KDR ) :: &
      RadiusZero = 0.0_KDR
    class ( Fluid_D_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      Initialize_C_F
    generic, public :: &
      Initialize => Initialize_C_F
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass ( C ) :: &
      Compute
    final :: &
      Finalize
  end type Coarsening_C_F_Form

    private :: &
      ComputeKernel, &
      ComputeBlocksKernel

    interface
      
      module subroutine ComputeKernel &
               ( FS_4D, nC, oC, iaB, nRZ, nPZ, iS_2, iS_3, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          FS_4D
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nC, &
          oC, &
          iaB
        integer ( KDI ), intent ( in ) :: &
          nRZ, nPZ, &
          iS_2, iS_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

      module subroutine ComputeBlocksKernel &
               ( FS, BP, BA, nBT, iS_2, iS_3, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          FS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          BP, &
          BA
        integer ( KDI ), intent ( in ) :: &
          nBT, &
          iS_2, iS_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeBlocksKernel

      module subroutine ComputeRadiusKernel &
               ( FS, R, RZ, iS_2, iS_3, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          FS
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          R
        real ( KDR ), intent ( in ) :: &
          RZ
        integer ( KDI ), intent ( in ) :: &
          iS_2, iS_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeRadiusKernel

    end interface


contains


  subroutine Initialize_C_F ( C, F, G )

    class ( Coarsening_C_F_Form ), intent ( inout ) :: &
      C
    class ( Fluid_D_Form ), intent ( in ), target :: &
      F
    class ( Geometry_F_Form ), intent ( in ) :: &
      G

    if ( C % Type == '' ) &
      C % Type  =  'a Coarsening_C_F' 
    
    call C % Coarsening_C_Form % Initialize ( G )

    C % Fluid  =>  F

    select type ( A  =>  C % Atlas )
      class is ( Atlas_SCG_C_Form )

    select type ( C_GS_C  =>  A % Chart_GS_C )
    class is ( Chart_GS_CE_Form )
      C % nRadiusZero  =  0
      C % nPolarZero   =  1
!      C % nBlocksThreshold  =  8
      C % RadiusZero  =  0.0_KDR
      call PROGRAM_HEADER % GetParameter ( C % RadiusZero, 'RadiusZero' )
      call PROGRAM_HEADER % GetParameter &
             ( C % nRadiusZero, 'nRadiusZero' )    
      call PROGRAM_HEADER % GetParameter &
             ( C % nPolarZero, 'nPolarZero' )    
!      call PROGRAM_HEADER % GetParameter &
!             ( C % nBlocksThreshold, 'nBlocksThreshold' )    
    class is ( Chart_GS_CC_Form )
      C % nRadiusZero  =  1
      C % nPolarZero   =  1
!      C % nBlocksThreshold  =  8
      C % RadiusZero  =  C_GS_C % RadiusCore  /  10.0_KDR    
      call PROGRAM_HEADER % GetParameter ( C % RadiusZero, 'RadiusZero' )
      call PROGRAM_HEADER % GetParameter &
             ( C % nRadiusZero, 'nRadiusZero' )    
      call PROGRAM_HEADER % GetParameter &
             ( C % nPolarZero, 'nPolarZero' )    
!      call PROGRAM_HEADER % GetParameter &
!             ( C % nBlocksThreshold, 'nBlocksThreshold' )    
    end select !-- C_GS_C

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_C_F', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize_C_F


  subroutine Show_FS ( FS )

    class ( Coarsening_C_F_Form ), intent ( in ) :: &
      FS

    call FS % Coarsening_C_Form % Show ( )

    call Show ( FS % nRadiusZero, 'nRadiusZero' )
    call Show ( FS % nPolarZero, 'nPolarZero' )
!    call Show ( FS % nBlocksThreshold, 'nBlocksThreshold' )
    
    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_C_Form )
    associate &
      ( C_GS_C  =>  A % Chart_GS_C )

    call Show ( FS % RadiusZero, C_GS_C % CoordinateUnit ( 1 ), 'RadiusZero' )

    end associate !-- C_GS_C
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Show_FS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Show_FS


  subroutine Compute ( FS, C )

    class ( FieldSetForm ), intent ( inout ) :: &
      FS
    class ( Coarsening_C_F_Form ), intent ( in ) :: &
      C

    integer ( KDI ) :: &
      iMomentum_2, &
      iMomentum_3
    real ( KDR ), dimension ( :, :, :, : ), pointer :: &
      FS_4D

    select type ( A  =>  FS % Atlas )
      class is ( Atlas_SCG_C_Form )
    associate &
      ( F       =>  C % Fluid, &
        G       =>  C % Geometry, &
        C_GS_C  =>  A % Chart_GS_C )

    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_2 )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_3 )

    call C % Coarsening_C_Form % Compute ( FS )
    
    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .false. )

    call C_GS_C % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )

    call ComputeKernel &
           ( FS_4D, &
             nC   = C_GS_C % nCells, &
             oC   = C_GS_C % nGhostLayers, &
             iaB  = C_GS_C % iaBrick, &
             nRZ  = C % nRadiusZero, &
             nPZ  = C % nPolarZero, &
             iS_2 = iMomentum_2, &
             iS_3 = iMomentum_3, &
             UseDeviceOption = FS % DeviceMemory )
    
    ! call ComputeBlocksKernel &
    !        (  FS  = FS % Storage_GS % Value, &
    !           BP  = C  % Storage_GS % Value ( :, C % N_BLOCKS_POLAR ), &
    !           BA  = C  % Storage_GS % Value ( :, C % N_BLOCKS_AZIMUTHAL ), &
    !          nBT  = C  % nBlocksThreshold, &
    !          iS_2 = iMomentum_2, &
    !          iS_3 = iMomentum_3, &
    !          UseDeviceOption = C % DeviceMemory )

    call ComputeRadiusKernel &
           (  FS  = FS % Storage_GS % Value, &
              R   = G  % Storage_GS % Value ( :, G % CENTER_U_1 ), &
              RZ  = C % RadiusZero, &
             iS_2 = iMomentum_2, &
             iS_3 = iMomentum_3, &
             UseDeviceOption = FS % DeviceMemory )

    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .true. )

    end associate !-- F, etc.
    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Compute


  impure elemental subroutine Finalize ( C )

    type ( Coarsening_C_F_Form ), intent ( inout ) :: &
      C

    nullify ( C % Fluid )

  end subroutine Finalize

  
end module Coarsening_C_F__Form
