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
      nPolarZero
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
      ComputeKernel

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

    C % Fluid        =>  F

    C % nRadiusZero  =  1
    C % nPolarZero   =  1
    call PROGRAM_HEADER % GetParameter ( C % nRadiusZero, 'nRadiusZero' )    
    call PROGRAM_HEADER % GetParameter ( C % nPolarZero,  'nPolarZero' )    

  end subroutine Initialize_C_F


  subroutine Show_FS ( FS )

    class ( Coarsening_C_F_Form ), intent ( in ) :: &
      FS

    call FS % Coarsening_C_Form % Show ( )

    call Show ( FS % nRadiusZero, 'nRadiusZero' )
    call Show ( FS % nPolarZero, 'nPolarZero' )

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
      class is ( Atlas_SCG_CC_Form )
    associate &
      ( F        =>  C % Fluid, &
        C_GS_CC  =>  A % Chart_GS_CC )

    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_2 )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_3 )

    call C % Coarsening_C_Form % Compute ( FS )
    
    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .false. )

    call C_GS_CC % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )

    call ComputeKernel &
           ( FS_4D, &
             nC   = C_GS_CC % nCells, &
             oC   = C_GS_CC % nGhostLayers, &
             iaB  = C_GS_CC % iaBrick, &
             nRZ  = C % nRadiusZero, &
             nPZ  = C % nPolarZero, &
             iS_2 = iMomentum_2, &
             iS_3 = iMomentum_3, &
             UseDeviceOption = C % DeviceMemory )
    
    call FS % Storage_GS % ReassociateHost &
           ( AssociateVariablesOption = .true. )


    end associate !-- F, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F_Form', 'module', CONSOLE % ERROR )
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
