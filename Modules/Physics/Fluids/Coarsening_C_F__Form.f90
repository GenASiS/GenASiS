module Coarsening_C_F__Form

  !-- Coarsening_Central_Fluid_Form

  use Basics
  use Mathematics
  use Fluid_D__Form

  implicit none
  private

  type, public, extends ( Coarsening_C_Form ) :: Coarsening_C_F_Form
    integer ( KDI ) :: &
      nRadiusZero
    integer ( KDI ), dimension ( : ), allocatable :: &
      nPolarZero
    class ( Fluid_D_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      Initialize_C_F
    generic, public :: &
      Initialize => Initialize_C_F
    procedure, public, pass ( C ) :: &
      Compute
    final :: &
      Finalize
  end type Coarsening_C_F_Form

    private :: &
      SetRadiusZero, &
      SetPolarZero

    private :: &
      ComputeKernel

    interface
      
      module subroutine ComputeKernel &
               ( FS_4D, nPZ, nC, oC, iaB, nRZ, iS_2, iS_3, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, :, :, : ), intent ( inout ) :: &
          FS_4D
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          nPZ, &
          nC, &
          oC, &
          iaB
        integer ( KDI ), intent ( in ) :: &
          nRZ, &
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

    C % Fluid  =>  F

    select type ( A  =>  G % Atlas )
    class is ( Atlas_SCG_CC_Form )

      call SetRadiusZero &
             (  C  = A % Chart_GS_CC, &
                CP = C % Storage_GS % Value ( :, C % COARSENING_POLAR ), &
               nRZ = C % nRadiusZero )
      call SetPolarZero &
             (  C  = A % Chart_GS_CC, &
                CA = C % Storage_GS % Value ( :, C % COARSENING_AZIMUTHAL ), &
               nPZ = C % nPolarZero )
!call Show ( C % nRadiusZero, '>>> nRadiusZero' )
!call Show ( C % nPolarZero, '>>> nPolarZero' )
!call A % Chart_GS_CC % Communicator % Synchronize ( )
!call PROGRAM_HEADER % Abort ( )

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Coarsening_C_F__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_C_F', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

  end subroutine Initialize_C_F


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

    call C_GS_CC % SetFieldPointer &
           ( FS % Storage_GS % Value, FS_4D )

    call ComputeKernel &
           ( FS_4D, &
             nPZ  = C % nPolarZero, &
             nC   = C_GS_CC % nCells, &
             oC   = C_GS_CC % nGhostLayers, &
             iaB  = C_GS_CC % iaBrick, &
             nRZ  = C % nRadiusZero, &
             iS_2 = iMomentum_2, &
             iS_3 = iMomentum_3, &
             UseDeviceOption = C % DeviceMemory )

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

    if ( allocated ( C % nPolarZero ) ) &
      deallocate ( C % nPolarZero )

  end subroutine Finalize


  subroutine SetRadiusZero ( C, CP, nRZ )

    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      CP  !-- CoarsenPolar
    integer ( KDI ), intent ( out ) :: &
      nRZ  !-- nRadiusZero

    integer ( KDI ) :: &
      iR, &  !-- iRadius
      nCP    !-- nCoarsenPolar
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      CP_3D

    nRZ  =  0

    if ( C % nDimensions  <  2 )  &
      return

    call C % SetFieldPointer ( CP, CP_3D )

    associate &
      ( nR   =>  C % nCellsBrick ( 1 ), &
        nTh  =>  C % nCellsBrick ( 2 ) )  !-- nCellsBrick ( 2 ) == nCells ( 2 )

    do iR  =  1,  nR
      nCP  =  CP_3D ( iR, 1, 1 )  +  0.5_KDR
      if ( nCP  ==  nTh ) &
        nRZ  =  nRZ + 1
    end do !-- iR

    end associate !-- nR, etc.

  end subroutine SetRadiusZero


  subroutine SetPolarZero ( C, CA, nPZ )

    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      CA  !-- CoarsenAzimuthal
    integer ( KDI ), dimension ( : ), intent ( out ), allocatable :: &
      nPZ  !-- nPolarZero

    integer ( KDI ) :: &
      iR, iTh, &  !-- iRadius, iTheta
      nCA    !-- nCoarsenAzimuthal
    integer ( KDI ), dimension ( : ), allocatable :: &
      nPZ_Temp
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      CA_3D

    if ( C % nDimensions  <  3 ) then
      allocate ( nPZ ( 0 ) )
      return
    end if

    call C % SetFieldPointer ( CA, CA_3D )

    associate &
      ( nR   =>  C % nCellsBrick ( 1 ), &
        nTh  =>  C % nCellsBrick ( 2 ), &  !-- nCellsBrick ( 2 ) == nCells ( 2 )
        nPh  =>  C % nCellsBrick ( 3 ) )   !-- nCellsBrick ( 3 ) == nCells ( 3 )

    allocate ( nPZ_Temp ( nR ) )
    nPZ_Temp  =  0

    do iR  =  1,  nR
      do iTh  =  1,  nTh / 2
        nCA  =  CA_3D ( iR, iTh, 1 )  +  0.5_KDR
        if ( nCA  ==  nPh ) &
          nPZ_Temp ( iR )  =  nPZ_Temp ( iR )  +  1
      end do !-- nTh / 2
    end do !-- iR

    allocate ( nPZ, source = pack ( nPZ_Temp, nPZ_Temp > 0 ) )

    end associate !-- nR, etc.

  end subroutine SetPolarZero


end module Coarsening_C_F__Form
