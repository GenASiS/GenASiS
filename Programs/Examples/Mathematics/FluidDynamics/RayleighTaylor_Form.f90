module RayleighTaylor_Form

  use Basics
  use Mathematics
  use Fluid_P_P__Form
  use Sources_F__Form
  use Fluid_ASC__Form
  use Tally_RT__Form
  
  implicit none
  private
  
  type, public, extends ( Integrator_C_PS_Form ) :: RayleighTaylorForm
  contains
    procedure, private, pass :: &
      Initialize_RT
    generic, public :: &
      Initialize => Initialize_RT
    final :: &
      Finalize  
  end type RayleighTaylorForm

    private :: &
      SetFluid, & 
      ApplySources
    
    real ( KDR ), private :: &
      DensityAbove, DensityBelow, &
      PressureBase, &   
      AdiabaticIndex, &
      Acceleration

    type ( UniverseHeaderForm ), private :: &
      Universe  !-- Non-functional dummy argument
      
    interface
    
      module subroutine ApplySourcesKernel &
               ( KVM, KVE, SVM, SVE, N, VY, dT, A, B, &
                 UseDeviceOption )
        use Basics
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          KVM, &
          KVE, &
          SVM, &
          SVE
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          N, &
          VY
        real ( KDR ), intent ( in ) :: &
          dT, &
          A, &
          B
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ApplySourcesKernel  
    
    end interface

contains

  
  subroutine Initialize_RT ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iB  !-- iBoundary
    class ( GeometryFlatForm ), pointer :: &
      G  

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: RT % PositionSpace )
    select type ( PS => RT % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    call PS % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 1 )
    call PS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )

    call PS % CreateChart &
           ( MinCoordinateOption = [ -0.25_KDR, -0.75_KDR ], &
             MaxCoordinateOption = [ +0.25_KDR, +0.75_KDR ], &
             nCellsOption = [ 64, 192 ] )

    call PS % SetGeometry ( UsePinnedMemoryOption = .true. )
    call PS % Geometry_ASC % AllocateDevice ( )
    G => PS % Geometry ( )
    call G % UpdateDevice ( )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: RT % Current_ASC )
    select type ( FA => RT % Current_ASC )
    class is ( Fluid_ASC_Form )

    allocate ( Tally_RT_Form :: FA % TallyInterior )
    allocate ( Tally_RT_Form :: FA % TallyTotal )
    allocate ( Tally_RT_Form :: FA % TallyChange )
    allocate ( FA % TallyBoundaryLocal ( PS % nBoundaries ) )
    allocate ( FA % TallyBoundaryGlobal ( PS % nBoundaries ) )
    do iB = 1, PS % nBoundaries 
      allocate ( Tally_RT_Form :: FA % TallyBoundaryLocal ( iB ) % Element )
      allocate ( Tally_RT_Form :: FA % TallyBoundaryGlobal ( iB ) % Element )
    end do !-- iB
       
    call FA % Initialize ( PS, 'POLYTROPIC', UsePinnedMemoryOption = .true. )
    call FA % AllocateDevice ( )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: RT % Step )
    select type ( S => RT % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( RT, FA, Name )
    S % ApplySources % Pointer => ApplySources
    end select !-- S

    !-- Problem definition

    DensityAbove   = 2.0_KDR
    DensityBelow   = 1.0_KDR
    PressureBase   = 2.5_KDR
    AdiabaticIndex = 1.4_KDR
    Acceleration   = 0.1_KDR

    call PROGRAM_HEADER % GetParameter ( DensityAbove, 'DensityAbove' )
    call PROGRAM_HEADER % GetParameter ( DensityBelow, 'DensityBelow' )
    call PROGRAM_HEADER % GetParameter ( AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( Acceleration, 'Acceleration' )
    
    !-- Set fluid and initialize template

    call SetFluid ( RT )
    call RT % Initialize ( Universe, Name, FinishTimeOption = 8.5_KDR )

    !-- Cleanup
    
    nullify ( G )
    
    end select !-- FA
    end select !-- PS

  end subroutine Initialize_RT

  
  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

  end subroutine Finalize


  subroutine SetFluid ( RT )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    integer ( KDI ) :: &
      iB  !-- iBoundary
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_P_Form ), pointer :: &
      F
      
    real ( KDR ) :: &
      Pi

    select type ( FA => RT % Current_ASC )
    class is ( Fluid_ASC_Form )

    !-- Initial conditions
       
    F => FA % Fluid_P_P ( )

    call F % SetAdiabaticIndex ( AdiabaticIndex )

    select type ( PS => RT % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      ( X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        N  => F % Value ( :, F % COMOVING_DENSITY ), &
        E  => F % Value ( :, F % INTERNAL_ENERGY ), &
        VY => F % Value ( :, F % VELOCITY_U ( 2 ) ) )
    
    !-- FIXME: Association to CONSTANT singleton gives wrong value
    !          with optimization on XL
    Pi = CONSTANT % PI 

    where ( Y > 0.0_KDR )     
      N  = DensityAbove
    elsewhere
      N  = DensityBelow
    end where

    E  =  ( PressureBase  -  N * Acceleration * Y ) &
          / ( AdiabaticIndex - 1.0_KDR )

    VY  =  ( 0.01_KDR / 4.0_KDR ) &
           *  ( 1.0_KDR  +  cos ( 4.0_KDR * Pi * X ) ) &
           *  ( 1.0_KDR  +  cos ( 3.0_KDR * Pi * Y ) )
    
    call F % UpdateDevice ( )
    call G % UpdateDevice ( )
    call F % ComputeFromPrimitive ( G )
    call F % UpdateHost ( )

    end associate !-- X, etc.
    end select !-- PS
    nullify ( F, G )

    !-- Tally

    select type ( TI => FA % TallyInterior )
    class is ( Tally_RT_Form )
       call TI % SetAcceleration ( Acceleration )
    end select !-- TI
    
    select type ( TT => FA % TallyTotal )
    class is ( Tally_RT_Form )
       call TT % SetAcceleration ( Acceleration )
    end select !-- TT
    
    select type ( TC => FA % TallyChange )
    class is ( Tally_RT_Form )
       call TC % SetAcceleration ( Acceleration )
    end select !-- TC
    
    do iB = 1, size ( FA % TallyBoundaryLocal )
      select type ( TB => FA % TallyBoundaryLocal ( iB ) % Element )
      class is ( Tally_RT_Form )
        call TB % SetAcceleration ( Acceleration )
      end select !-- TB
      select type ( TB => FA % TallyBoundaryGlobal ( iB ) % Element )
      class is ( Tally_RT_Form )
        call TB % SetAcceleration ( Acceleration )
      end select !-- TB
    end do !-- iB

    end select !-- FA
    
  end subroutine SetFluid

  
  subroutine ApplySources &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( StorageForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iMomentum_2, &
      iEnergy   
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( F => Fluid )
    class is ( Fluid_P_P_Form )

    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )

    associate &
      ( KVM => Increment % Value ( :, iMomentum_2 ), &
        KVE => Increment % Value ( :, iEnergy ), &
        N   => F % Value ( :, F % COMOVING_DENSITY ), &
        VY  => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        A   => Acceleration, &
        dT  => TimeStep )
    
    !KVM  =  KVM  -  dT * N * A
    !KVE  =  KVE  -  dT * N * A * VY

    select type ( FS => Sources_F )
    class is ( Sources_F_Form )
      associate &
        ( SVM => FS % Value ( :, FS % GRAVITATIONAL_S_D ( 2 ) ), &
          SVE => FS % Value ( :, FS % GRAVITATIONAL_G ) )
      !SVM  =  SVM  -  S % B ( iStage ) * N * A 
      !SVE  =  SVE  -  S % B ( iStage ) * N * A * VY 
      
      !call Increment % UpdateHost ( )
      !call Show &
      !      ( Increment % Value ( :, iMomentum_2 ), &
      !        '==== Increment iMomentum2 A ====' )
              
      call ApplySourcesKernel &
             ( Increment % Value ( :, iMomentum_2 ), &
               Increment % Value ( :, iEnergy ), &
               FS % Value ( :, FS % GRAVITATIONAL_S_D ( 2 ) ), &
               FS % Value ( :, FS % GRAVITATIONAL_G ), &
               F % Value ( :, F % COMOVING_DENSITY ), &
               F % Value ( :, F % VELOCITY_U ( 2 ) ), &
               TimeStep, Acceleration, S % B ( iStage ), &
               UseDeviceOption = F % AllocatedDevice )
      
      !call Increment % UpdateHost ( )
      !call Show &
      !      ( Increment % Value ( :, iMomentum_2 ), &
      !        '==== Increment iMomentum2 B ====' )
      
      end associate !-- SVM, etc.
    end select !-- FS

    end associate !-- KVM, etc.
    end select !-- F

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplySources

  
end module RayleighTaylor_Form
