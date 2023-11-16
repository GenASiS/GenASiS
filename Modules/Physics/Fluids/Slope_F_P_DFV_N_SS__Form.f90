module Slope_F_P_DFV_N_SS__Form

  !-- Slope_Fluid_Perfect_DivergenceFiniteVolume_Newton
  !   _SplitSource__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluid_P__Form
  use Slope_F_P_SS__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_F_P_DFV_N_SS_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F_P_DFV_N_SS
    generic, public :: &
      Initialize => InitializeAllocate_F_P_DFV_N_SS
    final :: &
      Finalize
  end type Slope_F_P_DFV_N_SS_Form


contains


  subroutine InitializeAllocate_F_P_DFV_N_SS ( S, RS, DF, DT, F )

    class ( Slope_F_P_DFV_N_SS_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F

    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_F_P_DFV_N_SS' 
    
    Name  =  trim ( F % Name ) // '_Slp_F_P_DFV_N_SS'

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY + 1 )

    associate ( nSC  =>  S % nComponents )

    !-- Slope component: Divergence + Newton gravity

    nSC  =  nSC + 1
    allocate ( Slope_DFV_N_Form :: S % Component ( nSC ) % Element )
    select type ( SD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_N_Form )

    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_B ( 1 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_B ( 2 ) )
    call Search ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_B ( 3 ) )
    call Search ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_B )

    call SD % Initialize &
           ( RS, DF, DT, &
             iVelocity_F   = F % VELOCITY_U, &
             iMomentum_B   = iMomentum_B, &
             iBaryonMass_F = F % BARYON_MASS, &
             iBaryonDensity_F = F % BARYON_DENSITY_B, &
             iEnergy_B = iEnergy_B )

    end select !-- SD

    !-- Slope component: Interactions

    nSC  =  nSC + 1
    allocate ( Slope_F_P_SS_Form :: S % Component ( nSC ) % Element )
    select type ( S_SS  =>  S % Component ( nSC ) % Element )
      class is ( Slope_F_P_SS_Form )

    call S_SS % Initialize ( F )

    end select !-- S_SS

    !-- Cleanup

    end associate !-- nSC

  end subroutine InitializeAllocate_F_P_DFV_N_SS


  impure elemental subroutine Finalize ( S )

    type ( Slope_F_P_DFV_N_SS_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_F_P_DFV_N_SS__Form
