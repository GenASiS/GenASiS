module Slope_DFV_N_F_P_HN__Form

  !-- Slope_DivergenceFiniteVolume_Newton_Fluid_Perfect_HeavyNucleus__Form

  use Basics
  use Mathematics
  use Gravitations
  use Slope_DFV_F_F_P_HN__Form

  implicit none
  private

  type, public, extends ( Slope_DFV_N_Form ) :: Slope_DFV_N_F_P_HN_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_N
    final :: &
      Finalize
  end type Slope_DFV_N_F_P_HN_Form


contains


  subroutine InitializeAllocate_N &
               ( S, RS, DP_1D, iVelocity_F, iMomentum_B, iBaryonMass_F, &
                 iBaryonDensity_F, iEnergy_B, SuffixOption )

    class ( Slope_DFV_N_F_P_HN_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    type ( DivergencePartElement ), dimension ( : ), intent ( in ) :: &
      DP_1D
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iVelocity_F, &
      iMomentum_B
    integer ( KDI ), intent ( in ) :: &
      iBaryonMass_F, &
      iBaryonDensity_F, &
      iEnergy_B
    character ( * ), intent ( in ), optional :: &
      SuffixOption

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_N_F_P_HN'

    call S % Slope_DFV_N_Form % Initialize &
           ( RS, DP_1D, iVelocity_F, iMomentum_B, iBaryonMass_F, &
             iBaryonDensity_F, iEnergy_B, SuffixOption )

    !-- Reallocate slope component: Flat_Fluid_Perfect_HeavyNucleus

    deallocate ( S % Component ( 1 ) % Element )
    allocate ( Slope_DFV_F_F_P_HN_Form :: S % Component ( 1 ) % Element )
    select type ( SF  =>  S % Component ( 1 ) % Element )
      class is ( Slope_DFV_F_F_P_HN_Form )

    call SF % Initialize ( RS, DP_1D, SuffixOption )

    end select !-- SF

  end subroutine InitializeAllocate_N


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_N_F_P_HN_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_N_F_P_HN__Form
