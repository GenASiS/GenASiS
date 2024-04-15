module Slope_NM_G_DFV_I__Form

  !-- Slope_RadiationMoments_DivergenceFiniteVolume_Interactions__Form

  use Basics
  use Mathematics
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form
  use Slope_NM_G_I__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_NM_G_DFV_I_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_NM_G_DFV_I
    generic, public :: &
      Initialize => InitializeAllocate_NM_G_DFV_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_NM_G_DFV_I_Form


contains


  subroutine InitializeAllocate_NM_G_DFV_I ( S, RS, DF, DT, R )

    class ( Slope_NM_G_DFV_I_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_NM_G_DFV_I' 
    
    Name  =  trim ( R % Name ) // '_Slp_NM_G_DFV_I'

    call S % Slope_H_Form % Initialize &
           ( R % Atlas, &
             FieldOption = R % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = R % DeviceMemory, &
             PinnedMemoryOption = R % PinnedMemory, &
             DevicesCommunicateOption = R % DevicesCommunicate, &
             nFieldsOption = R % nBalanced, &
             IgnorabilityOption = R % IGNORABILITY + 1 )

    associate ( nSC  =>  S % nComponents )

    !-- Slope component: Divergence

    nSC  =  nSC + 1
    allocate ( Slope_DFV_F_DT_Form :: S % Component ( nSC ) % Element )
    select type ( SD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_F_DT_Form )

    call SD % Initialize ( RS, DF, DT )

    !-- Slope component: Interactions

    nSC  =  nSC + 1
    allocate ( Slope_NM_G_I_Form :: S % Component ( nSC ) % Element )
    select type ( SI  =>  S % Component ( nSC ) % Element )
      class is ( Slope_NM_G_I_Form )

    call SI % Initialize ( R )

    end select !-- SD
    end select !-- SI

    !-- Cleanup

    end associate !-- nSC

  end subroutine InitializeAllocate_NM_G_DFV_I


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_NM_G_DFV_I_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( S_I  =>  S % Component ( 2 ) % Element )
      class is ( Slope_NM_G_I_Form )
    associate &
      ( R  =>  S_I % Radiation )

!    call R % SetFluidVelocity ( )

    end associate !-- R, F
    end select !-- S_I

    call S % Slope_H_Form % Compute ( dT, T_Option )

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_NM_G_DFV_I_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_NM_G_DFV_I__Form
