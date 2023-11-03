module Slope_F_P_DFV_SS__Form

  !-- Slope_RadiationMoments_DivergenceFiniteVolume_Interactions__Form

  use Basics
  use Mathematics
  use Fluid_P__Form
  use Slope_F_P_SS__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_F_P_DFV_SS_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F_P_DFV_SS
    generic, public :: &
      Initialize => InitializeAllocate_F_P_DFV_SS
    final :: &
      Finalize
  end type Slope_F_P_DFV_SS_Form


contains


  subroutine InitializeAllocate_F_P_DFV_SS ( S, RS, DF, DT, F )

    class ( Slope_F_P_DFV_SS_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_F_P_DFV_SS' 
    
    Name  =  trim ( F % Name ) // '_Slp_F_P_DFV_SS'

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

    !-- Slope component: Divergence

    nSC  =  nSC + 1
    allocate ( Slope_DFV_F_DT_Form :: S % Component ( nSC ) % Element )
    select type ( SD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_F_DT_Form )

    call SD % Initialize ( RS, DF, DT )

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

  end subroutine InitializeAllocate_F_P_DFV_SS


  impure elemental subroutine Finalize ( S )

    type ( Slope_F_P_DFV_SS_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_F_P_DFV_SS__Form
