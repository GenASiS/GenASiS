module Slope_RM_DFV_I__Form

  !-- Slope_RadiationMoments_Streaming_Interactions__Form

  use Basics
  use Mathematics
  use RadiationMoments_BM__Form
  use Slope_RM_I__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_RM_DFV_I_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM_DFV_I
    generic, public :: &
      Initialize => InitializeAllocate_RM_DFV_I
    final :: &
      Finalize
  end type Slope_RM_DFV_I_Form


contains


  subroutine InitializeAllocate_RM_DFV_I ( S, RS, DF, DT, RM )

    class ( Slope_RM_DFV_I_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    class ( RadiationMoments_BM_Form ), intent ( in ) :: &
      RM

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_RM_DFV_I' 
    
    Name  =  trim ( RM % Name ) // '_Slp_RM_DFV_I'

    call S % Slope_H_Form % Initialize &
           ( RM % Atlas, &
             FieldOption = RM % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = RM % DeviceMemory, &
             PinnedMemoryOption = RM % PinnedMemory, &
             DevicesCommunicateOption = RM % DevicesCommunicate, &
             nFieldsOption = RM % nBalanced, &
             IgnorabilityOption = RM % IGNORABILITY + 1 )

    associate ( nSC  =>  S % nComponents )

    !-- Slope component: Divergence

    nSC  =  nSC + 1
    allocate ( Slope_DFV_F_DT_Form :: S % Component ( nSC ) % Element )
    select type ( SD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_F_DT_Form )

    call SD % Initialize ( RS, DF, DT )

    end select !-- SD

    !-- Cleanup

    end associate !-- nSC

  end subroutine InitializeAllocate_RM_DFV_I


  impure elemental subroutine Finalize ( S )

    type ( Slope_RM_DFV_I_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_RM_DFV_I__Form
