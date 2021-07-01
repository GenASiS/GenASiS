module Slope_DFV_F__Form

  !-- Slope_DivergenceFiniteVolume_Flat__Form

  use Basics
  use RiemannSolver_HLL__Form
  use Slope_H__Form
  use Slope_DFV_PD__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_F_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Slope_DFV_F_Form


contains


  subroutine InitializeAllocate_F ( S, RS, SuffixOption )

    class ( Slope_DFV_F_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ), target :: &
      RS
    character ( * ), intent ( in ), optional :: &
      SuffixOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_F'

    Name  =  'Slp_DFV_F_' // trim ( RS % CurrentSet % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

    associate ( CS  =>  RS % CurrentSet )

    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

    end associate !-- CS

    !-- Slope component: Partial derivative

    associate ( nSC  =>  S % nComponents )
    nSC  =  nSC + 1
    allocate ( Slope_DFV_PD_Form :: S % Component ( nSC ) % Element )
    select type ( SPD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_PD_Form )

    call SPD % Initialize ( RS, SuffixOption )

    end select !-- SPD
    end associate !-- nSC

  end subroutine InitializeAllocate_F


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_F_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_F__Form
