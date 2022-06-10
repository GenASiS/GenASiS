module Slope_DFV_F_DT__Form

  !-- Slope_DivergenceFiniteVolume_Flat_DivergenceTotal__Form

  use Basics
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form
  use Slope_DFV_PD__Form
  use Slope_DFV_C_F__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_F_DT_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    final :: &
      Finalize
  end type Slope_DFV_F_DT_Form


contains


  subroutine InitializeAllocate_F ( S, RS, DT, IgnorabilityOption )

    class ( Slope_DFV_F_DT_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_F_DT'

    Name  =  trim ( RS % CurrentSet % Name ) // '_Slp_DFV_F_DT'

    associate ( CS  =>  RS % CurrentSet )
    call S % Slope_H_Form % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             PinnedMemoryOption = CS % PinnedMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = IgnorabilityOption )
    end associate !-- CS

    !-- Slope components: Partial derivative and flat connection

    associate ( nSC  =>  S % nComponents )

    nSC  =  nSC + 1
    allocate ( Slope_DFV_PD_Form :: S % Component ( nSC ) % Element )
    select type ( SPD  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_PD_Form )
    call SPD % Initialize ( RS, DT )
    end select !-- SPD

    nSC  =  nSC + 1
    allocate ( Slope_DFV_C_F_Form :: S % Component ( nSC ) % Element )
    select type ( SCF  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_C_F_Form )
    call SCF % Initialize ( DT )
    end select !-- SCF

    S % StreamComponents  =  .false.

    end associate !-- nSC

  end subroutine InitializeAllocate_F


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_F_DT_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_F_DT__Form
