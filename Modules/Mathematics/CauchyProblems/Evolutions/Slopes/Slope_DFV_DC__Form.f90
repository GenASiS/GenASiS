module Slope_DFV_DC__Form

  !-- Slope_DivergenceFiniteVolume_DivergenceContribution__Form

  use Basics
  use Fields
  use RiemannSolver_HLL__Form
  use Slope_H__Form
  use Slope_DFV_PD__Form
  use Slope_DFV_C_F__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_DFV_DC_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_DC
    generic, public :: &
      Initialize => InitializeAllocate_DC
    procedure, public, pass :: &
      ComputePartialDerivative
    procedure, public, pass :: &
      ComputeConnectionFlat
    final :: &
      Finalize
  end type Slope_DFV_DC_Form


contains


  subroutine InitializeAllocate_DC &
               ( S, RS, DC, SuffixOption, IgnorabilityOption )

    class ( Slope_DFV_DC_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DivergenceContribution_CS_Form ), intent ( in ) :: &
      DC
    character ( * ), intent ( in ), optional :: &
      SuffixOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_DFV_DC'

    associate ( CS  =>  RS % CurrentSet )

    if ( S % TimerName  ==  '' ) &
      S % TimerName  &
        =  'S_DFV_DC_' // trim ( DC % Name ) // '_' // trim ( CS % Name )

    Name  =  'S_DFV_DC_' // trim ( DC % Name ) // '_' // trim ( CS % Name )
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( Name ) // '_' // trim ( SuffixOption )

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
    call SPD % Initialize ( RS, DC, SuffixOption )
    end select !-- SPD

    nSC  =  nSC + 1
    allocate ( Slope_DFV_C_F_Form :: S % Component ( nSC ) % Element )
    select type ( SCF  =>  S % Component ( nSC ) % Element )
      class is ( Slope_DFV_C_F_Form )
    call SCF % Initialize ( DC, SuffixOption )
    end select !-- SCF

    end associate !-- nSC

  end subroutine InitializeAllocate_DC


  subroutine ComputePartialDerivative ( S, iC, iD, T_Option, iS_Option )

    class ( Slope_DFV_DC_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    select type ( S_PD  =>  S % Component ( 1 ) % Element )
      class is ( Slope_DFV_PD_Form )

    call S_PD % ComputeDimension ( iC, iD, T_Option, iS_Option )

    end select !-- S_PD

  end subroutine ComputePartialDerivative


  subroutine ComputeConnectionFlat ( S, iC, T_Option, iS_Option )

    class ( Slope_DFV_DC_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    select type ( S_C_F  =>  S % Component ( 2 ) % Element )
      class is ( Slope_DFV_C_F_Form )

    call S_C_F % ComputeChart ( iC, T_Option, iS_Option )

    end select !-- S_C_F

  end subroutine ComputeConnectionFlat


  impure elemental subroutine Finalize ( S )

    type ( Slope_DFV_DC_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_DFV_DC__Form
