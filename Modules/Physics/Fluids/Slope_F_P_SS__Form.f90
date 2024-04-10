module Slope_F_P_SS__Form

  !-- Slope_Fluid_Perfect_SplitSource__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_F_P_SS_Form
    class ( Fluid_P_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F_P_SS
    generic, public :: &
      Initialize => InitializeAllocate_F_P_SS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_F_P_SS_Form


contains


  subroutine InitializeAllocate_F_P_SS ( S, F )

    class ( Slope_F_P_SS_Form ), intent ( inout ) :: &
      S
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_F_P_SS' 
    
    Name  =  trim ( F % Name ) // '_Slp_F_P_SS'

    S % Fluid  =>  F

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY )

  end subroutine InitializeAllocate_F_P_SS


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_F_P_SS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate ( FSS  =>  S % Fluid % SplitSource )
    call FSS % Copy ( S )
    end associate !-- FS

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_F_P_SS_Form ), intent ( inout ) :: &
      S

    nullify ( S % Fluid )

  end subroutine Finalize


end module Slope_F_P_SS__Form
