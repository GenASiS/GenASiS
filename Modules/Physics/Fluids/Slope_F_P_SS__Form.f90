module Slope_F_P_S__Form

  !-- Slope_Fluid_Perfect_Source__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_F_P_S_Form
    class ( Fluid_P_Form ), pointer :: &
      Fluid => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F_P_S
    generic, public :: &
      Initialize => InitializeAllocate_F_P_S
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_F_P_S_Form


contains


  subroutine InitializeAllocate_F_P_S ( S, F )

    class ( Slope_F_P_S_Form ), intent ( inout ) :: &
      S
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_F_P_S' 
    
    Name  =  trim ( F % Name ) // '_Slp_F_P_S'

    S % Fluid  =>  F

    call S % Slope_H_Form % Initialize &
           ( F % Atlas, &
             FieldOption = F % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = F % DeviceMemory, &
             PinnedMemoryOption = F % PinnedMemory, &
             DevicesCommunicateOption = F % DevicesCommunicate, &
             nFieldsOption = F % nBalanced, &
             IgnorabilityOption = F % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_F_P_S


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_F_P_S_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( FSV  =>   S % Fluid % Source % Storage ( iC ) % Value, &
           SV  =>   S % Storage ( iC ) % Value )

      SV  =  FSV

      end associate !-- FSV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_F_P_S_Form ), intent ( inout ) :: &
      S

    nullify ( S % Fluid )

  end subroutine Finalize


end module Slope_F_P_S__Form
