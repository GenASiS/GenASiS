module Slope_F_RM_I__Form

  !-- Slope_Fluid_RadiationMoments_Interactions__Form

  use Basics
  use Mathematics
  use Fluids
  use RadiationMoments_BM__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_F_RM_I_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_F_RM_I
    generic, public :: &
      Initialize => InitializeAllocate_F_RM_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_F_RM_I_Form


contains


  subroutine InitializeAllocate_F_RM_I ( S, RM, F )

    class ( Slope_F_RM_I_Form ), intent ( inout ) :: &
      S
    class ( RadiationMoments_BM_Form ), intent ( in ) :: &
      RM
    class ( Fluid_P_Form ), intent ( in ) :: &
      F

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_F_RM_I' 
    
    Name  =  trim ( F % Name ) // '_Slp_F_RM_I'

  end subroutine InitializeAllocate_F_RM_I


  subroutine Compute ( S, T_Option )

    class ( Slope_F_RM_I_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      ! associate &
      !   ( SRV  =>   S % Source_R % Storage ( iC ) % Value, &
      !      SV  =>   S % Storage ( iC ) % Value )

      ! SV  =  - SRV

      ! end associate !-- SRV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_F_RM_I_Form ), intent ( inout ) :: &
      S

  end subroutine Finalize


end module Slope_F_RM_I__Form
