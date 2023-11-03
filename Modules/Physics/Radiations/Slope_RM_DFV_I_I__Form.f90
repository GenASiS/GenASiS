module Slope_RM_DFV_I_I__Form

  !-- Slope_RadiationMoments_DivergenceFiniteVolume_Interactions_Implicit__Form

  use Basics
  use Mathematics
  use RadiationMoments_BM__Form
  use Slope_RM_DFV_I__Form

  implicit none
  private

  type, public, extends ( Slope_RM_DFV_I_Form ) :: Slope_RM_DFV_I_I_Form
    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    class ( RadiationMoments_BM_Form ), pointer :: &
      RadiationMoments => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM_DFV_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_RM_DFV_I_I_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel &
               ( S_E, S_S_1, S_S_2, S_S_3, Chi_J, Chi_H, dT, UseDeviceOption ) 
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          S_E, &
          S_S_1, S_S_2, S_S_3
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          Chi_J, Chi_H
        real ( KDR ), intent ( in ) :: &
          dT
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_RM_DFV_I ( S, RS, DF, DT, RM )

    class ( Slope_RM_DFV_I_I_Form ), intent ( inout ) :: &
      S
    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      RS
    class ( DiffusionFactor_CS_Form ), intent ( in ) :: &
      DF
    class ( DivergencePart_CS_Form ), intent ( in ) :: &
      DT
    class ( RadiationMoments_BM_Form ), intent ( in ), target :: &
      RM

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_RM_DFV_I_I' 
    
    Name  =  trim ( RM % Name ) // '_Slp_RM_DFV_I_I'

    S % RadiationMoments  =>  RM

    call Search ( RM % iaBalanced, RM % ENERGY_DENSITY_B, &
                  S % iEnergy_B )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_1, &
                  S % iMomentum_B ( 1 ) )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_2, &
                  S % iMomentum_B ( 2 ) )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_3, &
                  S % iMomentum_B ( 3 ) )

    call S % Slope_RM_DFV_I_Form % Initialize &
           ( RS, DF, DT, RM )

  end subroutine InitializeAllocate_RM_DFV_I


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_RM_DFV_I_I_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    associate &
      ( I  =>  S % RadiationMoments % Interactions )

    call S % Slope_RM_DFV_I_Form % Compute ( dT, T_Option )

    do iC  =  1,  S % Atlas % nCharts

      associate &
        ( IV  =>   I % Storage ( iC ) % Value, &
          SV  =>   S % Storage ( iC ) % Value )
      associate &
        ( Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
          S_E    =>  SV ( :, S % iEnergy_B ), &
          S_S_1  =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2  =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3  =>  SV ( :, S % iMomentum_B ( 3 ) ) )

      call ComputeKernel &
             ( S_E, S_S_1, S_S_2, S_S_3, Chi_J, Chi_H, dT, &
               UseDeviceOption = S % DeviceMemory ) 

      end associate !-- IV, etc.
      end associate !-- Chi_J, etc.

    end do

    end associate !-- I

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_RM_DFV_I_I_Form ), intent ( inout ) :: &
      S

    nullify ( S % RadiationMoments )

  end subroutine Finalize


end module Slope_RM_DFV_I_I__Form
