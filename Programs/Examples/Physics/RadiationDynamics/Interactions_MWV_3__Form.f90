module Interactions_MWV_3__Form

  !-- Interactions_MarshakWaveVaytet_3__Form, Vaytet et al. 2011
  
  use GenASiS
  use Interactions_MWV_2__Form

  implicit none
  private

  type, public, extends ( Interactions_MWV_2_Form ) :: Interactions_MWV_3_Form
    real ( KDR ) :: &
      TemperatureScale = 1.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      SetTemperatureScale
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Interactions_MWV_3_Form

    real ( KDR ), private, parameter :: &
      PlanckRatio  =  3.83223_KDR  !-- P_4 / P_3

    private :: &
      ComputeKernel

  interface

    module subroutine ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, T, J_Eq, T_R, &
               Kappa, E_Max, T_0, Ratio_P, k_B, UseDeviceOption )
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
         Xi_J, &
        Chi_J, &
        Chi_H
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        T, &
        J_Eq, &
        T_R
      real ( KDR ), intent ( in ) :: &
        Kappa, &
        E_Max, &
        T_0, &
        Ratio_P, &
        k_B
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeKernel

  end interface


contains


  subroutine InitializeAllocate_I &
               ( I, R, Units_R, F, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_MWV_3_Form ), intent ( inout ) :: &
      I
    class ( RadiationMoments_BM_Form ), intent ( inout ), target :: &
      R
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( :, : ), intent ( in ), optional :: &
      UnitOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption, &
      IgnorabilityOption

    if ( I % Type  ==  '' ) &
      I % Type  =  'an Interactions_MWV_3' 
    
    call I % Interactions_BM_Form % Initialize &
           ( R, Units_R, F, &
             FieldOption = FieldOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetTemperatureScale ( I, TemperatureScale )

    class ( Interactions_MWV_3_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      TemperatureScale

    I % TemperatureScale  =  TemperatureScale

    call Show ( 'Setting TemperatureScale of an Interactions_MWV_3', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( I % EnergyMax, 'TemperatureScale', &
                I % IGNORABILITY + 1 )

  end subroutine SetTemperatureScale


  subroutine Compute ( I )

    class ( Interactions_MWV_3_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'Compute', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    select type ( R  =>  I % Radiation )
      class is ( PhotonMoments_G_Form )
    associate &
      ( F  =>  I % Fluid )

    do iC  =  1,  I % Atlas % nCharts
      associate &
        ( IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value, &
          FV  =>  F % Storage ( iC ) % Value )
      associate &
        (   M     =>  FV ( :, F % BARYON_MASS ), &
            N     =>  FV ( :, F % BARYON_DENSITY_C ), &
            T     =>  FV ( :, F % TEMPERATURE ), &
           Xi_J   =>  IV ( :, I % EMISSIVITY_J ), &
          Chi_J   =>  IV ( :, I % OPACITY_J ), &
          Chi_H   =>  IV ( :, I % OPACITY_H ), &
            J_Eq  =>  RV ( :, R % ENERGY_DENSITY_C_EQ ), &
            T_R   =>  RV ( :, R % TEMPERATURE_GREY ) )

      call ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, T, J_Eq, T_R, &
               Kappa = I % SpecificOpacity, E_Max = I % EnergyMax, &
               T_0 = I % TemperatureScale, Ratio_P = PlanckRatio, &
               k_B = CONSTANT % BOLTZMANN, UseDeviceOption = I % DeviceMemory )
             
      end associate !-- T, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end associate !-- F
    end select !-- R

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_3_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )

  end subroutine Finalize


end module Interactions_MWV_3__Form
