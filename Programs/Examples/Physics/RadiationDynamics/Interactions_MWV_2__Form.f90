module Interactions_MWV_2__Form

  !-- Interactions_MarshakWaveVaytet_2__Form, Vaytet et al. 2011
  
  use GenASiS
  use Interactions_MWV_1__Form

  implicit none
  private

  type, public, extends ( Interactions_MWV_1_Form ) :: Interactions_MWV_2_Form
    real ( KDR ) :: &
      SpecificOpacityMin = 0.0_KDR, &
      EnergyMax = 0.0_KDR
    class ( PhotonMoments_G_Form ), pointer :: &
      Radiation => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      SetEnergyMax
    procedure, public, pass :: &
      SetRadiation
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Interactions_MWV_2_Form

    real ( KDR ), private, parameter :: &
      PlanckRatio  =  3.83223_KDR  !-- P_4 / P_3

    private :: &
      ComputeKernel

  interface

    module subroutine ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, T, J_Eq, TP, &
               Kappa, E_Max, Ratio_P, k_B, UseDeviceOption )
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
        TP
      real ( KDR ), intent ( in ) :: &
        Kappa, &
        E_Max, &
        Ratio_P, &
        k_B
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeKernel

  end interface


contains


  subroutine InitializeAllocate_I &
               ( I, F, Units_R, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    class ( Fluid_P_Form ), intent ( in ), target :: &
      F
    class ( Units_R_Form ), dimension ( : ), intent ( in ) :: &
      Units_R
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
      I % Type  =  'an Interactions_MWV_2' 
    
    call I % Interactions_BM_Form % Initialize &
           ( F, Units_R, &
             FieldOption = FieldOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetEnergyMax ( I, EnergyMax )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      EnergyMax

    I % EnergyMax  =  EnergyMax

    call Show ( 'Setting EnergyMax of an Interactions_MWV_2', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( I % EnergyMax, 'EnergyMax', &
                I % IGNORABILITY + 1 )

  end subroutine SetEnergyMax


  subroutine SetRadiation ( I, R )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I
    class ( PhotonMoments_G_Form ), intent ( in ), target :: &
      R

    I % Radiation  =>  R

    call Show ( 'Setting Radiation of an Interactions_MWV_2', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( R % Name, 'Radiation', &
                I % IGNORABILITY + 1 )

  end subroutine SetRadiation


  subroutine Compute ( I )

    class ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'Compute', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    associate &
      ( F  =>  I % Fluid, &
        R  =>  I % Radiation )

    do iC  =  1,  I % Atlas % nCharts
      associate &
        ( FV  =>  F % Storage ( iC ) % Value, &
          IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value )
      associate &
        (   M    =>  FV ( :, F % BARYON_MASS ), &
            N    =>  FV ( :, F % BARYON_DENSITY_C ), &
            T    =>  FV ( :, F % TEMPERATURE ), &
           Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
           J_Eq  =>  IV ( :, I % EQUILIBRIUM_J ), &
          TP     =>  RV ( :, R % TEMPERATURE_PARAMETER ) )

      call I % Compute_J_Eq_Ph_G_Kernel &
             ( J_Eq, T, UseDeviceOption = I % DeviceMemory )

      call ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, T, J_Eq, TP, &
               Kappa = I % SpecificOpacity, E_Max = I % EnergyMax, &
               Ratio_P = PlanckRatio, k_B = CONSTANT % BOLTZMANN, &
               UseDeviceOption = I % DeviceMemory )
             
      end associate !-- T, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end associate !-- F

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_2_Form ), intent ( inout ) :: &
      I

    nullify ( I % Radiation )

  end subroutine Finalize


end module Interactions_MWV_2__Form
