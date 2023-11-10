module Interactions_MWV_1__Form

  !-- Interactions_MarshakWaveVaytet_1__Form, Vaytet et al. 2011
  
  use GenASiS

  implicit none
  private

  type, public, extends ( Interactions_BM_Form ) :: Interactions_MWV_1_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      SetSpecificOpacity
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Interactions_MWV_1_Form

    private :: &
      ComputeKernel

  interface

    module subroutine ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, J_Eq, Kappa, UseDeviceOption )
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
         Xi_J, &
        Chi_J, &
        Chi_H
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M, &
        N, &
        J_Eq
      real ( KDR ), intent ( in ) :: &
        Kappa
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeKernel

  end interface


contains


  subroutine InitializeAllocate_I &
               ( I, R, Units_R, F, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
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
      I % Type  =  'an Interactions_MWV_1' 
    
    call I % Interactions_BM_Form % Initialize &
           ( R, Units_R, F, &
             FieldOption = FieldOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetSpecificOpacity ( I, SpecificOpacity )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      SpecificOpacity

    I % SpecificOpacity  =  SpecificOpacity

    call Show ( 'Setting SpecificOpacity of an Interactions_MWV_1', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( I % SpecificOpacity, 'SpecificOpacity', &
                I % IGNORABILITY + 1 )

  end subroutine SetSpecificOpacity


  subroutine Compute ( I )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'Compute', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    associate &
      ( R  =>  I % Radiation, &
        F  =>  I % Fluid )

    do iC  =  1,  I % Atlas % nCharts
      associate &
        ( IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value, &
          FV  =>  F % Storage ( iC ) % Value )
      associate &
        (   M     =>  FV ( :, F % BARYON_MASS ), &
            N     =>  FV ( :, F % BARYON_DENSITY_C ), &
           Xi_J   =>  IV ( :, I % EMISSIVITY_J ), &
          Chi_J   =>  IV ( :, I % OPACITY_J ), &
          Chi_H   =>  IV ( :, I % OPACITY_H ), &
            J_Eq  =>  RV ( :, R % ENERGY_DENSITY_C_EQ ) )

      call ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, J_Eq, Kappa = I % SpecificOpacity, &
               UseDeviceOption = I % DeviceMemory )
             
      end associate !-- T, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end associate !-- R, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_MWV_1__Form
