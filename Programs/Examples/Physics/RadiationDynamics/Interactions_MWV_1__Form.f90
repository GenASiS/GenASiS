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
    procedure, private, pass :: &
      ComputeAll
    procedure, private, pass :: &
      ComputeSingle
    final :: &
      Finalize
  end type Interactions_MWV_1_Form

    private :: &
      ComputeAllKernel

  interface

    module subroutine ComputeAllKernel &
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
    end subroutine ComputeAllKernel

    module subroutine ComputeSingleKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, J_Eq, Kappa, iV, UseDeviceOption )
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
      integer ( KDI ), intent ( in ) :: &
        iV
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeSingleKernel

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


  subroutine ComputeAll ( I )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeAll', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    select type ( R  =>  I % Radiation )
      class is ( PhotonMoments_G_Form )
    associate &
      ( F  =>  I % Fluid )

    call R % ComputeSpectralParameters ( )
    call R % ComputeEquilibrium ( )

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

      call ComputeAllKernel &
             ( Xi_J, Chi_J, Chi_H, M, N, J_Eq, Kappa = I % SpecificOpacity, &
               UseDeviceOption = I % DeviceMemory )
             
      end associate !-- T, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end associate !-- F
    end select !-- R

  end subroutine ComputeAll


  subroutine ComputeSingle ( I, iC, iV )

    class ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

!    call Show ( 'ComputeAll', CONSOLE % INFO_6 )
!    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    select type ( R  =>  I % Radiation )
      class is ( PhotonMoments_G_Form )
    associate &
      ( F  =>  I % Fluid )

    associate &
      ( I_V  =>  I % Storage ( iC ) % Value, &
        R_V  =>  R % Storage ( iC ) % Value, &
        F_V  =>  F % Storage ( iC ) % Value )
    associate &
      (   M     =>  F_V ( :, F % BARYON_MASS ), &
          N     =>  F_V ( :, F % BARYON_DENSITY_C ), &
         Xi_J   =>  I_V ( :, I % EMISSIVITY_J ), &
        Chi_J   =>  I_V ( :, I % OPACITY_J ), &
        Chi_H   =>  I_V ( :, I % OPACITY_H ), &
          J_Eq  =>  R_V ( :, R % ENERGY_DENSITY_C_EQ ) )

    call R % ComputeSpectralParameters ( iC, iV )
    call R % ComputeEquilibrium ( iC, iV )

    call ComputeSingleKernel &
           ( Xi_J, Chi_J, Chi_H, M, N, J_Eq, I % SpecificOpacity, iV, &
             UseDeviceOption = I % DeviceMemory )
           
    end associate !-- T, etc.
    end associate !-- FV, etc.

    end associate !-- F
    end select !-- R

  end subroutine ComputeSingle


  impure elemental subroutine Finalize ( I )

    type ( Interactions_MWV_1_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_MWV_1__Form
