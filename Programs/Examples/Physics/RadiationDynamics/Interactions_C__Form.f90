module Interactions_C__Form

  !-- Interactions_Constant_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Interactions_BM_Form ) :: Interactions_C_Form
    real ( KDR ) :: &
      OpacityAbsorption = 0.0_KDR
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass :: &
      SetOpacityAbsorption
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Interactions_C_Form

    private :: &
      ComputeKernel

  interface

    module subroutine ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, J_Eq, Kappa_A, UseDeviceOption )
      use Basics
      implicit none
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
         Xi_J, &
        Chi_J, &
        Chi_H
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        J_Eq
      real ( KDR ), intent ( in ) :: &
        Kappa_A
      logical ( KDL ), intent ( in ), optional :: &
        UseDeviceOption
    end subroutine ComputeKernel

  end interface

contains


  subroutine InitializeAllocate_I &
               ( I, R, Units_R, F, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_C_Form ), intent ( inout ) :: &
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
      I % Type  =  'an Interactions_C' 
    
    call I % Interactions_BM_Form % Initialize &
           ( R, Units_R, F, &
             FieldOption = FieldOption, &
             NameOption = NameOption, &
             UnitOption = UnitOption, &
             nFieldsOption = nFieldsOption, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine InitializeAllocate_I


  subroutine SetOpacityAbsorption ( I, Kappa_A )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I
    real ( KDR ), intent ( in ) :: &
      Kappa_A

    I % OpacityAbsorption  =  Kappa_A

    call Show ( 'Setting OpacityAbsorption of an Interactions_C', &
                I % IGNORABILITY + 1 )
    call Show ( I % Name, 'Name', &
                I % IGNORABILITY + 1 )
    call Show ( I % OpacityAbsorption, 'OpacityAbsorption', &
                I % IGNORABILITY + 1 )

  end subroutine SetOpacityAbsorption


  subroutine Compute ( I )

    class ( Interactions_C_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'Compute', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    associate ( F  =>  I % Fluid )

    do iC  =  1,  I % Atlas % nCharts
      associate &
        ( FV  =>  F % Storage ( iC ) % Value, &
          IV  =>  I % Storage ( iC ) % Value )
      associate &
        (   T    =>  FV ( :, F % TEMPERATURE ), &
           Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
           J_Eq  =>  IV ( :, I % EQUILIBRIUM_J ) )

      call I % Compute_J_Eq_Ph_G_Kernel &
             ( J_Eq, T, UseDeviceOption = I % DeviceMemory )

      call ComputeKernel &
             ( Xi_J, Chi_J, Chi_H, J_Eq, Kappa_A = I % OpacityAbsorption, &
               UseDeviceOption = I % DeviceMemory )
             
      end associate !-- T, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end associate !-- F

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_C_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_C__Form
