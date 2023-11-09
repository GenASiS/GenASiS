module Slope_RM_I__Form

  !-- Slope_RadiationMoments_Interactions__Form

  use Basics
  use Mathematics
  use Interactions_BM__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_RM_I_Form
    integer ( KDI ) :: &
      iEnergy_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    class ( Interactions_BM_Form ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_RM_I
    generic, public :: &
      Initialize => InitializeAllocate_RM_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_RM_I_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel &
               ( ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, E, S_1, S_2, S_3, &
                 dT, S_E, S_S_1, S_S_2, S_S_3, UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
           Xi_J,  Xi_H, &
          Chi_J, Chi_H, &
          E, &
          S_1, S_2, S_3
        real ( KDR ), intent ( in ) :: &
          dT
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_E, &
          S_S_1, S_S_2, S_S_3
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_RM_I ( S, I )

    class ( Slope_RM_I_Form ), intent ( inout ) :: &
      S
    class ( Interactions_BM_Form ), intent ( in ), target :: &
      I

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_RM_I' 
    
    associate ( RM  =>  I % Radiation )

    Name  =  trim ( RM % Name ) // '_Slp_RM_I'

    S % Interactions  =>  I

    call Search ( RM % iaBalanced, RM % ENERGY_DENSITY_B, &
                  S % iEnergy_B )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_1, &
                  S % iMomentum_B ( 1 ) )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_2, &
                  S % iMomentum_B ( 2 ) )
    call Search ( RM % iaBalanced, RM % MOMENTUM_DENSITY_B_D_3, &
                  S % iMomentum_B ( 3 ) )

    call S % Slope_H_Form % Initialize &
           ( RM % Atlas, &
             FieldOption = RM % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = RM % DeviceMemory, &
             PinnedMemoryOption = RM % PinnedMemory, &
             DevicesCommunicateOption = RM % DevicesCommunicate, &
             nFieldsOption = RM % nBalanced, &
             IgnorabilityOption = RM % IGNORABILITY + 1 )

    end associate !-- RM

  end subroutine InitializeAllocate_RM_I


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_RM_I_Form ), intent ( inout ) :: &
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
      (  I  =>  S % Interactions, &
        RM  =>  S % Interactions % Radiation )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( IV  =>   I % Storage ( iC ) % Value, &
          RV  =>  RM % Storage ( iC ) % Value, &
          SV  =>   S % Storage ( iC ) % Value )
      associate &
        (  Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H  =>  IV ( :, I % EMISSIVITY_H ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
            E    =>  RV ( :, RM % ENERGY_DENSITY_B ), &
            S_1  =>  RV ( :, RM % MOMENTUM_DENSITY_B_D_1 ), &
            S_2  =>  RV ( :, RM % MOMENTUM_DENSITY_B_D_2 ), &
            S_3  =>  RV ( :, RM % MOMENTUM_DENSITY_B_D_3 ), &
          S_E    =>  SV ( :, S % iEnergy_B ), &
          S_S_1  =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2  =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3  =>  SV ( :, S % iMomentum_B ( 3 ) ) )

      call ComputeKernel &
             ( C % ProperCell, Xi_J, Xi_H, Chi_J, Chi_H, E, S_1, S_2, S_3, &
               dT, S_E, S_S_1, S_S_2, S_S_3, &
               UseDeviceOption = S % DeviceMemory )

      end associate !-- Xi_J, etc.
      end associate !-- IV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_RM_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- I, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_RM_I_Form ), intent ( inout ) :: &
      S

    nullify ( S % Interactions )
    
  end subroutine Finalize


end module Slope_RM_I__Form
