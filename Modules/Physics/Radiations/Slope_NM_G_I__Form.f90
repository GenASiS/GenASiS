module Slope_NM_G_I__Form

  !-- Slope_NeutrinoMoments_Interactions__Form

  use Basics
  use Mathematics
  use Interactions_BM__Form
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form

  implicit none
  private

  type, public, extends ( Slope_H_Form ) :: Slope_NM_G_I_Form
    integer ( KDI ) :: &
      iEnergy_B, &
      iNumber_B
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_B
    class ( Interactions_BM_Form ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_NM_I
    generic, public :: &
      Initialize => InitializeAllocate_NM_I
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Slope_NM_G_I_Form

    private :: &
      ComputeKernel

    interface

      module subroutine ComputeKernel &
               ( ProperCell, Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 E, S_1, S_2, S_3, D, dT, S_E, S_S_1, S_S_2, S_S_3, S_D, &
                 UseDeviceOption )
        use Basics
        implicit none
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N, &
          E, &
          S_1, S_2, S_3, &
          D
        real ( KDR ), intent ( in ) :: &
          dT
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_E, &
          S_S_1, S_S_2, S_S_3, &
          S_D
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_NM_I ( S, I )

    class ( Slope_NM_G_I_Form ), intent ( inout ) :: &
      S
    class ( Interactions_BM_Form ), intent ( in ), target :: &
      I

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_NM_G_I' 
    
    select type ( NM  =>  I % Radiation )
      class is ( NeutrinoMoments_G_Form )

    Name  =  trim ( NM % Name ) // '_Slp_NM_I'

    S % Interactions  =>  I

    call Search ( NM % iaBalanced, NM % ENERGY_DENSITY_B, &
                  S % iEnergy_B )
    call Search ( NM % iaBalanced, NM % MOMENTUM_DENSITY_B_D_1, &
                  S % iMomentum_B ( 1 ) )
    call Search ( NM % iaBalanced, NM % MOMENTUM_DENSITY_B_D_2, &
                  S % iMomentum_B ( 2 ) )
    call Search ( NM % iaBalanced, NM % MOMENTUM_DENSITY_B_D_3, &
                  S % iMomentum_B ( 3 ) )
    call Search ( NM % iaBalanced, NM % NUMBER_DENSITY_B, &
                  S % iNumber_B )

    call S % Slope_H_Form % Initialize &
           ( NM % Atlas, &
             FieldOption = NM % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = NM % DeviceMemory, &
             PinnedMemoryOption = NM % PinnedMemory, &
             DevicesCommunicateOption = NM % DevicesCommunicate, &
             nFieldsOption = NM % nBalanced, &
             IgnorabilityOption = NM % IGNORABILITY + 1 )

    end select !-- NM

  end subroutine InitializeAllocate_NM_I


  subroutine Compute ( S, dT, T_Option )

    class ( Slope_NM_G_I_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    select type ( I  =>  S % Interactions )
      class is ( Interactions_NM_G_Form )
    select type ( NM  =>  I % Radiation )
      class is ( NeutrinoMoments_G_Form )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( IV  =>   I % Storage ( iC ) % Value, &
          RV  =>  NM % Storage ( iC ) % Value, &
          SV  =>   S % Storage ( iC ) % Value )
      associate &
        (  Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H  =>  IV ( :, I % EMISSIVITY_H ), &
           Xi_N  =>  IV ( :, I % EMISSIVITY_N ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
          Chi_N  =>  IV ( :, I % OPACITY_N ), &
            E    =>  RV ( :, NM % ENERGY_DENSITY_B ), &
            S_1  =>  RV ( :, NM % MOMENTUM_DENSITY_B_D_1 ), &
            S_2  =>  RV ( :, NM % MOMENTUM_DENSITY_B_D_2 ), &
            S_3  =>  RV ( :, NM % MOMENTUM_DENSITY_B_D_3 ), &
            D    =>  RV ( :, NM % NUMBER_DENSITY_B ), &
          S_E    =>  SV ( :, S % iEnergy_B ), &
          S_S_1  =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2  =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3  =>  SV ( :, S % iMomentum_B ( 3 ) ), &
          S_D    =>  SV ( :, S % iNumber_B ) )

      call ComputeKernel &
             ( C % ProperCell, Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               E, S_1, S_2, S_3, D, dT, S_E, S_S_1, S_S_2, S_S_3, S_D, &
               UseDeviceOption = S % DeviceMemory )

      end associate !-- Xi_J, etc.
      end associate !-- IV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_NM_G_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end select !-- NM
    end select !-- I

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_NM_G_I_Form ), intent ( inout ) :: &
      S

    nullify ( S % Interactions )
    
  end subroutine Finalize


end module Slope_NM_G_I__Form
