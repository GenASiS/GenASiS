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
    class ( NeutrinoMoments_G_Form ), pointer :: &
      Radiation => null ( )
    class ( Interactions_NM_G_Form ), pointer :: &
      Interactions => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_NM_G_I
    generic, public :: &
      Initialize => InitializeAllocate_NM_G_I
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


  subroutine InitializeAllocate_NM_G_I ( S, R )

    class ( Slope_NM_G_I_Form ), intent ( inout ) :: &
      S
    class ( NeutrinoMoments_G_Form ), intent ( in ), target :: &
      R

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Slope_NM_G_I' 
    
    Name  =  trim ( R % Name ) // '_Slp_NM_I'

    S % Radiation  =>  R

    !-- FIXME: As a workaround because the correct type of R % Interactions 
    !          is not being recognized, this pointer is previously assigned 
    ! select type ( I  =>  R % Interactions )
    ! class is ( Interactions_NM_G_Form )
    !   S % Interactions  =>  I
    ! class default
    !   call Show ( 'Interactions type not recognized', CONSOLE % ERROR )
    !   call Show ( 'Slope_NM_G_I__Form', 'module', CONSOLE % ERROR )
    !   call Show ( 'InitializeAllocate_NM_G_I', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end select

    call Search ( R % iaBalanced, R % ENERGY_DENSITY_B, &
                  S % iEnergy_B )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, &
                  S % iMomentum_B ( 1 ) )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, &
                  S % iMomentum_B ( 2 ) )
    call Search ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, &
                  S % iMomentum_B ( 3 ) )
    call Search ( R % iaBalanced, R % NUMBER_DENSITY_B, &
                  S % iNumber_B )

    call S % Slope_H_Form % Initialize &
           ( R % Atlas, &
             FieldOption = R % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = R % DeviceMemory, &
             PinnedMemoryOption = R % PinnedMemory, &
             DevicesCommunicateOption = R % DevicesCommunicate, &
             nFieldsOption = R % nBalanced, &
             IgnorabilityOption = R % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_NM_G_I


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

    associate &
      (  I  =>  S % Interactions, &
         R  =>  S % Radiation )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      associate &
        ( IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value, &
          SV  =>  S % Storage ( iC ) % Value )
      associate &
        (  Xi_J  =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H  =>  IV ( :, I % EMISSIVITY_H ), &
           Xi_N  =>  IV ( :, I % EMISSIVITY_N ), &
          Chi_J  =>  IV ( :, I % OPACITY_J ), &
          Chi_H  =>  IV ( :, I % OPACITY_H ), &
          Chi_N  =>  IV ( :, I % OPACITY_N ), &
            E    =>  RV ( :, R % ENERGY_DENSITY_B ), &
            S_1  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_1 ), &
            S_2  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_2 ), &
            S_3  =>  RV ( :, R % MOMENTUM_DENSITY_B_D_3 ), &
            D    =>  RV ( :, R % NUMBER_DENSITY_B ), &
          S_E    =>  SV ( :, S % iEnergy_B ), &
          S_S_1  =>  SV ( :, S % iMomentum_B ( 1 ) ), &
          S_S_2  =>  SV ( :, S % iMomentum_B ( 2 ) ), &
          S_S_3  =>  SV ( :, S % iMomentum_B ( 3 ) ), &
          S_D    =>  SV ( :, S % iNumber_B ) )

!call Show ( Xi_J, '>>> Xi_J' )
      call ComputeKernel &
             ( C % ProperCell, Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               E, S_1, S_2, S_3, D, dT, S_E, S_S_1, S_S_2, S_S_3, S_D, &
               UseDeviceOption = S % DeviceMemory )
!call Show ( S_E, '>>> S_E' )

      end associate !-- Xi_J, etc.
      end associate !-- IV, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Slope_NM_G_I__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- I, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( S )

    type ( Slope_NM_G_I_Form ), intent ( inout ) :: &
      S

    nullify ( S % Interactions )
    
  end subroutine Finalize


end module Slope_NM_G_I__Form
