module DivergencePart_F_P_T__Form

  !-- DivergencePart_Fluid_Perfect_Total__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

  type, public, extends ( DivergencePart_CS_Form ) :: DivergencePart_F_P_T_Form
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass ( DP ) :: &
      ComputeFluxes
    procedure, public, pass ( DP ) :: &
      ComputeStresses
    final :: &
      Finalize
  end type DivergencePart_F_P_T_Form

    private :: &
      Compute_FS_G_Kernel, &
      Compute_S_UD_Kernel

    interface

      module subroutine Compute_FS_G_Kernel &
               ( D, S_1, S_2, S_3, G, P, V_Dim, iDim, &
                 F_D, F_S_1, F_S_2, F_S_3, F_G, UseDeviceOption )
        !-- Compute_FluxSet_Galileo_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          D, &
          S_1, S_2, S_3, &
          G, &
          P, &
          V_Dim
        integer ( KDI ), intent ( in ) :: &
          iDim
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          F_D, &
          F_S_1, F_S_2, F_S_3, &
          F_G
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_FS_G_Kernel

      module subroutine Compute_S_UD_Kernel &
               ( V_2, V_3, S_2, S_3, P, S_UD_22, S_UD_33, UseDeviceOption )
        !-- Compute_Stress_UD_Kernel
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          V_2, V_3, &
          S_2, S_3, &
          P
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          S_UD_22, S_UD_33
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_S_UD_Kernel

    end interface


contains


  subroutine Initialize ( DP, CS, NameOption, IgnorabilityOption )

    class ( DivergencePart_F_P_T_Form ), intent ( inout ) :: &
      DP
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption

    character ( LDL ) :: &
      Name

    if ( DP % Type  ==  '' ) &
      DP % Type  =  'a DivergencePart_F_P_T' 

    Name  =  'V'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call DP % DivergencePart_CS_Form % Initialize &
           ( CS, NameOption = Name, &
             IgnorabilityOption = IgnorabilityOption )

  end subroutine Initialize


  subroutine ComputeFluxes ( FS_F, DP, FS_CS, iC, iD )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      FS_F  !-- Fluxes
    class ( DivergencePart_F_P_T_Form ), intent ( in ) :: &
      DP
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      FS_CS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimension

    integer ( KDI ) :: &
      iDensity, &
      iEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum

    select type ( CS  =>  DP % CurrentSet )
      class is ( Fluid_P_Form )

    call Search &
           ( CS % iaBalanced, CS % BARYON_DENSITY_B, iDensity )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_1, iMomentum ( 1 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_2, iMomentum ( 2 ) )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_3, iMomentum ( 3 ) )
    call Search &
           ( CS % iaBalanced, CS % ENERGY_DENSITY_B, iEnergy )

    associate &
      ( FSV  =>  FS_F  % Storage ( iC ) % Value, &
        CSV  =>  FS_CS % Storage ( iC ) % Value )
    associate &
      ( F_D      =>  FSV ( :, iDensity ), &
        F_S_1    =>  FSV ( :, iMomentum ( 1 ) ), &
        F_S_2    =>  FSV ( :, iMomentum ( 2 ) ), &
        F_S_3    =>  FSV ( :, iMomentum ( 3 ) ), &
        F_G      =>  FSV ( :, iEnergy ), &
          D      =>  CSV ( :, CS % BARYON_DENSITY_B ), &
          S_1    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_1 ), &
          S_2    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
          S_3    =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
          G      =>  CSV ( :, CS % ENERGY_DENSITY_B ), &
          P      =>  CSV ( :, CS % PRESSURE ), &
          V_Dim  =>  CSV ( :, CS % VELOCITY_U ( iD ) ) )
 
    call Compute_FS_G_Kernel &
           ( D, S_1, S_2, S_3, G, P, V_Dim, iD, F_D, F_S_1, F_S_2, F_S_3, F_G, &
             UseDeviceOption = CS % DeviceMemory )
  
    end associate !-- F_D, etc.
    end associate !-- FSV, etc.
    end select !-- CS

  end subroutine ComputeFluxes


  subroutine ComputeStresses ( S_UD, DP, iC, iMomentum_1, iMomentum_2 )

    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      S_UD
    class ( DivergencePart_F_P_T_Form ), intent ( in ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iChart
    integer ( KDI ), intent ( out ) :: &
      iMomentum_1, iMomentum_2

    select type ( CS  =>  DP % CurrentSet )
      class is ( Fluid_P_Form )

    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_1, iMomentum_1 )
    call Search &
           ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D_2, iMomentum_2 )

    associate &
      (  S_UD_V  =>   S_UD % Storage ( iC ) % Value, &
        CSV      =>  CS    % Storage ( iC ) % Value )
    associate &
      ( S_UD_22  =>  S_UD_V ( :, 1 ), &
        S_UD_33  =>  S_UD_V ( :, 2 ), &
           V_2   =>  CSV ( :, CS % VELOCITY_U_2 ), &
           V_3   =>  CSV ( :, CS % VELOCITY_U_3 ), &
           S_2   =>  CSV ( :, CS % MOMENTUM_DENSITY_D_2 ), &
           S_3   =>  CSV ( :, CS % MOMENTUM_DENSITY_D_3 ), &
           P     =>  CSV ( :, CS % PRESSURE ) )
 
    call Compute_S_UD_Kernel &
           ( V_2, V_3, S_2, S_3, P, S_UD_22, S_UD_33, &
             UseDeviceOption = CS % DeviceMemory )

    end associate !-- S_UD_22, etc.
    end associate !-- S_UD_V, etc.
    end select !-- CS

  end subroutine ComputeStresses


  impure elemental subroutine Finalize ( DP )

    type ( DivergencePart_F_P_T_Form ), intent ( inout ) :: &
      DP

  end subroutine Finalize


end module DivergencePart_F_P_T__Form
