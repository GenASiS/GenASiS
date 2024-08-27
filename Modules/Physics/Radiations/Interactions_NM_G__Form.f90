module Interactions_NM_G__Form

  use Basics
  use Mathematics
  use Fluids
  use Units_R__Form
  use RadiationMoments_BM__Form
  use Interactions_BM__Form
  use NeutrinoMoments_G__Form

 implicit none
 private

     integer ( KDI ), private, parameter :: &
      N_FIELDS_NM_G = 10

  type, public, extends ( Interactions_BM_Form ) :: Interactions_NM_G_Form
    integer ( KDI ) :: &
      N_FIELDS_NM_G = N_FIELDS_NM_G
    integer ( KDI ) :: &
      EMISSIVITY_N = 0, &
         OPACITY_N = 0
    integer ( KDI ) :: &
      !-- Emission/absorption, nucleons
      EMISSIVITY_J_EA_N = 0, & 
         OPACITY_J_EA_N = 0, &
      !-- Emission/absorption, nuclei
      EMISSIVITY_J_EA_A = 0, & 
         OPACITY_J_EA_A = 0, &
      !-- Pairs, electron/positron
      EMISSIVITY_J_P_EP = 0, & 
         OPACITY_J_P_EP = 0, & 
      !-- Scattering on nucleons
      OPACITY_H_S_N = 0, &
      !-- Scattering on nuclei
      OPACITY_H_S_A = 0
    real ( KDR ) :: &
      DensityDetailedBalance
  contains
    procedure, private, pass :: &
      InitializeAllocate_I
    procedure, public, pass ( I ) :: &
      SetStream
    procedure, private, pass :: &
      ComputeAll
    procedure, private, pass :: &
      ComputeSingle
    final :: &
      Finalize
  end type Interactions_NM_G_Form
  
  public :: &
    Compute_EA_E_S_Kernel, &
    Compute_EA_EB_S_Kernel, &
    Compute_EA_HL_S_Kernel, &
    Compute_P_S_Kernel, & 
    Compute_S_S_Kernel 

    private :: &
      Compute_EA_E_A_Kernel, &
      Compute_EA_EB_A_Kernel, &
      Compute_EA_HL_A_Kernel, &

      Compute_P_A_Kernel, & 
      Compute_S_A_Kernel

    interface

      module subroutine Compute_EA_E_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 Rho_DB, UseDeviceOption )
        !-- Compute_EmissionAbsorption_Electron_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
          M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_E_A_Kernel

      module subroutine Compute_EA_E_S_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 Rho_DB, iV )
        !-- Compute_EmissionAbsorption_Electron_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
          M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_EA_E_S_Kernel

      module subroutine Compute_EA_EB_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, Mu_e, Mu_n_p, &
                 Rho_DB, UseDeviceOption )
        !-- Compute_EmissionAbsorption_ElectronBar_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
          M, N, T, X_n, X_p, Mu_e, Mu_n_p
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_EB_A_Kernel

      module subroutine Compute_EA_EB_S_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, Mu_e, Mu_n_p, &
                 Rho_DB, iV )
        !-- Compute_EmissionAbsorption_ElectronBar_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
          M, N, T, X_n, X_p, Mu_e, Mu_n_p
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_EA_EB_S_Kernel

      module subroutine Compute_EA_HL_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 Chi_H_S_N, Chi_H_S_A, &
                 UseDeviceOption )
        !-- Compute_EmissionAbsorption_HeavyLepton_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A, &
          Chi_H_S_N,  Chi_H_S_A
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_HL_A_Kernel

      module subroutine Compute_EA_HL_S_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 Chi_H_S_N, Chi_H_S_A, &
                 iV )
        !-- Compute_EmissionAbsorption_HeavyLepton_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_EA_N,  Xi_J_EA_A, &
          Chi_J_EA_N, Chi_J_EA_A, &
          Chi_H_S_N,  Chi_H_S_A
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_EA_HL_S_Kernel

      module subroutine Compute_P_A_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_P_EP, Chi_J_P_EP, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign, nSpecies, Rho_DB, UseDeviceOption )
        !-- Compute_Pair_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Xi_J, Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_P_EP, &
          Chi_J_P_EP
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, M, N, T, Mu_e
        integer ( KDI ), intent ( in ) :: &
          Sign, &
          nSpecies
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_P_A_Kernel

      module subroutine Compute_P_S_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_P_EP, Chi_J_P_EP, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign, nSpecies, Rho_DB, iV )
        !-- Compute_Pair_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Xi_J, Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J_P_EP, &
          Chi_J_P_EP
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, M, N, T, Mu_e
        real ( KDR ), intent ( in ) :: &
          Rho_DB
        integer ( KDI ), intent ( in ) :: &
          Sign, &
          nSpecies, &
          iV
      end subroutine Compute_P_S_Kernel

      module subroutine Compute_S_A_Kernel &
               ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
                 M, N, X_p, X_n, X_A, Z, A, UseDeviceOption )
        !-- Compute_Scattering_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H_S_N, Chi_H_S_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T_nu, Eta_nu, &
          M, N, X_p, X_n, X_A, Z, A
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_S_A_Kernel

      module subroutine Compute_S_S_Kernel &
               ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
                 M, N, X_p, X_n, X_A, Z, A, iV )
        !-- Compute_Scattering_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H_S_N, Chi_H_S_A
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T_nu, Eta_nu, &
          M, N, X_p, X_n, X_A, Z, A
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_S_S_Kernel

    end interface


contains


  subroutine InitializeAllocate_I &
               ( I, R, Units_R, F, FieldOption, NameOption, UnitOption, &
                 nFieldsOption, IgnorabilityOption )

    class ( Interactions_NM_G_Form ), intent ( inout ) :: &
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

    integer ( KDI ) :: &
      iC, &  !-- iChart
      oF, &  !-- oField
      nFields
    type ( QuantityForm ), dimension ( :, : ), allocatable :: &
      FieldUnit
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( I % Type  ==  '' ) &
      I % Type  =  'an Interactions_NM_G' 
    
    Name  =  'Interactions_' // trim ( R % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    !-- Field indices

    oF  =  I % N_FIELDS_I

    I % EMISSIVITY_N  =  oF + 1
    I %    OPACITY_N  =  oF + 2

    I % EMISSIVITY_J_EA_N  =  oF + 3 
    I %    OPACITY_J_EA_N  =  oF + 4

    I % EMISSIVITY_J_EA_A  =  oF + 5 
    I %    OPACITY_J_EA_A  =  oF + 6

    I % EMISSIVITY_J_P_EP  =  oF + 7 
    I %    OPACITY_J_P_EP  =  oF + 8

    I % OPACITY_H_S_N  =  oF +  9
    I % OPACITY_H_S_A  =  oF + 10

    nFields  =  oF  +  I % N_FIELDS_NM_G
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + I % N_FIELDS_NM_G ) &
      = [ 'Emissivity_N     ', &
          'Opacity_N        ', &
          'Emissivity_J_EA_N', &
          'Opacity_J_EA_N   ', &
          'Emissivity_J_EA_A', &
          'Opacity_J_EA_A   ', &
          'Emissivity_J_P_EP', &
          'Opacity_J_P_EP   ', &
          'Opacity_H_S_N    ', &
          'Opacity_H_S_A    ' ]
          
    !-- Units

    associate ( nC  =>  F % Atlas % nCharts )

    if ( present ( UnitOption ) ) then
      allocate ( FieldUnit, source = UnitOption )
    else
      allocate ( FieldUnit ( nFields, nC ) )
    end if !-- FieldOption

    do iC  =  1, nC

      FieldUnit ( I % EMISSIVITY_N, iC ) &
        =  Units_R ( iC ) % NumberDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_N, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)

      FieldUnit ( I % EMISSIVITY_J_EA_N, iC ) &
        =  Units_R ( iC ) % EnergyDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_J_EA_N, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)

      FieldUnit ( I % EMISSIVITY_J_EA_A, iC ) &
        =  Units_R ( iC ) % EnergyDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_J_EA_A, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)

      FieldUnit ( I % EMISSIVITY_J_P_EP, iC ) &
        =  Units_R ( iC ) % EnergyDensity  &
           *  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_J_P_EP, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)

      FieldUnit ( I % OPACITY_H_S_N, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)
      FieldUnit ( I % OPACITY_H_S_A, iC ) &
        =  ( UNIT % SPEED_OF_LIGHT * Units_R ( iC ) % Time ) ** (-1)

    end do !-- iC

    end associate !-- nC

    !-- Interactions_BM

    call I % Interactions_BM_Form % Initialize &
           ( R, Units_R, F, &
             FieldOption = Field, &
             NameOption = Name, &
             UnitOption = FieldUnit, &
             nFieldsOption = nFields, &
             IgnorabilityOption = IgnorabilityOption )

    call R % SetInteractions ( I )

    !-- Parameters

    I % DensityDetailedBalance  =  1.0e13_KDR * UNIT % MASS_DENSITY_CGS

  end subroutine InitializeAllocate_I


  subroutine SetStream ( S, I )

    class ( Stream_BM_Form ), intent ( inout ) :: &
      S
    class ( Interactions_NM_G_Form ), intent ( in ) :: &
      I

    call S % AddFieldSet &
           ( I, &
             iaSelectedOption &
               = [ I % EMISSIVITY_J, &
                   I % EMISSIVITY_H, &
                   I % EMISSIVITY_N, &
                   I % OPACITY_J, &
                   I % OPACITY_H, &
                   I % OPACITY_N, &
                   I % EMISSIVITY_J_EA_N, &
                   I % OPACITY_J_EA_N, &
                   I % EMISSIVITY_J_EA_A, &
                   I % OPACITY_J_EA_A, &
                   I % EMISSIVITY_J_P_EP, &
                   I % OPACITY_J_P_EP, &
                   I % OPACITY_H_S_N, &
                   I % OPACITY_H_S_A ] )

  end subroutine SetStream


  subroutine ComputeAll ( I )

    class ( Interactions_NM_G_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC

    call Show ( 'ComputeAll', CONSOLE % INFO_6 )
    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    select type ( R  =>  I % Radiation )
      class is ( NeutrinoMoments_G_Form )
    select type ( F  =>  I % Fluid )
      class is ( Fluid_P_HN_Form )

    call R % ComputeSpectralParameters ( )
    call R % ComputeEquilibrium ( )

    do iC  =  1,  I % Atlas % nCharts
      associate &
        ( IV  =>  I % Storage ( iC ) % Value, &
          RV  =>  R % Storage ( iC ) % Value, &
          FV  =>  F % Storage ( iC ) % Value )
      associate &
        (  Xi_J          =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H          =>  IV ( :, I % EMISSIVITY_H ), &
           Xi_N          =>  IV ( :, I % EMISSIVITY_N ), &
           Xi_J_EA_N     =>  IV ( :, I % EMISSIVITY_J_EA_N ), &
           Xi_J_EA_A     =>  IV ( :, I % EMISSIVITY_J_EA_A ), &
           Xi_J_P_EP     =>  IV ( :, I % EMISSIVITY_J_P_EP ), &
          Chi_J          =>  IV ( :, I % OPACITY_J ), &
          Chi_H          =>  IV ( :, I % OPACITY_H ), &
          Chi_N          =>  IV ( :, I % OPACITY_N ), &
          Chi_J_EA_N     =>  IV ( :, I % OPACITY_J_EA_N ), &
          Chi_J_EA_A     =>  IV ( :, I % OPACITY_J_EA_A ), &
          Chi_J_P_EP     =>  IV ( :, I % OPACITY_J_P_EP ), &
          Chi_H_S_N      =>  IV ( :, I % OPACITY_H_S_N ), &
          Chi_H_S_A      =>  IV ( :, I % OPACITY_H_S_A ), &
            T_Nu         =>  RV ( :, R % TEMPERATURE_GREY ), &
          Eta_Nu         =>  RV ( :, R % DEGENERACY_GREY ), &
            J_Eq         =>  RV ( :, R % ENERGY_DENSITY_C_EQ ), &
            N_Eq         =>  RV ( :, R % NUMBER_DENSITY_C_EQ ), &
            E_Ave        =>  RV ( :, R % ENERGY_AVERAGE ), &
            F_Ave        =>  RV ( :, R % OCCUPANCY_AVERAGE ), &
            M            =>  FV ( :, F % BARYON_MASS ), &
            N            =>  FV ( :, F % BARYON_DENSITY_C ), &
            T            =>  FV ( :, F % TEMPERATURE ), &
            X_p          =>  FV ( :, F % MASS_FRACTION_PROTON ), &
            X_n          =>  FV ( :, F % MASS_FRACTION_NEUTRON ), &
            X_A          =>  FV ( :, F % MASS_FRACTION_HEAVY ), &
            Z            =>  FV ( :, F % ATOMIC_NUMBER_HEAVY ), &
            A            =>  FV ( :, F % MASS_NUMBER_HEAVY ), &
           Mu_e          =>  FV ( :, F % CHEMICAL_POTENTIAL_E ), &
           Mu_n_p        =>  FV ( :, F % CHEMICAL_POTENTIAL_N_P ) )

      !-- Emission / Absorption

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )
        call Compute_EA_E_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 Rho_DB = I % DensityDetailedBalance, &
                 UseDeviceOption = I % DeviceMemory )
      case ( 'NEUTRINOS_EB' )
        call Compute_EA_EB_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
                 M, N, T, X_n, X_p, Mu_e, Mu_n_p, &
                 Rho_DB = I % DensityDetailedBalance, &
                 UseDeviceOption = I % DeviceMemory )
!      case ( 'NEUTRINOS_EB', 'NEUTRINOS_HL' )
      case ( 'NEUTRINOS_HL' )
        call Compute_EA_HL_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
                 Chi_H_S_N, Chi_H_S_A, &
                 UseDeviceOption = I % DeviceMemory )
      end select !-- RadiationType
             
      !-- Pair emission

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
        call Compute_P_A_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 Xi_J_P_EP, Chi_J_P_EP, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign = +1, nSpecies = 1, &
                 Rho_DB = I % DensityDetailedBalance, &
                 UseDeviceOption = I % DeviceMemory )
      ! case ( 'NEUTRINOS_HL' )
      !   call Compute_P_A_Kernel &
      !          ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
      !            Xi_J_P_EP, Chi_J_P_EP, &
      !            J_Eq, N_Eq, M, N, T, Mu_e, &
      !            Sign = -1, nSpecies = 4, &
      !            Rho_DB = I % DensityDetailedBalance, &
      !            UseDeviceOption = I % DeviceMemory )
      end select !-- RadiationType

      !-- Elastic scattering on nucleons and nuclei

      select case ( trim ( R % RadiationType ) )
!      case ( 'NEUTRINOS_E' )
      case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
!      case ( 'NEUTRINOS_E', 'NEUTRINOS_EB', 'NEUTRINOS_HL' )
        call Compute_S_A_Kernel &
               ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
                 M, N, X_p, X_n, X_A, Z, A, &
                 UseDeviceOption = I % DeviceMemory )
      end select !-- RadiationType

      end associate !-- Xi_J, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end select !-- F
    end select !-- R

  end subroutine ComputeAll


  subroutine ComputeSingle ( I, iC, iV )

    class ( Interactions_NM_G_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC, &
      iV

integer ( KDI ) :: &
  iV_S

!iV_S = 88
!if ( iV == iV_S ) then
!  call Show ( '>>> Interactions_NM_G % ComputeSingle', CONSOLE % INFO_6 )
!  call Show ( I % Name, '>>> Interactions', CONSOLE % INFO_6 )
!end if

    select type ( R  =>  I % Radiation )
      class is ( NeutrinoMoments_G_Form )
    select type ( F  =>  I % Fluid )
      class is ( Fluid_P_HN_Form )

    associate &
      ( I_V  =>  I % Storage ( iC ) % Value, &
        R_V  =>  R % Storage ( iC ) % Value, &
        F_V  =>  F % Storage ( iC ) % Value )
    associate &
      (  Xi_J          =>  I_V ( :, I % EMISSIVITY_J ), &
         Xi_H          =>  I_V ( :, I % EMISSIVITY_H ), &
         Xi_N          =>  I_V ( :, I % EMISSIVITY_N ), &
         Xi_J_EA_N     =>  I_V ( :, I % EMISSIVITY_J_EA_N ), &
         Xi_J_EA_A     =>  I_V ( :, I % EMISSIVITY_J_EA_A ), &
         Xi_J_P_EP     =>  I_V ( :, I % EMISSIVITY_J_P_EP ), &
        Chi_J          =>  I_V ( :, I % OPACITY_J ), &
        Chi_H          =>  I_V ( :, I % OPACITY_H ), &
        Chi_N          =>  I_V ( :, I % OPACITY_N ), &
        Chi_J_EA_N     =>  I_V ( :, I % OPACITY_J_EA_N ), &
        Chi_J_EA_A     =>  I_V ( :, I % OPACITY_J_EA_A ), &
        Chi_J_P_EP     =>  I_V ( :, I % OPACITY_J_P_EP ), &
        Chi_H_S_N      =>  I_V ( :, I % OPACITY_H_S_N ), &
        Chi_H_S_A      =>  I_V ( :, I % OPACITY_H_S_A ), &
          T_Nu         =>  R_V ( :, R % TEMPERATURE_GREY ), &
        Eta_Nu         =>  R_V ( :, R % DEGENERACY_GREY ), &
          J_Eq         =>  R_V ( :, R % ENERGY_DENSITY_C_EQ ), &
          N_Eq         =>  R_V ( :, R % NUMBER_DENSITY_C_EQ ), &
          E_Ave        =>  R_V ( :, R % ENERGY_AVERAGE ), &
          F_Ave        =>  R_V ( :, R % OCCUPANCY_AVERAGE ), &
          M            =>  F_V ( :, F % BARYON_MASS ), &
          N            =>  F_V ( :, F % BARYON_DENSITY_C ), &
          T            =>  F_V ( :, F % TEMPERATURE ), &
          X_p          =>  F_V ( :, F % MASS_FRACTION_PROTON ), &
          X_n          =>  F_V ( :, F % MASS_FRACTION_NEUTRON ), &
          X_A          =>  F_V ( :, F % MASS_FRACTION_HEAVY ), &
          Z            =>  F_V ( :, F % ATOMIC_NUMBER_HEAVY ), &
          A            =>  F_V ( :, F % MASS_NUMBER_HEAVY ), &
         Mu_e          =>  F_V ( :, F % CHEMICAL_POTENTIAL_E ), &
         Mu_n_p        =>  F_V ( :, F % CHEMICAL_POTENTIAL_N_P ) )

! if ( iV == iV_S ) then
!   call Show ( T_Nu ( iV_S ), '>>> T_Nu before ComputeSpectralParameters' )
!   call Show ( Eta_Nu ( iV_S ), '>>> Eta_Nu before ComputeSpectralParameters' )
! end if

    call R % ComputeSpectralParameters ( iC, iV )
    call R % ComputeEquilibrium ( iC, iV )

! if ( iV == iV_S ) then
!   call Show ( T_Nu ( iV_S ), '>>> T_Nu after ComputeSpectralParameters' )
!   call Show ( Eta_Nu ( iV_S ), '>>> Eta_Nu after ComputeSpectralParameters' )
! end if

    !-- Emission / Absorption

    select case ( trim ( R % RadiationType ) )
    case ( 'NEUTRINOS_E' )
      call Compute_EA_E_S_Kernel &
             ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
               J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
               M, N, T, X_n, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
               Rho_DB = I % DensityDetailedBalance, iV = iV )
    case ( 'NEUTRINOS_EB' )
      call Compute_EA_EB_S_Kernel &
             ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
               J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
               M, N, T, X_n, X_p, Mu_e, Mu_n_p, &
               Rho_DB = I % DensityDetailedBalance, iV = iV )
!      case ( 'NEUTRINOS_EB', 'NEUTRINOS_HL' )
    case ( 'NEUTRINOS_HL' )
      call Compute_EA_HL_S_Kernel &
             ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
               Chi_H_S_N, Chi_H_S_A, &
               iV = iV )
    end select !-- RadiationType
           
    !-- Pair emission

    select case ( trim ( R % RadiationType ) )
    case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
      call Compute_P_S_Kernel &
             ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
               Xi_J_P_EP, Chi_J_P_EP, &
               J_Eq, N_Eq, M, N, T, Mu_e, &
               Sign = +1, nSpecies = 1, &
               Rho_DB = I % DensityDetailedBalance, &
               iV = iV )
    ! case ( 'NEUTRINOS_HL' )
    !   call Compute_P_S_Kernel &
    !          ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
    !            Xi_J_P_EP, Chi_J_P_EP, &
    !            J_Eq, N_Eq, M, N, T, Mu_e, &
    !            Sign = -1, nSpecies = 4, &
    !            Rho_DB = I % DensityDetailedBalance, &
    !            iV = iV )
    end select !-- RadiationType

    !-- Elastic scattering on nucleons and nuclei

!if ( iV == iV_S ) then
!  call Show ( R % RadiationType, '>>> RadiationType' )
!  call Show ( Chi_H ( iV_S ), '>>> Chi_H before scattering' )
!  call Show ( Chi_H_S_N ( iV_S ), '>>> Chi_H_S_N' )
!  call Show ( Chi_H_S_A ( iV_S ), '>>> Chi_H_S_A' )
!end if

    select case ( trim ( R % RadiationType ) )
!    case ( 'NEUTRINOS_E' )
    case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
!    case ( 'NEUTRINOS_E', 'NEUTRINOS_EB', 'NEUTRINOS_HL' )
      call Compute_S_S_Kernel &
             ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
               M, N, X_p, X_n, X_A, Z, A, iV )
   end select !-- RadiationType

!if ( iV == iV_S ) then
!  call Show ( Chi_H ( iV_S ), '>>> Chi_H after scattering' )
!  call Show ( Chi_H_S_N ( iV_S ), '>>> Chi_H_S_N' )
!  call Show ( Chi_H_S_A ( iV_S ), '>>> Chi_H_S_A' )
!end if

    end associate !-- Xi_J, etc.
    end associate !-- FV, etc.

    end select !-- F
    end select !-- R

  end subroutine ComputeSingle


  impure elemental subroutine Finalize ( I )

    type ( Interactions_NM_G_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_NM_G__Form
