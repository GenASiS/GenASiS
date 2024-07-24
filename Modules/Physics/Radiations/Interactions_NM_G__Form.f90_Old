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
      N_FIELDS_NM_G = 2

  type, public, extends ( Interactions_BM_Form ) :: Interactions_NM_G_Form
    integer ( KDI ) :: &
      N_FIELDS_NM_G = N_FIELDS_NM_G
    integer ( KDI ) :: &
      EMISSIVITY_N = 0, &
      OPACITY_N    = 0
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

    private :: &
      Compute_EA_E_A_Kernel, &
      Compute_EA_E_S_Kernel, &
      Compute_EA_E_Bar_A_Kernel, &
      Compute_EA_E_Bar_S_Kernel, &
      Compute_S_N_A_A_Kernel, & 
      Compute_S_N_A_S_Kernel 

    interface

      module subroutine Compute_EA_E_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 UseDeviceOption )
        !-- Compute_EmissionAbsorption_Electron_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, F_Ave, &
          M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_E_A_Kernel

      module subroutine Compute_EA_E_S_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, iV )
        !-- Compute_EmissionAbsorption_Electron_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, F_Ave, &
          M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_EA_E_S_Kernel

      module subroutine Compute_EA_E_Bar_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_n, Mu_e, &
                 UseDeviceOption )
        !-- Compute_EmissionAbsorption_Electron_Bar_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, F_Ave, &
          M, N, T, X_n, Mu_e
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_E_Bar_A_Kernel

      module subroutine Compute_EA_E_Bar_S_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_n, Mu_e, iV )
        !-- Compute_EmissionAbsorption_Electron_Bar_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, F_Ave, &
          M, N, T, X_n, Mu_e
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_EA_E_Bar_S_Kernel

      module subroutine Compute_P_A_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign, nSpecies, UseDeviceOption )
        !-- Compute_Pair_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Xi_J, Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, M, N, T, Mu_e
        integer ( KDI ), intent ( in ) :: &
          Sign, &
          nSpecies
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_P_A_Kernel

      module subroutine Compute_P_S_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign, nSpecies, iV )
        !-- Compute_Pair_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Xi_J, Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_Eq, N_Eq, M, N, T, Mu_e
        integer ( KDI ), intent ( in ) :: &
          Sign, &
          nSpecies, &
          iV
      end subroutine Compute_P_S_Kernel

      module subroutine Compute_S_N_A_A_Kernel &
               ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, &
                 UseDeviceOption )
        !-- Compute_Scattering_Nucleons_Nuclei_All_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T_nu, Eta_nu, &
          M, N, X_p, X_n, X_A, Z, A
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_S_N_A_A_Kernel

      module subroutine Compute_S_N_A_S_Kernel &
               ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, iV )
        !-- Compute_Scattering_Nucleons_Nuclei_Single_Kernel
        implicit none
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T_nu, Eta_nu, &
          M, N, X_p, X_n, X_A, Z, A
        integer ( KDI ), intent ( in ) :: &
          iV
      end subroutine Compute_S_N_A_S_Kernel

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
    I % OPACITY_N     =  oF + 2

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
      = [ 'Emissivity_N', &
          'Opacity_N   ' ]
          
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
                   I % OPACITY_N ] )

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
        (  Xi_J    =>  IV ( :, I % EMISSIVITY_J ), &
           Xi_H    =>  IV ( :, I % EMISSIVITY_H ), &
           Xi_N    =>  IV ( :, I % EMISSIVITY_N ), &
          Chi_J    =>  IV ( :, I % OPACITY_J ), &
          Chi_H    =>  IV ( :, I % OPACITY_H ), &
          Chi_N    =>  IV ( :, I % OPACITY_N ), &
            T_Nu   =>  RV ( :, R % TEMPERATURE_GREY ), &
          Eta_Nu   =>  RV ( :, R % DEGENERACY_GREY ), &
            J_Eq   =>  RV ( :, R % ENERGY_DENSITY_C_EQ ), &
            N_Eq   =>  RV ( :, R % NUMBER_DENSITY_C_EQ ), &
            F_Ave  =>  RV ( :, R % OCCUPANCY_AVERAGE ), &
            M      =>  FV ( :, F % BARYON_MASS ), &
            N      =>  FV ( :, F % BARYON_DENSITY_C ), &
            T      =>  FV ( :, F % TEMPERATURE ), &
            X_p    =>  FV ( :, F % MASS_FRACTION_PROTON ), &
            X_n    =>  FV ( :, F % MASS_FRACTION_NEUTRON ), &
            X_A    =>  FV ( :, F % MASS_FRACTION_HEAVY ), &
            Z      =>  FV ( :, F % ATOMIC_NUMBER_HEAVY ), &
            A      =>  FV ( :, F % MASS_NUMBER_HEAVY ), &
           Mu_e    =>  FV ( :, F % CHEMICAL_POTENTIAL_E ), &
           Mu_n_p  =>  FV ( :, F % CHEMICAL_POTENTIAL_N_P ) )

      !-- Emission / Absorption

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )
        call Compute_EA_E_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 UseDeviceOption = I % DeviceMemory )
      case ( 'NEUTRINOS_EB' )
        call Compute_EA_E_Bar_A_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, F_Ave, M, N, T, X_n, Mu_e, &
                 UseDeviceOption = I % DeviceMemory )
      case ( 'NEUTRINOS_HL' )
        call I % Clear ( )
      end select !-- RadiationType
             
      !-- Pair emission

      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
        call Compute_P_A_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign = +1, nSpecies = 1, UseDeviceOption = I % DeviceMemory )
      case ( 'NEUTRINOS_HL' )
        call Compute_P_A_Kernel &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_Eq, N_Eq, M, N, T, Mu_e, &
                 Sign = -1, nSpecies = 4, UseDeviceOption = I % DeviceMemory )
      end select !-- RadiationType

      !-- Elastic scattering on nucleons and nuclei

      call Compute_S_N_A_A_Kernel &
             ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, &
               UseDeviceOption = I % DeviceMemory )

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
  iShow

!    call Show ( 'ComputeSingle', CONSOLE % INFO_6 )
!    call Show ( I % Name, 'Interactions', CONSOLE % INFO_6 )

    select type ( R  =>  I % Radiation )
      class is ( NeutrinoMoments_G_Form )
    select type ( F  =>  I % Fluid )
      class is ( Fluid_P_HN_Form )

    associate &
      ( I_V  =>  I % Storage ( iC ) % Value, &
        R_V  =>  R % Storage ( iC ) % Value, &
        F_V  =>  F % Storage ( iC ) % Value )
    associate &
      (  Xi_J    =>  I_V ( :, I % EMISSIVITY_J ), &
         Xi_H    =>  I_V ( :, I % EMISSIVITY_H ), &
         Xi_N    =>  I_V ( :, I % EMISSIVITY_N ), &
        Chi_J    =>  I_V ( :, I % OPACITY_J ), &
        Chi_H    =>  I_V ( :, I % OPACITY_H ), &
        Chi_N    =>  I_V ( :, I % OPACITY_N ), &
          T_Nu   =>  R_V ( :, R % TEMPERATURE_GREY ), &
        Eta_Nu   =>  R_V ( :, R % DEGENERACY_GREY ), &
          J_Eq   =>  R_V ( :, R % ENERGY_DENSITY_C_EQ ), &
          N_Eq   =>  R_V ( :, R % NUMBER_DENSITY_C_EQ ), &
          F_Ave  =>  R_V ( :, R % OCCUPANCY_AVERAGE ), &
          M      =>  F_V ( :, F % BARYON_MASS ), &
          N      =>  F_V ( :, F % BARYON_DENSITY_C ), &
          T      =>  F_V ( :, F % TEMPERATURE ), &
          X_p    =>  F_V ( :, F % MASS_FRACTION_PROTON ), &
          X_n    =>  F_V ( :, F % MASS_FRACTION_NEUTRON ), &
          X_A    =>  F_V ( :, F % MASS_FRACTION_HEAVY ), &
          Z      =>  F_V ( :, F % ATOMIC_NUMBER_HEAVY ), &
          A      =>  F_V ( :, F % MASS_NUMBER_HEAVY ), &
         Mu_e    =>  F_V ( :, F % CHEMICAL_POTENTIAL_E ), &
         Mu_n_p  =>  F_V ( :, F % CHEMICAL_POTENTIAL_N_P ) )

! iShow = 31
! if ( iV == iShow ) then
!   call Show ( Eta_Nu ( iV ), '>>> Eta_Nu' )
! end if

! call Show ( R % RadiationType, '>>> Type' )

! call Show ( Eta_Nu ( iV ), '>>> Eta_Nu entry' ) 
! call Show ( T_Nu ( iV ), '>>> T_Nu entry' )

    call R % ComputeSpectralParameters ( iC, iV )
    call R % ComputeEquilibrium ( iC, iV )

! associate &
!   ( E_Nu  =>  R_V ( :, R % ENERGY_DENSITY_C ), &
!     N_Nu  =>  R_V ( :, R % NUMBER_DENSITY_C ) )
! call Show ( E_Nu ( iV ), '>>> E_Nu' )
! call Show ( N_Nu ( iV ), '>>> N_Nu' )
! call Show ( Eta_Nu ( iV ), '>>> Eta_Nu' ) 
! call Show ( T_Nu ( iV ), '>>> T_Nu' )
! end associate !-- E_Nu, etc.

    !-- Emission / Absorption

    select case ( trim ( R % RadiationType ) )
    case ( 'NEUTRINOS_E' )
      call Compute_EA_E_S_Kernel &
             ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               J_Eq, N_Eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, iV )
    case ( 'NEUTRINOS_EB' )
      call Compute_EA_E_Bar_S_Kernel &
             ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
               J_Eq, N_Eq, F_Ave, M, N, T, X_n, Mu_e, iV )
    case ( 'NEUTRINOS_HL' )
      call I % Clear ( )
    end select !-- RadiationType
! call Show ( Xi_J ( iV ), '>>> Xi_J after EA' )
! call Show ( Xi_N ( iV ), '>>> Xi_N after EA' )
! call Show ( Chi_J ( iV ), '>>> Chi_J after EA' )
! call Show ( Chi_N ( iV ), '>>> Chi_N after EA' )
! call Show ( Chi_H ( iV ), '>>> Chi_H after EA' )
           
    !-- Pair emission

    select case ( trim ( R % RadiationType ) )
    case ( 'NEUTRINOS_E', 'NEUTRINOS_EB' )
      call Compute_P_S_Kernel &
             ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
               J_Eq, N_Eq, M, N, T, Mu_e, &
               Sign = +1, nSpecies = 1, iV = iV )
    case ( 'NEUTRINOS_HL' )
      call Compute_P_S_Kernel &
             ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
               J_Eq, N_Eq, M, N, T, Mu_e, &
               Sign = -1, nSpecies = 4, iV = iV )
! call Show ( Xi_J ( iV ), '>>> Xi_J after P' )
! call Show ( Xi_N ( iV ), '>>> Xi_N after P' )
! call Show ( Chi_J ( iV ), '>>> Chi_J after P' )
! call Show ( Chi_N ( iV ), '>>> Chi_N after P' )
! call Show ( Chi_H ( iV ), '>>> Chi_H after P' )
    end select !-- RadiationType

    !-- Elastic scattering on nucleons and nuclei

    call Compute_S_N_A_S_Kernel &
           ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, iV )
! call Show ( Chi_H ( iV ), '>>> Chi_H after ISO' )

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
