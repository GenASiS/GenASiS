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
    final :: &
      Finalize
  end type Interactions_NM_G_Form

    private :: &
      Compute_EA_E_Kernel, &
      Compute_EA_E_Bar_Kernel

    interface

      module subroutine Compute_EA_E_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_eq, N_eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 UseDeviceOption )
        !-- Compute_EmissionAbsorption_Electron_Kernel
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_eq, N_eq, F_Ave, &
          M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_E_Kernel

      module subroutine Compute_EA_E_Bar_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_eq, N_eq, F_Ave, M, N, T, X_n, Mu_e, &
                 UseDeviceOption )
        !-- Compute_EmissionAbsorption_Electron_Bar_Kernel
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
           Xi_J,  Xi_H,  Xi_N, &
          Chi_J, Chi_H, Chi_N
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          J_eq, N_eq, F_Ave, &
          M, N, T, X_n, Mu_e
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_EA_E_Bar_Kernel

      module subroutine Compute_S_N_A_Kernel &
               ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, &
                 UseDeviceOption )
        !-- Compute_Scattering_Nucleons_Nuclei
        real ( KDR ), dimension ( : ), intent ( inout ) :: &
          Chi_H
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          T_nu, Eta_nu, &
          M, N, X_p, X_n, X_A, Z, A
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_S_N_A_Kernel

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
    
    Name  =  'Interactions'
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
        call Compute_EA_E_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_eq, N_eq, F_Ave, M, N, T, X_p, X_A, Z, A, Mu_e, Mu_n_p, &
                 UseDeviceOption = I % DeviceMemory )
      case ( 'NEUTRINOS_E_BAR' )
        call Compute_EA_E_Bar_Kernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
                 J_eq, N_eq, F_Ave, M, N, T, X_n, Mu_e, &
                 UseDeviceOption = I % DeviceMemory )
      end select !-- RadiationType
             
      !-- Elastic scattering on nucleons and nuclei

      call Compute_S_N_A_Kernel &
             ( Chi_H, T_nu, Eta_nu, M, N, X_p, X_n, X_A, Z, A, &
               UseDeviceOption = I % DeviceMemory )

      end associate !-- Xi_J, etc.
      end associate !-- FV, etc.
    end do !-- iC

    end select !-- F
    end select !-- R

  end subroutine ComputeAll


  impure elemental subroutine Finalize ( I )

    type ( Interactions_NM_G_Form ), intent ( inout ) :: &
      I

  end subroutine Finalize


end module Interactions_NM_G__Form
