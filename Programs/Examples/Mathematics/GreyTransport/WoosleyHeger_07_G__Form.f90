module WoosleyHeger_07_G__Form

  !-- WoosleyHeger_07_Grey__Form

  use Mathematics
  use Fluid_ASC__Form
  use WoosleyHeger_07_Header_Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use Interactions_Template
  use Interactions_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ) :: &
    WoosleyHeger_07_G_Form
      integer ( KDI ) :: &
        NEUTRINOS_E_NU             = 1, &
        NEUTRINOS_E_NU_BAR         = 2, &
        NEUTRINOS_MU_TAU_NU_NU_BAR = 3, &
        FLUID                      = 4         
      type ( WoosleyHeger_07_HeaderForm ), allocatable :: &
        Header
      type ( Interactions_ASC_Form ), allocatable :: &
        Interactions_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type WoosleyHeger_07_G_Form

    private :: &
      ApplySources_Radiation, &
      ApplySources_Fluid, &
      ApplySourcesKernel

    type ( VariableGroupForm ), private, pointer :: &
      Increment_E_Nu           => null ( ), &
      Increment_E_NuBar        => null ( ), &
      Increment_MuTau_Nu_NuBar => null ( )
    class ( InteractionsTemplate ), private, pointer :: &
      Interactions => null ( )

contains


  subroutine Initialize ( WH, Name )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ), target :: &
      WH
    character ( * ), intent ( in )  :: &
      Name

    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit

    allocate ( WH % Header )
    associate ( WHH => WH % Header )

    !-- PositionSpace

    call WHH % InitializePositionSpace ( WH )

    select type ( PS => WH % PositionSpace )
    class is ( Atlas_SC_Form )

    !-- Prepare for Currents

    WH % N_CURRENTS_PS = 4
    allocate ( WH % Current_ASC_1D ( WH % N_CURRENTS_PS ) )
    allocate ( WH % TimeStepLabel ( WH % N_CURRENTS_PS ) )
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU ) &
      = 'Neutrinos_E_Nu'
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU_BAR ) &
      = 'Neutrinos_E_NuBar'
    WH % TimeStepLabel ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) &
      = 'Neutrinos_MuTau_NuNuBar'
    WH % TimeStepLabel ( WH % FLUID ) &
      = 'Fluid'

    !-- FIXME: This does not account for curvilinear coordinates
    MomentumDensity_U_Unit  =  WHH % MomentumUnit
    MomentumDensity_D_Unit  =  WHH % MomentumUnit

    !-- Electron Neutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU ) % Element )
    select type ( RA_E => WH % Current_ASC_1D &
                               ( WH % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_E % Initialize &
           ( PS, 'NEUTRINOS_E_NU', NameShortOption = 'Neutrinos_E_Nu', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )

    !-- Electron Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    select type ( RA_E_BAR => WH % Current_ASC_1D &
                                ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_E_BAR % Initialize &
           ( PS, 'NEUTRINOS_E_NU_BAR', NameShortOption = 'Neutrinos_E_Nu_Bar', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )

    !-- Mu and Tau Neutrinos and Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    select type ( RA_MU_TAU => WH % Current_ASC_1D &
                                 ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_MU_TAU % Initialize &
           ( PS, 'NEUTRINOS_MU_TAU_NU_NU_BAR', &
             NameShortOption = 'Neutrinos_MuTau_NuNuBar', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )

    !-- Fluid

    allocate ( Fluid_ASC_Form :: WH % Current_ASC_1D ( WH % FLUID ) % Element )
    select type ( FA => WH % Current_ASC_1D ( WH % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    call WHH % InitializeFluid ( FA )
    end select !-- FA

    !-- Interactions

    allocate ( WH % Interactions_ASC )
    associate ( IA => WH % Interactions_ASC )
    call IA % Initialize &
           ( PS, 'NEUTRINO_MOMENTS_GREY_1', &
             LengthUnitOption = WHH % CoordinateUnit ( 1 ), &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit )
    call RA_E % SetInteractions ( IA )
    call RA_E_BAR % SetInteractions ( IA )
    call RA_MU_TAU % SetInteractions ( IA )
    end associate !-- IA

    !-- Step

    allocate ( Step_RK2_C_ASC_1D_Form :: WH % Step )
    select type ( S => WH % Step )
    class is ( Step_RK2_C_ASC_1D_Form )

    call S % Initialize ( WH % Current_ASC_1D, Name )

    S % ApplySources_1D ( WH % NEUTRINOS_E_NU ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_E_NU ) % Pointer &
      =>  ApplyRelaxation_Interactions

    S % ApplySources_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplyRelaxation_Interactions

    S % ApplySources_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ApplyRelaxation_Interactions

    S % ApplySources_1D ( WH % FLUID ) % Pointer &
      =>  ApplySources_Fluid

    end select !-- S

    !-- Cleanup

    end select !-- RA_MU_TAU
    end select !-- RA_E_BAR
    end select !-- RA_E
    end select !-- PS
    end associate !-- WHH

  end subroutine Initialize


  impure elemental subroutine Finalize ( WH )

    type ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    if ( allocated ( WH % Interactions_ASC ) ) &
      deallocate ( WH % Interactions_ASC )
    if ( allocated ( WH % Header ) ) &
      deallocate ( WH % Header )

    call WH % FinalizeTemplate_C_1D_PS ( )

  end subroutine Finalize


  subroutine ApplySources_Radiation ( S, Increment, Radiation, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep

    !-- No sources applied here; just an occasion to compute interactions
    !   to be used in relaxation, and set a pointer to the Radiation increment.

    call Interactions % Compute ( Radiation )

    select case ( trim ( Radiation % Type ) )
    case ( 'NEUTRINOS_E_NU' )
      Increment_E_Nu => Increment 
    case ( 'NEUTRINOS_E_NU_BAR' )
      Increment_E_NuBar => Increment 
    case ( 'NEUTRINOS_MU_TAU_NU_NU_BAR' )
      Increment_MuTau_Nu_NuBar => Increment
    case default
      call Show ( 'Radiation Type not recognized', CONSOLE % ERROR )
      call Show ( 'WoosleyHeger_07_G__Form', 'module', CONSOLE % ERROR )
      call Show ( 'ApplySources_Radiation', 'subroutine', CONSOLE % ERROR )
    end select !-- Radiation % Type

  end subroutine ApplySources_Radiation

  
  subroutine ApplySources_Fluid ( S, Increment, Fluid, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep

    ! integer ( KDI ) :: &
    !   iEnergy_F, &
    !   iEnergy_R, &
    !   iMomentum_1_F, &
    !   iMomentum_1_R

    call ApplySourcesGravity ( S, Increment, Fluid, TimeStep )

    ! select type ( F => Fluid )
    ! class is ( Fluid_P_NR_Form )

    ! select type ( Chart => S % Chart )
    ! class is ( Chart_SL_Template )

    ! associate &
    !   ( I => Interactions, &
    !     R => Radiation )

    ! call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy_F )
    ! call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1_F )
    ! call Search ( R % iaConserved, R % CONSERVED_ENERGY_DENSITY, iEnergy_R )
    ! call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), &
    !               iMomentum_1_R )

    ! !-- Taking shortcuts on conserved vs. comoving here
    ! call ApplySourcesKernel &
    !        ( Increment % Value ( :, iEnergy_F ), &
    !          Increment % Value ( :, iMomentum_1_F ), &
    !          Chart % IsProperCell, &
    !          I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
    !          I % Value ( :, I % EFFECTIVE_OPACITY ), &
    !          I % Value ( :, I % TRANSPORT_OPACITY ), &
    !          R % Value ( :, R % COMOVING_ENERGY_DENSITY ), &
    !          RadiationIncrement % Value ( :, iEnergy_R ), &
    !          R % Value ( :, R % CONSERVED_MOMENTUM_DENSITY_D ( 1 ) ), &
    !          RadiationIncrement % Value ( :, iMomentum_1_R ), &
    !          CONSTANT % SPEED_OF_LIGHT, TimeStep ) 

    ! end associate !-- I, etc.
    ! end select !-- Chart
    ! end select !-- F

  end subroutine ApplySources_Fluid

  
  subroutine ApplySourcesKernel &
               ( KVE, KVM_1, IsProperCell, ED, EO, TO, J, dJ, H, dH, c, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVE, &
      KVM_1
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      ED, &
      EO, &
      TO, &
      J,  &
      dJ, &
      H,  &
      dH
    real ( KDR ) :: &
      c, &
      dT

    ! integer ( KDI ) :: &
    !   iV, &
    !   nV

    ! nV = size ( KVE )

    ! !$OMP parallel do private ( iV )
    ! do iV = 1, nV
    !   if ( .not. IsProperCell ( iV ) ) &
    !     cycle
    !   KVE ( iV )  &
    !     =  KVE ( iV )  -  c * dT  *  EO ( iV )  &
    !                       *  ( ED ( iV )  -  ( J ( iV ) + dJ ( iV ) ) ) 
    !   KVM_1 ( iV )  &
    !     =  KVM_1 ( iV )  +  c * dT  *  TO ( iV )  *  ( H ( iV ) + dH ( iV ) )
    ! end do
    ! !$OMP end parallel do

  end subroutine ApplySourcesKernel


end module WoosleyHeger_07_G__Form
