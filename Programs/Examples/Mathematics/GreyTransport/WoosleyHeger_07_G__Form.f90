module WoosleyHeger_07_G__Form

  !-- WoosleyHeger_07_Grey__Form

  use Mathematics
  use Fluid_ASC__Form
  use WoosleyHeger_07_Header_Form
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
      ApplySources_Radiation

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
    select type ( RA => WH % Current_ASC_1D &
                          ( WH % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA % Initialize &
           ( PS, 'NEUTRINOS_E_NU', NameShortOption = 'Neutrinos_E_Nu', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )
    end select !-- RA

    !-- Electron Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    select type ( RA => WH % Current_ASC_1D &
                          ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA % Initialize &
           ( PS, 'NEUTRINOS_E_NU_BAR', NameShortOption = 'Neutrinos_E_Nu_Bar', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )
    end select !-- RA

    !-- Mu and Tau Neutrinos and Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    select type ( RA => WH % Current_ASC_1D &
                          ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA % Initialize &
           ( PS, 'NEUTRINOS_MU_TAU_NU_NU_BAR', &
             NameShortOption = 'Neutrinos_MuTau_NuNuBar', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )
    end select !-- RA

    !-- Fluid

    allocate ( Fluid_ASC_Form :: WH % Current_ASC_1D ( WH % FLUID ) % Element )
    select type ( FA => WH % Current_ASC_1D ( WH % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    call WHH % InitializeFluid ( FA )
    end select !-- FA

    !-- Step

    allocate ( Step_RK2_C_ASC_1D_Form :: WH % Step )
    select type ( S => WH % Step )
    class is ( Step_RK2_C_ASC_1D_Form )

    call S % Initialize ( WH % Current_ASC_1D, Name )

    S % ApplySources_1D ( WH % NEUTRINOS_E_NU ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplySources_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplySources_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    ! S % ApplySources_1D ( WH % FLUID ) % Pointer &
    !   =>  ApplySources_Fluid
    ! S % ApplyRelaxation_1D ( WH % RADIATION ) % Pointer &
    !   =>  ApplyRelaxation_Interactions

    end select !-- S

    !-- Cleanup

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

  
end module WoosleyHeger_07_G__Form
