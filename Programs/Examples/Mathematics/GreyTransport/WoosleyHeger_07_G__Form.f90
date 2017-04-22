module WoosleyHeger_07_G__Form

  !-- WoosleyHeger_07_Grey__Form

  use Basics
  use Mathematics
  use Fluid_P_MHN__Form
  use Fluid_ASC__Form
  use WoosleyHeger_07_Header_Form
  use RadiationMoments_Form
  use NeutrinoMoments_Form
  use RadiationMoments_ASC__Form
  use Interactions_Template
  use Interactions_NM_G_1__Form
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
    procedure, public, pass :: &
      SetWriteTimeInterval
    final :: &
      Finalize
  end type WoosleyHeger_07_G_Form

    private :: &
      PrepareInteractions, &
      SetRadiation, &
      ApplySources_Radiation, &
      ApplySources_Fluid, &
      ComputeFluidSource_Radiation

      private :: &
        ApplySources_Fluid_Kernel, &
        ComputeFluidSource_G_S_Radiation_Kernel, &
        ComputeFluidSource_DP_Radiation_Kernel

    type ( VariableGroupForm ), private, allocatable :: &
      FluidSource_Radiation
    class ( Fluid_P_MHN_Form ), pointer :: &
      Fluid => null ( )
    class ( Interactions_NM_G_1_Form ), private, pointer :: &
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
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU ) = 'Nu_E'
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU_BAR ) = 'NuBar_E'
    WH % TimeStepLabel ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) = 'Nu_X'
    WH % TimeStepLabel ( WH % FLUID ) = 'Fluid'

    !-- FIXME: This does not account for curvilinear coordinates
    MomentumDensity_U_Unit  =  WHH % EnergyDensityUnit / UNIT % SPEED_OF_LIGHT
    MomentumDensity_D_Unit  =  WHH % EnergyDensityUnit / UNIT % SPEED_OF_LIGHT

    !-- Electron Neutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU ) % Element )
    select type ( RA_E => WH % Current_ASC_1D &
                            ( WH % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_E % Initialize &
           ( PS, 'NEUTRINOS_E_NU', NameShortOption = 'Nu_E', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit, &
             DegeneracyInfinityOption = 3.0_KDR )

    !-- Electron Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    select type ( RA_E_Bar => WH % Current_ASC_1D &
                                ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_E_Bar % Initialize &
           ( PS, 'NEUTRINOS_E_NU_BAR', NameShortOption = 'NuBar_E', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit, &
             DegeneracyInfinityOption = 2.0_KDR )

    !-- Mu and Tau Neutrinos and Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    select type &
      ( RA_MuTau => WH % Current_ASC_1D &
                       ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_MuTau % Initialize &
           ( PS, 'NEUTRINOS_MU_TAU_NU_NU_BAR', &
             NameShortOption = 'Nu_X', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit, &
             DegeneracyInfinityOption = 0.0_KDR )

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
    call RA_E_Bar % SetInteractions ( IA )
    call RA_MuTau % SetInteractions ( IA )
    call PrepareInteractions ( WH )
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
    S % HarvestIncrement_1D ( WH % NEUTRINOS_E_NU ) % Pointer &
      =>  ComputeFluidSource_Radiation

    S % ApplySources_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplyRelaxation_Interactions
    S % HarvestIncrement_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ComputeFluidSource_Radiation

    S % ApplySources_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ApplyRelaxation_Interactions
    S % HarvestIncrement_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
      =>  ComputeFluidSource_Radiation

    S % ApplySources_1D ( WH % FLUID ) % Pointer &
      =>  ApplySources_Fluid

    end select !-- S

    !-- Initial conditions

    call WHH % SetFluid ( )
    call SetRadiation ( WH )

    !-- Integrator

    call WHH % SetIntegratorParameters ( )

    call WH % InitializeTemplate_C_1D_PS &
           ( Name, TimeUnitOption = WHH % TimeUnit, &
             FinishTimeOption = WHH % FinishTime, &
             CourantFactorOption = WHH % CourantFactor, &
             nWriteOption = WHH % nWrite )

    !-- Cleanup

    end select !-- RA_MuTau
    end select !-- RA_E_Bar
    end select !-- RA_E
    end select !-- PS
    end associate !-- WHH

  end subroutine Initialize


  subroutine SetWriteTimeInterval ( I )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      I

    call I % Header % SetWriteTimeInterval ( )

  end subroutine SetWriteTimeInterval


  impure elemental subroutine Finalize ( WH )

    type ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    if ( allocated ( WH % Interactions_ASC ) ) &
      deallocate ( WH % Interactions_ASC )
    if ( allocated ( WH % Header ) ) &
      deallocate ( WH % Header )

    call WH % FinalizeTemplate_C_1D_PS ( )

  end subroutine Finalize


  subroutine PrepareInteractions ( WH )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    class ( Fluid_P_MHN_Form ), pointer :: &
      F
    class ( NeutrinoMomentsForm ), pointer :: &
      R_E, &
      R_E_Bar, &
      R_MuTau
    class ( InteractionsTemplate ), pointer :: &
      I

    associate ( IA => WH % Interactions_ASC )

    select type &
      ( RA_E => WH % Current_ASC_1D &
                  ( WH % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )

    select type &
      ( RA_E_Bar => WH % Current_ASC_1D &
                      ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )

    select type &
      ( RA_MuTau => WH % Current_ASC_1D &
                       ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )

    select type ( FA => WH % Current_ASC_1D ( WH % FLUID ) % Element )
    class is ( Fluid_ASC_Form )

    I        =>  IA % Interactions ( )
    R_E      =>  RA_E % NeutrinoMoments ( )
    R_E_Bar  =>  RA_E_Bar % NeutrinoMoments ( )
    R_MuTau  =>  RA_MuTau % NeutrinoMoments ( )
    F        =>  FA % Fluid_P_MHN ( )

    select type ( I )
    type is ( Interactions_NM_G_1_Form )
      call I % Set ( R_E, R_E_Bar, R_MuTau, F )
      Interactions => I  !-- Module variable for computing fluid source
    class default
      call Show ( 'Interactions type not recognized', CONSOLE % ERROR )
      call Show ( 'WoosleyHeger_07_G__Form', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareInteractions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I

    !-- Module variable for computing fluid source
    Fluid => F
    
    !-- Storage for fluid source terms

    allocate ( FluidSource_Radiation )
    call FluidSource_Radiation % Initialize &
           ( [ F % nValues, F % N_CONSERVED ], ClearOption = .true. )

    end select !-- FA
    end select !-- RA_MuTau
    end select !-- RA_E_Bar
    end select !-- RA_E
    end associate !-- IA
    nullify ( F, R_E, R_E_Bar, R_MuTau, I )

  end subroutine PrepareInteractions


  subroutine SetRadiation ( WH )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( NeutrinoMomentsForm ), pointer :: &
      R_E, &
      R_E_Bar, &
      R_MuTau    

    select type ( PS => WH % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type &
      ( RA_E => WH % Current_ASC_1D &
                  ( WH % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )

    select type &
      ( RA_E_Bar => WH % Current_ASC_1D &
                      ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )

    select type &
      ( RA_MuTau => WH % Current_ASC_1D &
                       ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )

    R_E      =>  RA_E % NeutrinoMoments ( )
    R_E_Bar  =>  RA_E_Bar % NeutrinoMoments ( )
    R_MuTau  =>  RA_MuTau % NeutrinoMoments ( )

    !-- No initial radiation, but set eigenspeeds
    call R_E % ComputeFromPrimitive ( G )
    call R_E_Bar % ComputeFromPrimitive ( G )
    call R_MuTau % ComputeFromPrimitive ( G )

    end select !-- RA_MuTau
    end select !-- RA_E_Bar
    end select !-- RA_E
    end select !-- PS
    nullify ( G, R_E, R_E_Bar, R_MuTau )

  end subroutine SetRadiation


  subroutine ApplySources_Radiation ( S, Increment, Radiation, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep

    call Show ( 'ApplySources_Radiation', CONSOLE % INFO_4 )

    call ApplySourcesCurvilinear_RadiationMoments &
           ( S, Increment, Radiation, TimeStep )

    !-- An occasion to compute the interactions to be used in relaxation
    call Interactions % Compute ( Radiation )

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

    integer ( KDI ) :: &
      iEnergy, &
      iMomentum_1, &
      iMomentum_2, &
      iMomentum_3, &
      iProton, &
      iEntropy
integer ( KDI ) :: &
  dG_G_maxloc, &
  dDP_DP_maxloc, &
  dDS_DS_maxloc
real ( KDR ) :: &
  dG_G_maxval, &
  dDP_DP_maxval, &
  dDS_DS_maxval, &
  Threshold
real ( KDR ), dimension ( Fluid % nValues ) :: &
  dG_G, &
  dDP_DP, &
  dDS_DS

    call Show ( 'ApplySources_Fluid', CONSOLE % INFO_4 )

    call ApplySourcesGravity ( S, Increment, Fluid, TimeStep )

    select type ( F => Fluid )
    class is ( Fluid_P_MHN_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 3 ), iMomentum_3 )
    call Search ( F % iaConserved, F % CONSERVED_PROTON_DENSITY, iProton )
    call Search ( F % iaConserved, F % CONSERVED_ENTROPY, iEntropy )

    call ApplySources_Fluid_Kernel &
           ( Increment % Value ( :, iEnergy ), &
             Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Increment % Value ( :, iMomentum_3 ), &
             Increment % Value ( :, iProton ), &
             Increment % Value ( :, iEntropy ), &
             Chart % IsProperCell, &
             FluidSource_Radiation % Value ( :, iEnergy ), &
             FluidSource_Radiation % Value ( :, iMomentum_1 ), &
             FluidSource_Radiation % Value ( :, iMomentum_2 ), &
             FluidSource_Radiation % Value ( :, iMomentum_3 ), &
             FluidSource_Radiation % Value ( :, iProton ), &
             F % Value ( :, F % TEMPERATURE ) )

dG_G   = Increment % Value ( :, iEnergy ) &
         /  max ( abs ( F % Value ( :, F % CONSERVED_ENERGY ) ), &
            tiny ( 0.0_KDR ) )
dDP_DP = Increment % Value ( :, iProton ) &
         /  max ( abs ( F % Value ( :, F % CONSERVED_PROTON_DENSITY ) ), &
                  tiny ( 0.0_KDR ) )
dDS_DS = Increment % Value ( :, iEntropy ) &
         /  max ( abs ( F % Value ( :, F % CONSERVED_ENTROPY ) ), &
                  tiny ( 0.0_KDR ) )
!call Show ( dG_G, '>>> dG_G' )
!call Show ( dDP_DP, '>>> dDP_DP' )
!call Show ( dDS_DS, '>>> dDS_DS' )

Threshold = 1.0e-1_KDR
dG_G_maxval   = maxval ( abs ( dG_G ), mask = Chart % IsProperCell )
dDP_DP_maxval = maxval ( abs ( dDP_DP ), mask = Chart % IsProperCell )
dDS_DS_maxval = maxval ( abs ( dDS_DS ), mask = Chart % IsProperCell )
if ( dG_G_maxval > Threshold ) then
  dG_G_maxloc = maxloc ( abs ( dG_G ), dim = 1, mask = Chart % IsProperCell )
  call Show ( dG_G_maxval, '>>> dG_G_maxval', CONSOLE % ERROR )
  call Show ( dG_G_maxloc, '>>> dG_G_maxloc', CONSOLE % ERROR )
  call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
              CONSOLE % ERROR )
end if
if ( dDP_DP_maxval > Threshold ) then
  dDP_DP_maxloc &
    = maxloc ( abs ( dDP_DP ), dim = 1, mask = Chart % IsProperCell )
  call Show ( dDP_DP_maxval, '>>> dDP_DP_maxval', CONSOLE % ERROR )
  call Show ( dDP_DP_maxloc, '>>> dDP_DP_maxloc', CONSOLE % ERROR )
  call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
              CONSOLE % ERROR )
end if
if ( dDS_DS_maxval > Threshold ) then
  dDS_DS_maxloc &
    = maxloc ( abs ( dDS_DS ), dim = 1, mask = Chart % IsProperCell )
  call Show ( dDS_DS_maxval, '>>> dDS_DS_maxval', CONSOLE % ERROR )
  call Show ( dDS_DS_maxloc, '>>> dDS_DS_maxloc', CONSOLE % ERROR )
  call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
              CONSOLE % ERROR )
end if
 
    call Clear ( FluidSource_Radiation % Value )

    end select !-- Chart
    end select !-- F

  end subroutine ApplySources_Fluid

  
  subroutine ComputeFluidSource_Radiation ( S, Increment, Radiation, TimeStep )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iEnergy_F, &
      iMomentum_1_F, &
      iMomentum_2_F, &
      iMomentum_3_F, &
      iProton_F, &
      iEnergy_R, &
      iMomentum_1_R, &
      iMomentum_2_R, &
      iMomentum_3_R

    select type ( R => Radiation )
    class is ( NeutrinoMomentsForm )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    associate &
      ( I => Interactions, &
        F => Fluid )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy_F )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), &
                  iMomentum_1_F )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), &
                  iMomentum_2_F )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 3 ), &
                  iMomentum_3_F )
    call Search ( F % iaConserved, F % CONSERVED_PROTON_DENSITY, iProton_F )
    call Search ( R % iaConserved, R % CONSERVED_ENERGY_DENSITY, iEnergy_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), &
                  iMomentum_1_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_DENSITY_D ( 2 ), &
                  iMomentum_2_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_DENSITY_D ( 3 ), &
                  iMomentum_3_R )

    !-- Taking shortcuts on conserved vs. comoving here

    call ComputeFluidSource_G_S_Radiation_Kernel &
           ( FluidSource_Radiation % Value ( :, iEnergy_F ), & 
             FluidSource_Radiation % Value ( :, iMomentum_1_F ), &
             FluidSource_Radiation % Value ( :, iMomentum_2_F ), &
             FluidSource_Radiation % Value ( :, iMomentum_3_F ), &
             Chart % IsProperCell, &
             I % Value ( :, I % EMISSIVITY ), &
             I % Value ( :, I % EFFECTIVE_OPACITY ), &
             I % Value ( :, I % TRANSPORT_OPACITY ), &
             R % Value ( :, R % COMOVING_ENERGY_DENSITY ), &
             Increment % Value ( :, iEnergy_R ), &
             R % Value ( :, R % CONSERVED_MOMENTUM_DENSITY_D ( 1 ) ), &
             R % Value ( :, R % CONSERVED_MOMENTUM_DENSITY_D ( 2 ) ), &
             R % Value ( :, R % CONSERVED_MOMENTUM_DENSITY_D ( 3 ) ), &
             Increment % Value ( :, iMomentum_1_R ), &
             Increment % Value ( :, iMomentum_2_R ), &
             Increment % Value ( :, iMomentum_3_R ), &
             CONSTANT % SPEED_OF_LIGHT, TimeStep ) 

    select case ( trim ( Radiation % Type ) )
    case ( 'NEUTRINOS_E_NU' )
      call ComputeFluidSource_DP_Radiation_Kernel &
             ( FluidSource_Radiation % Value ( :, iProton_F ), & 
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_NUMBER ), &
               I % Value ( :, I % EFFECTIVE_OPACITY_NUMBER ), &
               R % Value ( :, R % COMOVING_ENERGY_DENSITY ), &
               Increment % Value ( :, iEnergy_R ), &
               CONSTANT % SPEED_OF_LIGHT, TimeStep, Sign = +1.0_KDR )
    case ( 'NEUTRINOS_E_NU_BAR' )
      call ComputeFluidSource_DP_Radiation_Kernel &
             ( FluidSource_Radiation % Value ( :, iProton_F ), & 
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_NUMBER ), &
               I % Value ( :, I % EFFECTIVE_OPACITY_NUMBER ), &
               R % Value ( :, R % COMOVING_ENERGY_DENSITY ), &
               Increment % Value ( :, iEnergy_R ), &
               CONSTANT % SPEED_OF_LIGHT, TimeStep, Sign = -1.0_KDR )
    end select !-- Radiation % Type

    end associate !-- I, etc.
    end select !-- Chart
    end select !-- R

  end subroutine ComputeFluidSource_Radiation


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, K_S_2, K_S_3, K_DP, K_DS, IsProperCell, FS_R_G, &
                 FS_R_S_1, FS_R_S_2, FS_R_S_3, FS_R_DP, T )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_G, &
      K_S_1, K_S_2, K_S_3, &
      K_DP, &
      K_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3, &
      FS_R_DP, &
      T

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU

    nV = size ( K_G )

    AMU = CONSTANT % ATOMIC_MASS_UNIT

! call Show ( K_G, '>>> K_G' )
! call Show ( FS_R_G, '>>> FS_R_G' )

! call Show ( K_S_1, '>>> K_S_1' )
! call Show ( FS_R_S_1, '>>> FS_R_S_1' )

! call Show ( K_DP, '>>> K_DP' )
! call Show ( FS_R_DP, '>>> FS_R_DP' )

! call Show ( K_DS, '>>> K_DS' )
! call Show ( AMU * FS_R_G / T, '>>> FS_R_DS' )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )
      K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )
      K_S_2 ( iV )  =  K_S_2 ( iV )  +  FS_R_S_2 ( iV )
      K_S_3 ( iV )  =  K_S_3 ( iV )  +  FS_R_S_3 ( iV )
      K_DP ( iV )   =  K_DP ( iV )   +  FS_R_DP ( iV )
      K_DS ( iV )   =  K_DS ( iV )   +  AMU * FS_R_G ( iV ) / T ( iV )
    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


  subroutine ComputeFluidSource_G_S_Radiation_Kernel &
               ( FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, IsProperCell, &
                 E, EO, TO, J, dJ, H_1, H_2, H_3, dH_1, dH_2, dH_3, &
                 c, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      E, &
      EO, &
      TO, &
      J,  &
      dJ, &
      H_1, H_2, H_3, &
      dH_1, dH_2, dH_3
    real ( KDR ) :: &
      c, &
      dT

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( FS_R_G )

! call Show ( FS_R_G, '>>> FS_R_G' )
! call Show ( E, '>>> E' )
! call Show ( c * dT, '>>> c * dT' )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      FS_R_G ( iV )  &
        =  FS_R_G ( iV )  &
           -  c * dT  *  ( E ( iV )  -  EO ( iV ) * ( J ( iV ) + dJ ( iV ) ) ) 
      FS_R_S_1 ( iV )  &
        =  FS_R_S_1 ( iV )  +  c * dT  *  TO ( iV )  &
                               *  ( H_1 ( iV ) + dH_1 ( iV ) )
      FS_R_S_2 ( iV )  &
        =  FS_R_S_2 ( iV )  +  c * dT  *  TO ( iV )  &
                               *  ( H_2 ( iV ) + dH_2 ( iV ) )
      FS_R_S_3 ( iV )  &
        =  FS_R_S_3 ( iV )  +  c * dT  *  TO ( iV )  &
                               *  ( H_3 ( iV ) + dH_3 ( iV ) )
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_G_S_Radiation_Kernel


  subroutine ComputeFluidSource_DP_Radiation_Kernel &
               ( FS_R_DP, IsProperCell, EN, EON, J, dJ, c, dT, Sign )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_DP
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      EN, &
      EON, &
      J,  &
      dJ
    real ( KDR ) :: &
      c, &
      dT, &
      Sign

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU

    nV = size ( FS_R_DP )

    AMU = CONSTANT % ATOMIC_MASS_UNIT

! call Show ( FS_R_DP, '>>> FS_R_DP' )
! call Show ( AMU * EN, '>>> AMU * EN' )
! call Show ( c * dT, '>>> c * dT' )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      FS_R_DP ( iV )  &
        =  FS_R_DP ( iV )  &
           -  Sign * c * dT * AMU &
              *  ( EN ( iV )  -  EON ( iV ) * ( J ( iV ) + dJ ( iV ) ) ) 
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_DP_Radiation_Kernel


end module WoosleyHeger_07_G__Form
