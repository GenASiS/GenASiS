module WoosleyHeger_07_G__Form

  !-- WoosleyHeger_07_Grey__Form

  use Basics
  use Mathematics
  use Fluid_P_MHN__Form
  use Fluid_ASC__Form
  use WoosleyHeger_07_Header_Form
  use RadiationMoments_Form
  use NeutrinoMoments_G__Form
  use RadiationMoments_ASC__Form
  use ApplyCurvilinear_RM__Command
  use Interactions_Template
  use Interactions_NM_1_G__Form
  use Interactions_ASC__Form
  use ApplyRelaxation_NM_G__Command

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ) :: &
    WoosleyHeger_07_G_Form
      integer ( KDI ) :: &
        NEUTRINOS_E_NU             = 1, &
        NEUTRINOS_E_NU_BAR         = 2, &
        ! NEUTRINOS_MU_TAU_NU_NU_BAR = 3, &
        ! FLUID                      = 4         
        FLUID                      = 3         
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
    procedure, private, pass :: &
      ComputeTimeStepLocal
  end type WoosleyHeger_07_G_Form

    private :: &
      PrepareInteractions, &
      SetRadiation, &
      ApplySources_Radiation, &
      ApplySources_Fluid, &
      ComputeFluidSource_Radiation, &
      SetRadiation_FluidVelocity

      private :: &
        ComputeTimeStep_SB_J_Kernel, &
        ComputeTimeStep_SB_N_Kernel, &
        ComputeTimeStep_S_H_Kernel, &
        ComputeTimeStep_YE_N_Kernel, &
        ComputeTimeStep_SB_Div_J_Kernel, &
        ComputeTimeStep_SB_Div_N_Kernel, &
        ApplySources_Fluid_Kernel, &
        ImposeBetaEquilibrium_Kernel, &
        ComputeFluidSource_G_S_Radiation_Kernel, &
        ComputeFluidSource_DP_Radiation_Kernel

    logical ( KDL ), private :: &
      ReachedEquilibrium = .false.
    type ( VariableGroupForm ), private, allocatable :: &
      FluidSource_Radiation
    class ( Fluid_P_MHN_Form ), private, pointer :: &
      Fluid => null ( )
    class ( NeutrinoMoments_G_Form ), private, pointer :: &
      Radiation_E     => null ( ), &
      Radiation_E_Bar => null ( )!, &
!      Radiation_MuTau => null ( )   
    class ( Interactions_NM_1_G_Form ), private, pointer :: &
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

!    WH % N_CURRENTS_PS = 4
    WH % N_CURRENTS_PS = 3
    allocate ( WH % Current_ASC_1D ( WH % N_CURRENTS_PS ) )
    allocate ( WH % TimeStepLabel ( WH % N_CURRENTS_PS + 3 ) )
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU ) = 'Nu_E'
    WH % TimeStepLabel ( WH % NEUTRINOS_E_NU_BAR ) = 'NuBar_E'
!    WH % TimeStepLabel ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) = 'Nu_X'
    WH % TimeStepLabel ( WH % FLUID ) = 'Fluid'
    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 1 ) = 'Interactions_SB_J'
    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 2 ) = 'Interactions_SB_N'
!    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 3 ) = 'Interactions_S_H'
    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 3 ) = 'Interactions_YE_N'
!    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 5 ) = 'Interactions_SB_Div_J'
!    WH % TimeStepLabel ( WH % N_CURRENTS_PS + 6 ) = 'Interactions_SB_Div_N'

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
             UseLimiterOption = .true., LimiterParameterOption = 1.0_KDR, &
             Velocity_U_UnitOption = WHH % VelocityUnit, &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )

    !-- Electron Antineutrinos

    allocate &
      ( RadiationMoments_ASC_Form :: &
          WH % Current_ASC_1D ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    select type ( RA_E_Bar => WH % Current_ASC_1D &
                                ( WH % NEUTRINOS_E_NU_BAR ) % Element )
    class is ( RadiationMoments_ASC_Form )
    call RA_E_Bar % Initialize &
           ( PS, 'NEUTRINOS_E_NU_BAR', NameShortOption = 'NuBar_E', &
             UseLimiterOption = .true., LimiterParameterOption = 1.0_KDR, &
             Velocity_U_UnitOption = WHH % VelocityUnit, &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
             TemperatureUnitOption = WHH % TemperatureUnit )

    ! !-- Mu and Tau Neutrinos and Antineutrinos

    ! allocate &
    !   ( RadiationMoments_ASC_Form :: &
    !       WH % Current_ASC_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    ! select type &
    !   ( RA_MuTau => WH % Current_ASC_1D &
    !                    ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    ! class is ( RadiationMoments_ASC_Form )
    ! call RA_MuTau % Initialize &
    !        ( PS, 'NEUTRINOS_MU_TAU_NU_NU_BAR', &
    !          NameShortOption = 'Nu_X', &
    !          MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
    !          MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
    !          EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
    !          TemperatureUnitOption = WHH % TemperatureUnit )

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
           ( PS, 'NEUTRINO_MOMENTS_1_GREY', &
             LengthUnitOption = WHH % CoordinateUnit ( 1 ), &
             EnergyDensityUnitOption = WHH % EnergyDensityUnit )
    call RA_E % SetInteractions ( IA )
    call RA_E_Bar % SetInteractions ( IA )
!    call RA_MuTau % SetInteractions ( IA )
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
      =>  ApplyRelaxation_NM_G
    S % HarvestIncrement_1D ( WH % NEUTRINOS_E_NU ) % Pointer &
      =>  ComputeFluidSource_Radiation

    S % ApplySources_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplySources_Radiation
    S % ApplyRelaxation_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ApplyRelaxation_NM_G
    S % HarvestIncrement_1D ( WH % NEUTRINOS_E_NU_BAR ) % Pointer &
      =>  ComputeFluidSource_Radiation

    ! S % ApplySources_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
    !   =>  ApplySources_Radiation
    ! S % ApplyRelaxation_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
    !   =>  ApplyRelaxation_NM_G
    ! S % HarvestIncrement_1D ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Pointer &
    !   =>  ComputeFluidSource_Radiation

    S % ApplySources_1D ( WH % FLUID ) % Pointer &
      =>  ApplySources_Fluid
    S % HarvestCurrent_1D ( WH % FLUID ) % Pointer &
      =>  SetRadiation_FluidVelocity

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

!    end select !-- RA_MuTau
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


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iTSC, &
      iEnergy, &
      iNumber
    class ( NeutrinoMoments_G_Form ), pointer :: &
      R_E

    call I % ComputeTimeStepLocalTemplate ( TimeStepCandidate )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( CSL => PS % Chart )
    class is ( Chart_SL_Template )

    select type &
      ( RA_E => I % Current_ASC_1D &
                  ( I % NEUTRINOS_E_NU ) % Element )
    class is ( RadiationMoments_ASC_Form )

    R_E  =>  RA_E % NeutrinoMoments_G ( )

    associate &
      ( I_R_E => Interactions, &
        F => Fluid )

    call I_R_E % Compute ( R_E )

    call ComputeTimeStep_SB_J_Kernel &
           ( CSL % IsProperCell, &
             I_R_E % Value ( :, I_R_E % EMISSIVITY_J ), &
             I_R_E % Value ( :, I_R_E % OPACITY_J ), &
             R_E % Value ( :, R_E % COMOVING_ENERGY ), &
             F % Value ( :, F % TEMPERATURE ), &
             F % Value ( :, F % ENTROPY_PER_BARYON ), &
             F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             TimeStepCandidate ( I % N_CURRENTS_PS + 1 ) )

    call ComputeTimeStep_SB_N_Kernel &
           ( CSL % IsProperCell, &
             I_R_E % Value ( :, I_R_E % EMISSIVITY_N ), &
             I_R_E % Value ( :, I_R_E % OPACITY_N ), &
             R_E % Value ( :, R_E % COMOVING_NUMBER ), &
             F % Value ( :, F % TEMPERATURE ), &
             F % Value ( :, F % ENTROPY_PER_BARYON ), &
             F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
             TimeStepCandidate ( I % N_CURRENTS_PS + 2 ) )

    ! call ComputeTimeStep_S_H_Kernel &
    !        ( CSL % IsProperCell, &
    !          I_R_E % Value ( :, I_R_E % EMISSIVITY_J ), &
    !          I_R_E % Value ( :, I_R_E % OPACITY_J ), &
    !          I_R_E % Value ( :, I_R_E % OPACITY_H ), &
    !          R_E % Value ( :, R_E % COMOVING_ENERGY ), &
    !          R_E % Value ( :, R_E % COMOVING_MOMENTUM_U ( 1 ) ), &
    !          R_E % Value ( :, R_E % FLUID_VELOCITY_U ( 1 ) ), &
    !          F % Value ( :, F % MOMENTUM_DENSITY_D ( 1 ) ), &
    !          TimeStepCandidate ( I % N_CURRENTS_PS + 3 ) )

    call ComputeTimeStep_YE_N_Kernel &
           ( CSL % IsProperCell, &
             I_R_E % Value ( :, I_R_E % EMISSIVITY_N ), &
             I_R_E % Value ( :, I_R_E % OPACITY_N ), &
             R_E % Value ( :, R_E % COMOVING_NUMBER ), &
             F % Value ( :, F % ELECTRON_FRACTION ), &
             F % Value ( :, F % BARYON_MASS ), &
             F % Value ( :, F % COMOVING_DENSITY ), &
             TimeStepCandidate ( I % N_CURRENTS_PS + 3 ) )

    ! call Search ( R_E % iaConserved, R_E % CONSERVED_ENERGY, iEnergy )
    ! call Search ( R_E % iaConserved, R_E % CONSERVED_NUMBER, iNumber )

    ! call ComputeTimeStep_SB_Div_J_Kernel &
    !        ( CSL % IsProperCell, &
    !          R_E % Sources % Value ( :, iEnergy ), &
    !          F % Value ( :, F % TEMPERATURE ), &
    !          F % Value ( :, F % ENTROPY_PER_BARYON ), &
    !          F % Value ( :, F % BARYON_MASS ), &
    !          F % Value ( :, F % COMOVING_DENSITY ), &
    !          TimeStepCandidate ( I % N_CURRENTS_PS + 5 ) )

    ! call ComputeTimeStep_SB_Div_N_Kernel &
    !        ( CSL % IsProperCell, &
    !          R_E % Sources % Value ( :, iNumber ), &
    !          F % Value ( :, F % TEMPERATURE ), &
    !          F % Value ( :, F % ENTROPY_PER_BARYON ), &
    !          F % Value ( :, F % BARYON_MASS ), &
    !          F % Value ( :, F % COMOVING_DENSITY ), &
    !          F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
    !          F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
    !          TimeStepCandidate ( I % N_CURRENTS_PS + 6 ) )

    if ( ( mod ( I % iCycle, 1000 ) == 0 &
           .or. any ( TimeStepCandidate < 1.0e-12_KDR * UNIT % SECOND ) ) &
         .and. &
         minloc ( TimeStepCandidate, dim = 1 ) > I % N_CURRENTS_PS ) &
    then
      call Show ( I % iCycle, '>>> iCycle', I % IGNORABILITY )
      do iTSC = 1, I % nTimeStepCandidates
        call Show ( TimeStepCandidate ( iTSC ), I % TimeUnit, &
                    trim ( I % TimeStepLabel ( iTSC ) ) // ' TimeStep', &
                    I % IGNORABILITY )
      end do !-- iTSC
    end if

    end associate !-- I_R_E, etc.
    end select !-- RA_E
    end select !-- CSL
    end select !-- PS
    nullify ( R_E )

  end subroutine ComputeTimeStepLocal


  subroutine PrepareInteractions ( WH )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    class ( Fluid_P_MHN_Form ), pointer :: &
      F
    class ( NeutrinoMoments_G_Form ), pointer :: &
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

    ! select type &
    !   ( RA_MuTau => WH % Current_ASC_1D &
    !                    ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    ! class is ( RadiationMoments_ASC_Form )

    select type ( FA => WH % Current_ASC_1D ( WH % FLUID ) % Element )
    class is ( Fluid_ASC_Form )

    I        =>  IA % Interactions ( )
    R_E      =>  RA_E % NeutrinoMoments_G ( )
    R_E_Bar  =>  RA_E_Bar % NeutrinoMoments_G ( )
    R_MuTau  =>  null ( )
!    R_MuTau  =>  RA_MuTau % NeutrinoMoments_G ( )
    F        =>  FA % Fluid_P_MHN ( )

    select type ( I )
    type is ( Interactions_NM_1_G_Form )
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
!    end select !-- RA_MuTau
    end select !-- RA_E_Bar
    end select !-- RA_E
    end associate !-- IA
!    nullify ( F, R_E, R_E_Bar, R_MuTau, I )
    nullify ( F, R_E, R_E_Bar, I )

  end subroutine PrepareInteractions


  subroutine SetRadiation ( WH )

    class ( WoosleyHeger_07_G_Form ), intent ( inout ) :: &
      WH

    class ( GeometryFlatForm ), pointer :: &
      G

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

    ! select type &
    !   ( RA_MuTau => WH % Current_ASC_1D &
    !                    ( WH % NEUTRINOS_MU_TAU_NU_NU_BAR ) % Element )
    ! class is ( RadiationMoments_ASC_Form )

    Radiation_E      =>  RA_E % NeutrinoMoments_G ( )
    Radiation_E_Bar  =>  RA_E_Bar % NeutrinoMoments_G ( )
!    Radiation_MuTau  =>  RA_MuTau % NeutrinoMoments_G ( )

    !-- No initial radiation, but set eigenspeeds
    call Radiation_E % ComputeFromPrimitive ( G )
    call Radiation_E_Bar % ComputeFromPrimitive ( G )
!    call Radiation_MuTau % ComputeFromPrimitive ( G )

!    end select !-- RA_MuTau
    end select !-- RA_E_Bar
    end select !-- RA_E
    end select !-- PS
    nullify ( G )

  end subroutine SetRadiation


  subroutine ApplySources_Radiation &
               ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call Show ( 'ApplySources_Radiation', CONSOLE % INFO_4 )

    call ApplyCurvilinear_RM &
           ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

    !-- An occasion to compute the interactions to be used in relaxation
    call Interactions % Compute ( Radiation )
!    call Interactions % Regulate ( Radiation, TimeStep )

  end subroutine ApplySources_Radiation

  
  subroutine ApplySources_Fluid &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

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

    call ApplyGravity ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

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
             FluidSource_Radiation % Value ( :, iEntropy ), &
             F % Value ( :, F % TEMPERATURE ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_E ) )

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
    class ( CurrentTemplate ), intent ( inout ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep

    integer ( KDI ) :: &
      iEnergy_F, &
      iMomentum_1_F, &
      iMomentum_2_F, &
      iMomentum_3_F, &
      iProton_F, &
      iEntropy_F, &
      iEnergy_R, &
      iMomentum_1_R, &
      iMomentum_2_R, &
      iMomentum_3_R, &
      iNumber_R
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( R => Radiation )
    class is ( NeutrinoMoments_G_Form )

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
    call Search ( F % iaConserved, F % CONSERVED_ENTROPY, iEntropy_F )
    call Search ( R % iaConserved, R % CONSERVED_ENERGY, iEnergy_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_D ( 2 ), &
                  iMomentum_2_R )
    call Search ( R % iaConserved, R % CONSERVED_MOMENTUM_D ( 3 ), &
                  iMomentum_3_R )
    call Search ( R % iaConserved, R % CONSERVED_NUMBER, iNumber_R )

    G => Chart % Geometry ( )

    !-- Taking shortcuts on conserved vs. comoving here

    ! if ( trim ( Radiation % Type ) == 'NEUTRINOS_E_NU' ) &
    !   call ImposeBetaEquilibrium_Kernel &
    !          ( R % Value ( :, R % BETA_EQUILIBRIUM ), &
    !            Increment % Value ( :, iEnergy_R ), &
    !            Increment % Value ( :, iNumber_R ), &
    !            R, &
    !            F % Value ( :, F % BARYON_MASS ), &
    !            F % Value ( :, F % COMOVING_DENSITY ), &
    !            F % Value ( :, F % TEMPERATURE ), &
    !            F % Value ( :, F % ENTROPY_PER_BARYON ), &
    !            F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
    !            F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
    !            R % Value ( :, R % COMOVING_ENERGY ), &
    !            R % Value ( :, R % COMOVING_ENERGY_EQ ), &
    !            R % Value ( :, R % COMOVING_NUMBER ), &
    !            R % Value ( :, R % COMOVING_NUMBER_EQ ) )

    call ComputeFluidSource_G_S_Radiation_Kernel &
           ( FluidSource_Radiation % Value ( :, iEnergy_F ), & 
             FluidSource_Radiation % Value ( :, iMomentum_1_F ), &
             FluidSource_Radiation % Value ( :, iMomentum_2_F ), &
             FluidSource_Radiation % Value ( :, iMomentum_3_F ), &
             FluidSource_Radiation % Value ( :, iEntropy_F ), &
             Chart % IsProperCell, &
             I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             R % Value ( :, R % COMOVING_ENERGY ), &
             Increment % Value ( :, iEnergy_R ), &
             R % Value ( :, R % COMOVING_MOMENTUM_U ( 1 ) ), &
             R % Value ( :, R % COMOVING_MOMENTUM_U ( 2 ) ), &
             R % Value ( :, R % COMOVING_MOMENTUM_U ( 3 ) ), &
             Increment % Value ( :, iMomentum_1_R ), &
             Increment % Value ( :, iMomentum_2_R ), &
             Increment % Value ( :, iMomentum_3_R ), &
             R % Value ( :, R % BETA_EQUILIBRIUM ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 1 ) ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 2 ) ), &
             R % Value ( :, R % FLUID_VELOCITY_U ( 3 ) ), &
             F % Value ( :, F % TEMPERATURE ), &
             G % Value ( :, G % METRIC_DD_22 ), &
             G % Value ( :, G % METRIC_DD_33 ), &
             CONSTANT % SPEED_OF_LIGHT, TimeStep ) 

    select case ( trim ( Radiation % Type ) )
    case ( 'NEUTRINOS_E_NU' )
      call ComputeFluidSource_DP_Radiation_Kernel &
             ( FluidSource_Radiation % Value ( :, iProton_F ), & 
               FluidSource_Radiation % Value ( :, iEntropy_F ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_N ), &
               I % Value ( :, I % OPACITY_N ), &
               R % Value ( :, R % COMOVING_NUMBER ), &
               Increment % Value ( :, iNumber_R ), &
               R % Value ( :, R % BETA_EQUILIBRIUM ), &
               F % Value ( :, F % TEMPERATURE ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
               CONSTANT % SPEED_OF_LIGHT, TimeStep, Sign = +1.0_KDR )
    case ( 'NEUTRINOS_E_NU_BAR' )
      call ComputeFluidSource_DP_Radiation_Kernel &
             ( FluidSource_Radiation % Value ( :, iProton_F ), & 
               FluidSource_Radiation % Value ( :, iEntropy_F ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_N ), &
               I % Value ( :, I % OPACITY_N ), &
               R % Value ( :, R % COMOVING_NUMBER ), &
               Increment % Value ( :, iNumber_R ), &
               R % Value ( :, R % BETA_EQUILIBRIUM ), &
               F % Value ( :, F % TEMPERATURE ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
               CONSTANT % SPEED_OF_LIGHT, TimeStep, Sign = -1.0_KDR )
    end select !-- Radiation % Type

    end associate !-- I, etc.
    end select !-- Chart
    end select !-- R

    nullify ( G )
    
  end subroutine ComputeFluidSource_Radiation


  subroutine SetRadiation_FluidVelocity ( S, Fluid )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid

    select type ( F => Fluid )
    class is ( Fluid_P_MHN_Form )

    associate &
      ( R_E     => Radiation_E, &
        R_E_Bar => Radiation_E_Bar )!, &
!        R_MuTau => Radiation_MuTau )
    call Copy ( F % Value ( :, F % VELOCITY_U ( 1 ) : F % VELOCITY_U ( 3 ) ), &
                R_E % Value ( :, R_E % FLUID_VELOCITY_U ( 1 ) &
                               : R_E % FLUID_VELOCITY_U ( 3 ) ) )
    call Copy ( F % Value ( :, F % VELOCITY_U ( 1 ) : F % VELOCITY_U ( 3 ) ), &
                R_E_Bar % Value ( :, R_E_Bar % FLUID_VELOCITY_U ( 1 ) &
                                  : R_E_Bar % FLUID_VELOCITY_U ( 3 ) ) )
!    call Copy ( F % Value ( :, F % VELOCITY_U ( 1 ) : F % VELOCITY_U ( 3 ) ), &
!                R_MuTau % Value ( :, R_MuTau % FLUID_VELOCITY_U ( 1 ) &
!                                  : R_MuTau % FLUID_VELOCITY_U ( 3 ) ) )
    end associate !-- R_E, etc.

    end select !-- F

  end subroutine SetRadiation_FluidVelocity


  subroutine ComputeTimeStep_SB_J_Kernel &
               ( IsProperCell, Xi_J, Chi_J, J, T, SB, MB, NB, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, Chi_J, J, &
      T, SB, MB, NB
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Xi_J )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  abs ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) )  &
                  /  ( T ( iV )  *  SB ( iV )  &
                       *  MB ( iV )  *  NB ( iV )  /  AMU  ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_SB_J_Kernel


  subroutine ComputeTimeStep_SB_N_Kernel &
               ( IsProperCell, Xi_N, Chi_N, N, T, SB, MB, NB, Mu_e, Mu_n_p, &
                 TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_N, Chi_N, N, &
      T, SB, MB, NB, Mu_e, Mu_n_p
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Xi_N )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  ( Mu_e ( iV )  -  Mu_n_p ( iV ) )  &
                  *  abs ( Xi_N ( iV )  -  Chi_N ( iV ) * N ( iV ) )  &
                  /  ( T ( iV )  *  SB ( iV )  &
                       *  MB ( iV )  *  NB ( iV )  /  AMU  ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_SB_N_Kernel


  subroutine ComputeTimeStep_S_H_Kernel &
               ( IsProperCell, Xi_J, Chi_J, Chi_H, J, H_U_1, V_1, S_D_1, &
                 TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, Chi_J, Chi_H, &
      J, H_U_1, &
      V_1, &
      S_D_1
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Chi_H )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  abs ( - Chi_H ( iV ) * H_U_1 ( iV )  &
                        +  V_1 ( iV )  &
                           *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )  &
                  /  abs ( S_D_1 ( iV ) ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_S_H_Kernel


  subroutine ComputeTimeStep_YE_N_Kernel &
               ( IsProperCell, Xi_N, Chi_N, N, Y_e, MB, NB, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_N, Chi_N, N, &
      Y_e, MB, NB
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Xi_N )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  abs ( Xi_N ( iV )  -  Chi_N ( iV ) * N ( iV ) )  &
                  /  ( Y_e ( iV )  *  MB ( iV )  *  NB ( iV )  /  AMU  ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_YE_N_Kernel


  subroutine ComputeTimeStep_SB_Div_J_Kernel &
               ( IsProperCell, Div_J, T, SB, MB, NB, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Div_J, &
      T, SB, MB, NB
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Div_J )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  Div_J ( iV )  &
                  /  ( T ( iV )  *  SB ( iV )  &
                       *  MB ( iV )  *  NB ( iV )  /  AMU  ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_SB_Div_J_Kernel


  subroutine ComputeTimeStep_SB_Div_N_Kernel &
               ( IsProperCell, Div_N, T, SB, MB, NB, Mu_e, Mu_n_p, TimeStep )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Div_N, &
      T, SB, MB, NB, Mu_e, Mu_n_p
    real ( KDR ), intent ( inout ) :: &
      TimeStep    

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      f, &  !-- factor
      TimeStepInverse

    nV = size ( Div_N )

    AMU  =  CONSTANT % ATOMIC_MASS_UNIT
    f    =  0.01_KDR

    TimeStepInverse = - huge ( 0.0_KDR )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( IsProperCell ( iV ) ) then
        TimeStepInverse &
          = max ( TimeStepInverse, &
                  ( Mu_e ( iV )  -  Mu_n_p ( iV ) ) * Div_N ( iV )  &
                  /  ( T ( iV )  *  SB ( iV )  &
                       *  MB ( iV )  *  NB ( iV )  /  AMU  ) )
      end if
    end do
    !$OMP end parallel do

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, f * 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStep_SB_Div_N_Kernel


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, K_S_2, K_S_3, K_DP, K_DS, IsProperCell, FS_R_G, &
                 FS_R_S_1, FS_R_S_2, FS_R_S_3, FS_R_DP, FS_R_DS, &
                 T, Mu_n_p, Mu_e )

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
      FS_R_DS, &
      T, &
      Mu_n_p, Mu_e

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
! call Show ( AMU * FS_R_G / T, '>>> FS_R_DS_G' )
! call Show ( - ( Mu_e  -  Mu_n_p ) * FS_R_DP / T, '>>> FS_R_DS_DP' )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )
      K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )
      K_S_2 ( iV )  =  K_S_2 ( iV )  +  FS_R_S_2 ( iV )
      K_S_3 ( iV )  =  K_S_3 ( iV )  +  FS_R_S_3 ( iV )
      K_DP ( iV )   =  K_DP ( iV )   +  FS_R_DP ( iV )
      K_DS ( iV )   =  K_DS ( iV )   +  FS_R_DS ( iV )

      ! K_DS ( iV )   &
      !   =  K_DS ( iV )   &
      !      +  ( AMU * FS_R_G ( iV ) &
      !           -  ( Mu_e ( iV )  -  Mu_n_p ( iV ) )  &
      !              *  FS_R_DP ( iV ) ) &
      !         /  T ( iV )

    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


  subroutine ImposeBetaEquilibrium_Kernel &
               ( Beta_EQ, dJ, dN, NM, M_F, N_F, T_F, S_F, Mu_e, Mu_n_p, &
                 J, J_EQ, N, N_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Beta_EQ, &
      dJ, &
      dN
    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      NM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M_F, N_F, T_F, S_F, Mu_E, Mu_N_P, &
      J, J_EQ, &
      N, N_EQ

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU, &
      Rho_EQ, &
      f, &         !-- fraction
      r_J, r_N, &  !-- regulation factor
      Delta_J, Delta_N

    nV = size ( Beta_EQ )

    AMU     =  CONSTANT % ATOMIC_MASS_UNIT
    Rho_EQ  =  1.0e13_KDR  *  UNIT % GRAM  *  UNIT % CENTIMETER ** (-3)
    f       =  0.01_KDR

    !$OMP parallel do private ( iV, Delta_J, Delta_N, r_J, r_N )
    do iV = 1, nV
      if ( M_F ( iV )  *  N_F ( iV )  >  Rho_EQ ) then

        if ( .not. ReachedEquilibrium ) then
          ReachedEquilibrium = .true.
call Show ( '>>> Reached beta equilibrium', CONSOLE % ERROR )
call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
        end if

        Beta_EQ ( iV )  =  1.0_KDR

        Delta_J  =  J_EQ ( iV )  -  J ( iV )
        Delta_N  =  N_EQ ( iV )  -  N ( iV )

        r_J  =  min ( 1.0_KDR,  &
                    f  *  T_F ( iV )  *  S_F ( iV )  &
                       *  M_F ( iV )  *  N_F ( iV )  /  AMU  &
                    /  abs ( Delta_J ) )

        ! r_N  =  min ( 1.0_KDR,  &
        !             f  *  T_F ( iV )  *  S_F ( iV )  &
        !                *  M_F ( iV )  *  N_F ( iV )  /  AMU  &
        !             /  ( ( Mu_e ( iV )  -  Mu_n_p ( iV ) ) &
        !                  * abs ( Delta_N ) ) )
        r_N  =  1.0_KDR

        dJ ( iV )  =  r_J  *  Delta_J
        dN ( iV )  =  r_N  *  Delta_N

if ( r_J < 1.0_KDR .or. r_N < 1.0_KDR ) then
call Show ( '>>> Regulating approach to beta equilibrium', CONSOLE % ERROR )
call Show ( NM % Name, '>>> Species', CONSOLE % ERROR )
call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
call Show ( iV, '>>> iV', CONSOLE % ERROR )
call Show ( r_J, '>>> RegulationFactor_J', CONSOLE % ERROR )
call Show ( r_N, '>>> RegulationFactor_N', CONSOLE % ERROR )
end if

      end if
    end do
    !$OMP end parallel do

  end subroutine ImposeBetaEquilibrium_Kernel


  subroutine ComputeFluidSource_G_S_Radiation_Kernel &
               ( FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, FS_R_DS, IsProperCell, &
                 Xi_J, Chi_J, Chi_H, J, dE, H_1, H_2, H_3, dH_1, dH_2, dH_3, &
                 Beta_EQ, V_1, V_2, V_3, T, M_DD_22, M_DD_33, c, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3, &
      FS_R_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, Chi_J, &
      Chi_H, &
      J,  &
      dE, &
      H_1, H_2, H_3, &
      dH_1, dH_2, dH_3, &
      Beta_EQ, &
      V_1, V_2, V_3, &
      T, &
      M_DD_22, M_DD_33
    real ( KDR ) :: &
      c, &
      dT

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AMU

    nV = size ( FS_R_G )

    AMU = CONSTANT % ATOMIC_MASS_UNIT

! call Show ( FS_R_G, '>>> FS_R_G' )
! call Show ( E, '>>> E' )
! call Show ( c * dT, '>>> c * dT' )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
!      if ( Beta_EQ ( iV ) > 0.0_KDR ) then
!        FS_R_G ( iV )  =  - dJ ( iV )
!      else

        FS_R_G ( iV )  &
          =  FS_R_G ( iV )  &
             -  c * dT  &
                *  ( Xi_J ( iV )  &
                     -  Chi_J ( iV ) * ( J ( iV ) + dE ( iV ) )  &
                     -  Chi_H ( iV )  &
                        *  (                     V_1 ( iV ) * H_1 ( iV )  &
                             +  M_DD_22 ( iV ) * V_2 ( iV ) * H_2 ( iV )  &
                             +  M_DD_33 ( iV ) * V_3 ( iV ) * H_3 ( iV ) ) )
                     
!      end if

      FS_R_S_1 ( iV )  &
        =  FS_R_S_1 ( iV )  &
           -  c * dT  &  
              * ( -  Chi_H ( iV )  *  H_1 ( iV )  &
                  +  V_1 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      FS_R_S_2 ( iV )  &
        =  FS_R_S_2 ( iV )  &
           -  c * dT  *  M_DD_22 ( iV ) &  
              * ( -  Chi_H ( iV )  *  H_2 ( iV )  &
                  +  V_2 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      FS_R_S_3 ( iV )  &
        =  FS_R_S_3 ( iV )  &
           -  c * dT  *  M_DD_33 ( iV ) &  
              * ( -  Chi_H ( iV )  *  H_3 ( iV )  &
                  +  V_3 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      FS_R_DS ( iV )  &
        =  FS_R_DS ( iV )  &
           -  c * dT * AMU  /  T ( iV ) &
              *  ( Xi_J ( iV )  -  Chi_J ( iV ) * ( J ( iV ) + dE ( iV ) ) )

    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_G_S_Radiation_Kernel


  subroutine ComputeFluidSource_DP_Radiation_Kernel &
               ( FS_R_DP, FS_R_DS, IsProperCell, Xi_N, Chi_N, N, dD, Beta_EQ, &
                 T, Mu_n_p, Mu_e, c, dT, Sign )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_DP, &
      FS_R_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_N, Chi_N, &
      N, dD, &
      Beta_EQ, &
      T, Mu_n_p, Mu_e
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
!      if ( Beta_EQ ( iV ) > 0.0_KDR ) then
!        FS_R_DP ( iV )  =  - AMU * dJN ( iV )
!      else

        FS_R_DP ( iV )  &
          =  FS_R_DP ( iV )  &
             -  Sign * c * dT * AMU &
                *  ( Xi_N ( iV )  -  Chi_N ( iV ) * ( N ( iV ) + dD ( iV ) ) )

        FS_R_DS ( iV )  &
          =  FS_R_DS ( iV )  &
             +  Sign * c * dT * AMU &
                *  ( Mu_e ( iV )  -  Mu_n_p ( iV ) )  /  T ( iV )  &
                *  ( Xi_N ( iV )  -  Chi_N ( iV ) * ( N ( iV ) + dD ( iV ) ) )

!      end if
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_DP_Radiation_Kernel


end module WoosleyHeger_07_G__Form
