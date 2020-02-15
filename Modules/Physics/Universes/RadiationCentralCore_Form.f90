module RadiationCentralCore_Form

  use Basics
  use Mathematics
  use StressEnergies
  use ComputeGravity_Command
  use ApplyGravity_F__Command
  use FluidCentralCore_Form

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: RadiationCentralCoreForm
    real ( KDR ) :: &
      InteractionFactor
    character ( LDL ) :: &
      MomentsType = ''
    class ( Interactions_ASC_Template ), allocatable :: &
      Interactions_ASC
    class ( Interactions_BSLL_ASC_CSLD_Template ), allocatable :: &
      Interactions_BSLL_ASC_CSLD
    class ( Relaxation_RM_ASC_Form ), allocatable :: &
      Relaxation_RM_ASC
    class ( Relaxation_RM_BSLL_ASC_CSLD_Form ), allocatable :: &
      Relaxation_RM_BSLL_ASC_CSLD
  contains
    procedure, private, pass :: &
      Initialize_RCC
    generic, public :: &
      Initialize => Initialize_RCC
    final :: &
      Finalize
    procedure, private, pass :: &
      AllocateIntegrator_RCC
    generic, public :: &
      AllocateIntegrator => AllocateIntegrator_RCC
    procedure, public, pass :: &
      InitializeMomentumSpace
    procedure, public, pass :: &
      InitializeRadiation
    procedure, public, pass :: &
      InitializeSteps
  end type RadiationCentralCoreForm

    private :: &
      ComputeTimeStepLocal, &
      PrepareStep, &
      IntegrateSources, &
      ApplySources_Radiation, &
      ApplySources_Fluid

      private :: &
        ComputeTimeStepInteractions, &
        ComputeFluidSource_G_S_Radiation_Kernel, &
        ComputeFluidSource_DE_Radiation_Kernel, &
        ApplySources_Fluid_Kernel

contains


  subroutine Initialize_RCC &
               ( RCC, RadiationName, RadiationType, MomentsType, FluidType, &
                 GeometryType, Name, EnergySpacingOption, FinishTimeOption, &
                 CourantFactorOption, GravityFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 MinEnergyOption, MaxEnergyOption, MinWidthEnergyOption, &
                 EnergyScaleOption, RadiusMaxOption, RadiusCoreOption, &
                 RadialRatioOption, nCellsEnergyOption, nCellsPolarOption, &
                 nWriteOption )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      MomentsType, &
      FluidType, &
      GeometryType, &
      Name
    character ( * ), intent ( in ), optional :: &
      EnergySpacingOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      GravityFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      MinEnergyOption, &
      MaxEnergyOption, &
      MinWidthEnergyOption, &
      EnergyScaleOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsEnergyOption, &
      nCellsPolarOption, &
      nWriteOption

    real ( KDR ) :: &
      FinishTime

    if ( RCC % Type == '' ) &
      RCC % Type = 'a RadiationCentralCore'
    
    call RCC % InitializeTemplate ( Name )

    RCC % MomentsType  =  MomentsType

    call RCC % AllocateIntegrator &
           ( RadiationName )
    call RCC % InitializePositionSpace &
           ( GeometryType, RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption )
    call RCC % InitializeMomentumSpace &
           ( EnergySpacingOption = EnergySpacingOption, &
             MinEnergyOption = MinEnergyOption, &
             MaxEnergyOption = MaxEnergyOption, &
             MinWidthEnergyOption = MinWidthEnergyOption, &
             EnergyScaleOption = EnergyScaleOption, &
             nCellsEnergyOption = nCellsEnergyOption )
    call RCC % InitializeFluid &
           ( FluidType, LimiterParameterOption = LimiterParameterOption, &
             ShockThresholdOption = ShockThresholdOption )
    call RCC % InitializeRadiation &
           ( RadiationName, RadiationType )
    call RCC % InitializeSteps &
           ( Name )
    call RCC % InitializeDiagnostics &
           ( )

    call RCC % SetCentralTemplate ( )
    call RCC % SetCentralCore ( GeometryType, GravityFactorOption )

    FinishTime  =  1.0_KDR  *  RCC % Units % Time
    if ( present ( FinishTimeOption ) ) &
      FinishTime = FinishTimeOption

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_PS_Form )
      I % ComputeTimeStepLocal => ComputeTimeStepLocal
      call I % Initialize &
             ( RCC, Name, TimeUnitOption = RCC % Units % Time, &
               FinishTimeOption = FinishTime, &
               CourantFactorOption = CourantFactorOption, &
               nWriteOption = nWriteOption )
    end select !-- I

    call Show ( 'RadiationCentralCore parameter', RCC % IGNORABILITY )
    RCC % InteractionFactor  =  1.0e-2_KDR
    call PROGRAM_HEADER % GetParameter &
           ( RCC % InteractionFactor, 'InteractionFactor' )
    call Show ( RCC % InteractionFactor, 'InteractionFactor', &
                RCC % IGNORABILITY )

  end subroutine Initialize_RCC


  impure elemental subroutine Finalize ( RCC )

    type ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC

    if ( allocated ( RCC % Relaxation_RM_BSLL_ASC_CSLD ) ) &
      deallocate ( RCC % Relaxation_RM_BSLL_ASC_CSLD )
    if ( allocated ( RCC % Relaxation_RM_ASC ) ) &
      deallocate ( RCC % Relaxation_RM_ASC )

    if ( allocated ( RCC % Interactions_BSLL_ASC_CSLD ) ) &
      deallocate ( RCC % Interactions_BSLL_ASC_CSLD )
    if ( allocated ( RCC % Interactions_ASC ) ) &
      deallocate ( RCC % Interactions_ASC )

  end subroutine Finalize


  subroutine AllocateIntegrator_RCC ( RCC, RadiationName )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName

    integer ( KDI ) :: &
      iC  !-- iCurrent

    select case ( trim ( RCC % MomentsType ) )
    case ( 'NONE' )
      allocate ( Integrator_C_PS_Form :: RCC % Integrator )
    case ( 'GREY' )
      allocate ( Integrator_C_1D_PS_C_PS_Form :: RCC % Integrator )
    case ( 'SPECTRAL' )
      allocate ( Integrator_C_1D_MS_C_PS_Form :: RCC % Integrator )
    case default
      call Show ( 'MomentsType not recognized', CONSOLE % ERROR )
      call Show ( RCC % MomentsType, 'MomentsType', CONSOLE % ERROR )
      call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator_RCC', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- MomentsType

    select type ( I => RCC % Integrator )
    type is ( Integrator_C_PS_Form )
      allocate ( I % TimeStepLabel ( 1 + 1 ) )
    class is ( Integrator_C_1D_C_PS_Template )
      I % N_CURRENTS_1D  =  size ( RadiationName )
      allocate ( I % TimeStepLabel &
                   ( 1  +  1  +  I % N_CURRENTS_1D  +  I % N_CURRENTS_1D ) )
    end select !-- I

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_PS_Form )
      I % TimeStepLabel ( 1 )      =  'Gravity'
      I % TimeStepLabel ( 1 + 1 )  =  'Fluid Advection'
    end select !-- I

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )
      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( 1 + 1 + iC )  &
          =  trim ( RadiationName ( iC ) ) // ' Streaming'
      end do !-- iC
      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( 1  +  1  +  I % N_CURRENTS_1D  +  iC )  &
          =  trim ( RadiationName ( iC ) ) // ' Interactions'
      end do !-- iC
    end select !-- I

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )
      allocate ( I % Current_ASC_1D ( I % N_CURRENTS_1D ) )
    class is ( Integrator_C_1D_MS_C_PS_Form )
      allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENTS_1D ) )
    end select !-- I

  end subroutine AllocateIntegrator_RCC


  subroutine InitializeMomentumSpace &
               ( RCC, EnergySpacingOption, MinEnergyOption, MaxEnergyOption, &
                 MinWidthEnergyOption, EnergyScaleOption, nCellsEnergyOption )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), intent ( in ), optional :: &
      EnergySpacingOption
    real ( KDR ), intent ( in ), optional :: &
      MinEnergyOption, &
      MaxEnergyOption, &
      MinWidthEnergyOption, &
      EnergyScaleOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsEnergyOption

    integer ( KDI ) :: &
      nCellsEnergy
    real ( KDR ) :: &
      MinEnergy, &
      MaxEnergy, &
      MinWidthEnergy, &
      EnergyScale
    character ( LDL ) :: &
      EnergySpacing

    if ( .not. RCC % Dimensionless ) then
      RCC % Units % Coordinate_MS        =  UNIT % IDENTITY
      RCC % Units % Coordinate_MS ( 1 )  =  UNIT % MEGA_ELECTRON_VOLT
    end if

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_MS_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Bundle_SLL_ASC_CSLD_Form :: I % MomentumSpace )
    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, 'MomentumSpace' )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    nCellsEnergy = 16
    if ( present ( nCellsEnergyOption ) ) &
      nCellsEnergy = nCellsEnergyOption
    call PROGRAM_HEADER % GetParameter ( nCellsEnergy, 'nCellsEnergy' )

    EnergySpacing = 'GEOMETRIC'
    if ( present ( EnergySpacingOption ) ) &
      EnergySpacing = EnergySpacingOption
    call PROGRAM_HEADER % GetParameter ( EnergySpacing, 'EnergySpacing' )

    select case ( trim ( EnergySpacing ) )
    case ( 'GEOMETRIC' )

      MinEnergy       =    0.0_KDR
      MaxEnergy       =  100.0_KDR
      MinWidthEnergy  =    0.1_KDR
      if ( present ( MinEnergyOption ) ) &
        MinEnergy = MinEnergyOption
      if ( present ( MaxEnergyOption ) ) &
        MaxEnergy = MaxEnergyOption
      if ( present ( MinWidthEnergyOption ) ) &
        MinWidthEnergy = MinWidthEnergyOption
      call PROGRAM_HEADER % GetParameter ( MinEnergy, 'MinEnergy' )
      call PROGRAM_HEADER % GetParameter ( MaxEnergy, 'MaxEnergy' )
      call PROGRAM_HEADER % GetParameter ( MinWidthEnergy, 'MinWidthEnergy' )

      call MS % CreateChart &
             ( SpacingOption = [ 'GEOMETRIC' ], &
               CoordinateSystemOption = 'SPHERICAL', &
               CoordinateUnitOption = RCC % Units % Coordinate_MS, &
               MinCoordinateOption = [ MinEnergy ], &
               MaxCoordinateOption = [ MaxEnergy ], &
               ScaleOption = [ MinWidthEnergy ], &
               nCellsOption = [ nCellsEnergy ], &
               nGhostLayersOption = [ 0, 0, 0 ] )

    case ( 'COMPACTIFIED' )

      EnergyScale = 10.0_KDR
      if ( present ( EnergyScaleOption ) ) &
        EnergyScale = EnergyScaleOption
      call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

      call MS % CreateChart &
             ( SpacingOption = [ 'COMPACTIFIED' ], &
               CoordinateSystemOption = 'SPHERICAL', &
               CoordinateUnitOption = RCC % Units % Coordinate_MS, &
               ScaleOption = [ EnergyScale ], &
               nCellsOption = [ nCellsEnergy ], &
               nGhostLayersOption = [ 0, 0, 0 ] )

    case default

      call Show ( 'EnergySpacing not recognized', CONSOLE % ERROR )
      call Show ( EnergySpacing, 'EnergySpacing', CONSOLE % ERROR )
      call Show ( 'RadiationBox_Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeMomentumSpace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )

    end select !-- Spacing

    end select !-- MS
    end select !-- PS
    end select !-- I

  end subroutine InitializeMomentumSpace


  subroutine InitializeRadiation ( RCC, RadiationName, RadiationType )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType

    integer ( KDI ) :: &
      iC  !-- iCurrent

    select type ( I_1D => RCC % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    do iC = 1, I_1D % N_CURRENTS_1D

      select type ( I => I_1D )
      class is ( Integrator_C_1D_PS_C_PS_Form )

        select type ( PS => I % PositionSpace )
        class is ( Atlas_SC_Form )

        select case ( trim ( RadiationType ( iC ) ) )
        case ( 'NEUTRINOS_E', 'NEUTRINOS_E_BAR', 'NEUTRINOS_X' )
          allocate &
            ( NeutrinoMoments_ASC_Form :: I % Current_ASC_1D ( iC ) % Element )
        case default
          call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
          call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
          call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
          call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- RadiationType
       
        select type ( RA => I % Current_ASC_1D ( iC ) % Element )
        class is ( RadiationMoments_ASC_Form )
          call RA % Initialize &
                 ( PS, RadiationType = RadiationType ( iC ), &
                   MomentsType = 'GREY', Units = RCC % Units, &
                   NameShortOption = RadiationName ( iC ) )
        end select !-- RA        
        end select !-- PS

      class is ( Integrator_C_1D_MS_C_PS_Form )

        select type ( MS => I % MomentumSpace )
        class is ( Bundle_SLL_ASC_CSLD_Form )

        select case ( trim ( RadiationType ( iC ) ) )
        case ( 'NEUTRINOS_E', 'NEUTRINOS_E_BAR', 'NEUTRINOS_X' )
          allocate &
            ( NeutrinoMoments_BSLL_ASC_CSLD_Form :: &
                I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        case default
          call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
          call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
          call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
          call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- RadiationType

        select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
          call RMB % Initialize &
                 ( MS, RadiationType = RadiationType ( iC ), &
                   Units = RCC % Units, NameShortOption = RadiationName ( iC ) )
        end select !-- RMB
        end select !-- MS

      end select !--I

    end do !-- iC

    end select !-- I_1D

  end subroutine InitializeRadiation


  subroutine InitializeSteps ( RCC, Name )

    class ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iC  !-- iCurrent


    !-- Relaxation

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )

      allocate ( RCC % Relaxation_RM_ASC )
      associate ( R => RCC % Relaxation_RM_ASC )
      select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )
        call R % Initialize ( RMA, Name = RCC % Name )
      end select !-- RMA
      end associate !-- R

    class is ( Integrator_C_1D_MS_C_PS_Form )

      allocate ( RCC % Relaxation_RM_BSLL_ASC_CSLD )
      associate ( R => RCC % Relaxation_RM_BSLL_ASC_CSLD )
      select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
        call R % Initialize ( RMB, Name = RCC % Name )
      end select !-- RMB
      end associate !-- R

    end select !-- I


    !-- Radiation step

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      allocate ( Step_RK2_C_ASC_1D_Form :: I % Step_1D )
      select type ( S_1D => I % Step_1D )
      class is ( Step_RK2_C_ASC_1D_Form )
      call S_1D % Initialize ( I, I % Current_ASC_1D, Name )

      do iC = 1, I % N_CURRENTS_1D
        S_1D % ApplySources_1D ( iC ) % Pointer  &
          =>  ApplySources_Radiation
        S_1D % ApplyRelaxation_1D ( iC ) % Pointer  &
          =>  RCC % Relaxation_RM_ASC % Apply 
      end do !-- iC

      end select !-- S_1D

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      allocate ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form :: I % Step_1D )
      select type ( S_1D => I % Step_1D )
      class is ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form )
      call S_1D % Initialize ( I, I % Current_BSLL_ASC_CSLD_1D, Name )

      do iC = 1, I % N_CURRENTS_1D
        S_1D % ApplySources_S ( iC ) % Pointer  &
          =>  ApplySources_Radiation
        S_1D % ApplyRelaxation_F ( iC ) % Pointer  &
          =>  RCC % Relaxation_RM_BSLL_ASC_CSLD % Apply 
      end do !-- iC

      end select !-- S_1D

      I % PrepareStep => PrepareStep

    end select !-- I


    !-- Fluid

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_PS_Form )

      allocate ( Step_RK2_C_ASC_Form :: I % Step )
      select type ( S => I % Step )
      class is ( Step_RK2_C_ASC_Form )
        call S % Initialize ( I, I % Current_ASC, 'Fluid' )
        S % ComputeConstraints % Pointer  =>  ComputeGravity
        S % ApplySources % Pointer  =>  ApplySources_Fluid
      end select !-- S

    end select !-- I

  end subroutine InitializeSteps


  subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

    class ( IntegratorTemplate ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      oC  !-- oCandidate

    TimeStepCandidate  =  huge ( 1.0_KDR )

    !-- Gravity
    select type ( RCC => I % Universe )
    class is ( RadiationCentralCoreForm )
      call RCC % ComputeTimeStep_G_ASC ( TimeStepCandidate ( 1 ) )
    end select !-- I

    !-- Fluid advection, no radiation
    select type ( I )
    type is ( Integrator_C_PS_Form )
      call I % ComputeTimeStep_C_ASC &
             ( TimeStepCandidate ( 1 + 1 ), I % Current_ASC )
    end select !-- I

    select type ( I )
    class is ( Integrator_C_1D_C_PS_Template )
      !-- Fluid advection and Radiation streaming
      call I % ComputeTimeStepLocalTemplate &
             ( TimeStepCandidate ( 1 + 1  :  1  +  1  +  I % N_CURRENTS_1D ) )
      !-- Interactions
      oC  =  1  +  1  +  I % N_CURRENTS_1D
      call ComputeTimeStepInteractions &
             ( I, TimeStepCandidate ( oC + 1  :  oC  +  I % N_CURRENTS_1D ) )
    end select !-- I

  end subroutine ComputeTimeStepLocal


  subroutine PrepareStep ( I )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iC  !-- iCurrent
    class ( Fluid_P_HN_Form ), pointer :: &
      F
    class ( RadiationMomentsForm ), pointer :: &
      RM

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
      F => FA % Fluid_P_HN ( )
    end select !-- FA

    select type ( SF => F % Sources )
    class is ( Sources_F_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( PSC => PS % Chart )
    class is ( Chart_SL_Template )

    call Clear ( SF % Value )

    do iC = 1, I % N_CURRENTS_1D

      select type ( I )
      class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

        select type ( RMA => I % Current_ASC_1D ( iC ) % Element )
        class is ( RadiationMoments_ASC_Form )
          RM => RMA % RadiationMoments ( )
        end select !-- RMA

      class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

        call IntegrateSources ( I, iC )

        select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
        select type ( RMA => RMB % BundleIntegral )
        class is ( RadiationMoments_ASC_Form )
          RM => RMA % RadiationMoments ( )
        end select !-- RMA
        end select !-- RMB

      end select !-- I

      select type ( SR => RM % Sources )
      class is ( Sources_RM_Form )

      call ComputeFluidSource_G_S_Radiation_Kernel &
             ( SF % Value ( :, SF % RADIATION_G ), & 
               SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
               SF % Value ( :, SF % RADIATION_S_D ( 2 ) ), &
               SF % Value ( :, SF % RADIATION_S_D ( 3 ) ), &
               PSC % IsProperCell, &
               SR % Value ( :, SR % INTERACTIONS_J ), &
               SR % Value ( :, SR % INTERACTIONS_H_D ( 1 ) ), &
               SR % Value ( :, SR % INTERACTIONS_H_D ( 2 ) ), &
               SR % Value ( :, SR % INTERACTIONS_H_D ( 3 ) ) )

      select case ( trim ( RM % RadiationType ) )
      case ( 'NEUTRINOS_E' )
        call ComputeFluidSource_DE_Radiation_Kernel &
               ( SF % Value ( :, SF % RADIATION_DE ), & 
                 SF % Value ( :, SF % RADIATION_DS ), &
                 PSC % IsProperCell, &
                 SR % Value ( :, SR % INTERACTIONS_N ), &
                 F % Value ( :, F % TEMPERATURE ), &
                 F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
                 F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
                 Sign  =  + 1.0_KDR )
      case ( 'NEUTRINOS_E_BAR' )
        call ComputeFluidSource_DE_Radiation_Kernel &
               ( SF % Value ( :, SF % RADIATION_DE ), & 
                 SF % Value ( :, SF % RADIATION_DS ), &
                 PSC % IsProperCell, &
                 SR % Value ( :, SR % INTERACTIONS_N ), &
                 F % Value ( :, F % TEMPERATURE ), &
                 F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
                 F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
                 Sign  =  - 1.0_KDR )
      end select !-- RadiationType

      end select !-- SR

    end do !-- iC

    end select !-- PSC
    end select !-- PS
    end select !-- SF
    nullify ( F, RM )

  end subroutine PrepareStep


  subroutine IntegrateSources ( I, iC )

    type ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ), target :: &
      I
    integer ( KDI ), intent ( in ) :: &
      iC  !-- iCurrent

    integer ( KDI ) :: &
      iI, &  !-- iIntegral
      iF, &  !-- iFiber
      nIntegrals
    real ( KDR ), dimension ( : ), allocatable :: &
      Integral
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      Integrand
    type ( VolumeIntegralForm ) :: &
      VI
    class ( RadiationMomentsForm ), pointer :: &
      RMEI, &
      RMF

    select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
    class is ( RadiationMoments_BSLL_ASC_CSLD_Form )

    select type ( RMA => RMB % BundleIntegral )
    class is ( RadiationMoments_ASC_Form )

    RMEI => RMA % RadiationMoments ( )

    nIntegrals = 5

    allocate ( Integral  ( nIntegrals ) )
    allocate ( Integrand ( nIntegrals ) )
    do iI = 1, nIntegrals
      call Integrand ( iI ) % Initialize ( RMB % nEnergyValues )
    end do !-- iC

    associate ( MS => RMB % Bundle_SLL_ASC_CSLD )
    do iF = 1, MS % nFibers
      associate &
        ( iBC => MS % iaBaseCell ( iF ), &
          CF  => MS % Fiber_CSLL )

      RMF => RMB % RadiationMoments ( iF )
      select type ( SF => RMF % Sources )
      class is ( Sources_RM_Form )

        Integrand ( 1 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_J )
        Integrand ( 2 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_H_D ( 1 ) )
        Integrand ( 3 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_H_D ( 2 ) )
        Integrand ( 4 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_H_D ( 3 ) )

        SF % Value ( :, SF % INTERACTIONS_N ) &
          =  SF % Value ( :, SF % INTERACTIONS_J )  /  RMB % Energy
        Integrand ( 5 ) % Value  &
          =  SF % Value ( :, SF % INTERACTIONS_N )

      end select !-- SF

      call VI % Compute ( CF, Integrand, Integral )

      select type ( SEI => RMEI % Sources )
      class is ( Sources_RM_Form )
        SEI % Value ( iBC, SEI % INTERACTIONS_J )         = Integral ( 1 ) 
        SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 1 ) ) = Integral ( 2 ) 
        SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 2 ) ) = Integral ( 3 ) 
        SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 3 ) ) = Integral ( 4 ) 
        SEI % Value ( iBC, SEI % INTERACTIONS_N )         = Integral ( 5 ) 
      end select !-- SEI
      
      end associate !-- iBC, etc.
    end do !-- iF
    end associate !-- MS

    end select !-- RMA
    end select !-- RMB

    nullify ( RMEI, RMF )
 
  end subroutine IntegrateSources


  subroutine ApplySources_Radiation &
               ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( StorageForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    call ApplyCurvilinear_RM &
           ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

  end subroutine ApplySources_Radiation


  subroutine ApplySources_Fluid &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( StorageForm ), intent ( inout ), target :: &
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
      iEntropy, &
      iElectron

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyGravity_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    select type ( I => S % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    select type ( SF => Sources_F )
    class is ( Sources_F_Form )

    select type ( F => Fluid )
    class is ( Fluid_P_HN_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 3 ), iMomentum_3 )
    call Search ( F % iaConserved, F % CONSERVED_ENTROPY, iEntropy )
    call Search ( F % iaConserved, F % CONSERVED_ELECTRON_DENSITY, iElectron )

    call ApplySources_Fluid_Kernel &
           ( Increment % Value ( :, iEnergy ), &
             Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Increment % Value ( :, iMomentum_3 ), &
             Increment % Value ( :, iElectron ), &
             Increment % Value ( :, iEntropy ), &
             Chart % IsProperCell, &
             SF % Value ( :, SF % RADIATION_G ), &
             SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
             SF % Value ( :, SF % RADIATION_S_D ( 2 ) ), &
             SF % Value ( :, SF % RADIATION_S_D ( 3 ) ), &
             SF % Value ( :, SF % RADIATION_DE ), &
             SF % Value ( :, SF % RADIATION_DS ), &
             TimeStep )

    end select !-- Chart
    end select !-- F
    end select !-- SF
    end select !-- I

  end subroutine ApplySources_Fluid


  subroutine ComputeTimeStepInteractions ( I, TimeStepCandidate )

    class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      TimeStepCandidate

    integer ( KDI ) :: &
      iC     !-- iCurrent

    select type ( RCC => I % Universe )
    class is ( RadiationCentralCoreForm )

    select type ( I )
    class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

      if ( .not. allocated ( RCC % Interactions_ASC ) ) then
        call Show ( 'Interactions_ASC not allocated', CONSOLE % ERROR )
        call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeTimeStepInteractions', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      associate ( IA => RCC % Interactions_ASC )  
      do iC = 1, I % N_CURRENTS_1D  
        associate ( RMA => I % Current_ASC_1D ( iC ) % Element )
        call IA % ComputeTimeScale ( RMA, TimeStepCandidate ( iC ) )
        TimeStepCandidate ( iC )  &
          =  RCC % InteractionFactor  *  TimeStepCandidate ( iC )
        end associate !-- RA
      end do !-- iC
      end associate !-- IA

    class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

      if ( .not. allocated ( RCC % Interactions_BSLL_ASC_CSLD ) ) then
        call Show ( 'Interactions_BSLL_ASC_CSLD not allocated', &
                    CONSOLE % ERROR )
        call Show ( 'RadiationCentralCore_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeTimeStepInteractions', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if

      associate ( IB => RCC % Interactions_BSLL_ASC_CSLD )
      do iC = 1, I % N_CURRENTS_1D  
        associate ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        call IB % ComputeTimeScale ( RMB, TimeStepCandidate ( iC ) )
        TimeStepCandidate ( iC )  &
          =  RCC % InteractionFactor  *  TimeStepCandidate ( iC )
        end associate !-- RCC
      end do !-- iC
      end associate !-- IB

    end select !-- I
    end select !-- RCC

  end subroutine ComputeTimeStepInteractions


  subroutine ComputeFluidSource_G_S_Radiation_Kernel &
               ( FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, IsProperCell, &
                 Q, A_1, A_2, A_3 )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Q, &
      A_1, A_2, A_3

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( FS_R_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      FS_R_G ( iV )    =    FS_R_G ( iV )  -    Q ( iV ) 
      FS_R_S_1 ( iV )  =  FS_R_S_1 ( iV )  -  A_1 ( iV )
      FS_R_S_2 ( iV )  =  FS_R_S_2 ( iV )  -  A_2 ( iV )
      FS_R_S_3 ( iV )  =  FS_R_S_3 ( iV )  -  A_3 ( iV )
    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_G_S_Radiation_Kernel


  subroutine ComputeFluidSource_DE_Radiation_Kernel &
               ( FS_R_DE, FS_R_DS, IsProperCell, R, T, Mu_e, Mu_n_p, Sign )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      FS_R_DE, &
      FS_R_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      T, Mu_e, Mu_n_p
    real ( KDR ), intent ( in ) :: &
      Sign

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( FS_R_DE )

    !$OMP parallel do private ( iV )
    do iV = 1, nV

      if ( .not. IsProperCell ( iV ) ) &
        cycle

      FS_R_DE ( iV )  &
        =  FS_R_DE ( iV )  &
           -  Sign  *  R ( iV )

      FS_R_DS ( iV )  &
        =  FS_R_DS ( iV )  &
           +  Sign  *  ( Mu_e ( iV )  -  Mu_n_p ( iV ) )  /  T ( iV )  &
                    *  R ( iV )

    end do
    !$OMP end parallel do

  end subroutine ComputeFluidSource_DE_Radiation_Kernel


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, K_S_2, K_S_3, K_DE, K_DS, IsProperCell, FS_R_G, &
                 FS_R_S_1, FS_R_S_2, FS_R_S_3, FS_R_DE, FS_R_DS, dt )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_G, &
      K_S_1, K_S_2, K_S_3, &
      K_DE, &
      K_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3, &
      FS_R_DE, &
      FS_R_DS
    real ( KDR ), intent ( in ) :: &
      dt

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( K_G )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle
      K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )    *  dt
      K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )  *  dt
      K_S_2 ( iV )  =  K_S_2 ( iV )  +  FS_R_S_2 ( iV )  *  dt
      K_S_3 ( iV )  =  K_S_3 ( iV )  +  FS_R_S_3 ( iV )  *  dt
      K_DE ( iV )   =  K_DE ( iV )   +  FS_R_DE ( iV )   *  dt
      K_DS ( iV )   =  K_DS ( iV )   +  FS_R_DS ( iV )   *  dt
    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


end module RadiationCentralCore_Form
