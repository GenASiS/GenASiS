module RadiationCentralCore_Form

  use Basics
  use Mathematics
  use StressEnergies
  use FluidCentralCore_Form

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: RadiationCentralCoreForm
    character ( LDL ) :: &
      MomentsType = ''
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
!      ComputeTimeStepLocal, &
!      PrepareStep, &
!      IntegrateSources, &
      ApplySources_Radiation, &
      ApplySources_Fluid

      private :: &
!        ComputeTimeStepInteractions, &
!        ComputeFluidSource_G_S_Radiation_Kernel, &
        ApplySources_Fluid_Kernel

contains


  subroutine Initialize_RCC &
               ( RCC, RadiationName, RadiationType, MomentsType, FluidType, &
                 GeometryType, Name, FinishTimeOption, CourantFactorOption, &
                 GravityFactorOption, LimiterParameterOption, &
                 ShockThresholdOption, MinEnergyOption, MaxEnergyOption, &
                 MinWidthEnergyOption, EnergyScaleOption, RadiusMaxOption, &
                 RadiusCoreOption, RadialRatioOption, nCellsEnergyOption, &
                 nCellsPolarOption, nWriteOption )

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
           ( MinEnergyOption = MinEnergyOption, &
             MaxEnergyOption = MaxEnergyOption, &
             MinWidthEnergyOption = MinWidthEnergyOption, &
             EnergyScaleOption = EnergyScaleOption, &
             nCellsEnergyOption = nCellsEnergyOption )
    call RCC % InitializeRadiation &
           ( RadiationName, RadiationType )
    call RCC % InitializeFluid &
           ( FluidType, LimiterParameterOption = LimiterParameterOption, &
             ShockThresholdOption = ShockThresholdOption )
    call RCC % InitializeSteps &
           ( Name )

  end subroutine Initialize_RCC


  impure elemental subroutine Finalize ( RCC )

    type ( RadiationCentralCoreForm ), intent ( inout ) :: &
      RCC

    if ( allocated ( RCC % Relaxation_RM_BSLL_ASC_CSLD ) ) &
      deallocate ( RCC % Relaxation_RM_BSLL_ASC_CSLD )
    if ( allocated ( RCC % Relaxation_RM_ASC ) ) &
      deallocate ( RCC % Relaxation_RM_ASC )

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
    class is ( Integrator_C_1D_C_PS_Template )

      I % N_CURRENTS_1D  =  size ( RadiationName )
      allocate ( I % TimeStepLabel &
                   ( I % N_CURRENTS_1D  +  1  +  1  +  I % N_CURRENTS_1D ) )

      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( iC )  &
          =  trim ( RadiationName ( iC ) ) // ' Streaming'
      end do !-- iC

      I % TimeStepLabel ( I % N_CURRENTS_1D  +  1 )  &
          =  'Fluid Advection'

      I % TimeStepLabel ( I % N_CURRENTS_1D  +  1  +  1 )  &
          =  'Gravity'

      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( I % N_CURRENTS_1D  +  1  +  iC )  &
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
    character ( LDL ) :: &
      RadiationTypeLocal

    select type ( I_1D => RCC % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    do iC = 1, I_1D % N_CURRENTS_1D

      select type ( I => I_1D )
      class is ( Integrator_C_1D_PS_C_PS_Form )

        select type ( PS => I % PositionSpace )
        class is ( Atlas_SC_Form )

        select case ( trim ( RadiationType ( iC ) ) )
        case ( 'NEUTRINOS' )
          allocate &
            ( NeutrinoMoments_ASC_Form :: I % Current_ASC_1D ( iC ) % Element )
          RadiationTypeLocal = 'NEUTRINOS_GREY'
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
                 ( PS, RadiationTypeLocal, RCC % Units, &
                   NameShortOption = RadiationName ( iC ) )
        end select !-- RA        
        end select !-- PS

      class is ( Integrator_C_1D_MS_C_PS_Form )

        select type ( MS => I % MomentumSpace )
        class is ( Bundle_SLL_ASC_CSLD_Form )

        select case ( trim ( RadiationType ( iC ) ) )
        case ( 'NEUTRINOS' )
          allocate &
            ( NeutrinoMoments_BSLL_ASC_CSLD_Form :: &
                I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
          RadiationTypeLocal = 'NEUTRINOS_SPECTRAL'
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
                 ( MS, RadiationTypeLocal, RCC % Units, &
                   NameShortOption = RadiationName ( iC ) )
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

    end select !-- I


    !-- Fluid

    select type ( I => RCC % Integrator )
    class is ( Integrator_C_PS_Form )

      allocate ( Step_RK2_C_ASC_Form :: I % Step )
      select type ( S => I % Step )
      class is ( Step_RK2_C_ASC_Form )
        call S % Initialize ( I, I % Current_ASC, 'Fluid' )
        S % ApplySources % Pointer  =>  ApplySources_Fluid
      end select !-- S

!    I % PrepareStep => PrepareStep

    end select !-- I

  end subroutine InitializeSteps


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
      iMomentum_3

    call ApplyCurvilinear_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    call ApplyGravity_F &
           ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    select type ( I => S % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    select type ( SF => Sources_F )
    class is ( Sources_F_Form )

    select type ( F => Fluid )
    class is ( Fluid_P_I_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
    call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 3 ), iMomentum_3 )

    call ApplySources_Fluid_Kernel &
           ( Increment % Value ( :, iEnergy ), &
             Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             Increment % Value ( :, iMomentum_3 ), &
             Chart % IsProperCell, &
             SF % Value ( :, SF % RADIATION_G ), &
             SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
             SF % Value ( :, SF % RADIATION_S_D ( 2 ) ), &
             SF % Value ( :, SF % RADIATION_S_D ( 3 ) ), &
             TimeStep )

    end select !-- Chart
    end select !-- F
    end select !-- SF
    end select !-- I

  end subroutine ApplySources_Fluid


  subroutine ApplySources_Fluid_Kernel &
               ( K_G, K_S_1, K_S_2, K_S_3, IsProperCell, &
                 FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, dt )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_G, &
      K_S_1, K_S_2, K_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FS_R_G, &
      FS_R_S_1, FS_R_S_2, FS_R_S_3
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
    end do
    !$OMP end parallel do

  end subroutine ApplySources_Fluid_Kernel


end module RadiationCentralCore_Form
