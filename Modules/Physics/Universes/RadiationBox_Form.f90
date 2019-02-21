module RadiationBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use Universe_Template

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RadiationBoxForm
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
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RadiationBoxForm

    private :: &
      AllocateIntegrator, &
      InitializePositionSpace, &
      InitializeMomentumSpace, &
      InitializeRadiation, &
      InitializeFluid, &
      InitializeSteps

contains

  
  subroutine Initialize &
               ( RB, RadiationName, RadiationType, MomentsType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption, MinCoordinateOption, MaxCoordinateOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 EnergyScaleOption, nCellsPositionOption, nCellsEnergyOption, &
                 nWriteOption )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      MomentsType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption, &
      EvolveFluidOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      EnergyScaleOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsPositionOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsEnergyOption, &
      nWriteOption

    if ( RB % Type == '' ) &
      RB % Type = 'a RadiationBox'
    
    call RB % InitializeTemplate ( Name )

    RB % MomentsType  =  MomentsType

    call AllocateIntegrator &
           ( RB, RadiationName )
    call InitializePositionSpace &
           ( RB, MinCoordinateOption, MaxCoordinateOption, &
             nCellsPositionOption )
    call InitializeMomentumSpace &
           ( RB, EnergyScaleOption, nCellsEnergyOption )
    call InitializeRadiation &
           ( RB, RadiationName, RadiationType )
    call InitializeFluid &
           ( RB )
    call InitializeSteps &
           ( RB, Name, ApplyStreamingOption, ApplyInteractionsOption, &
             EvolveFluidOption )

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )
      call I % Initialize &
             ( Name, TimeUnitOption = TimeUnitOption, &
               FinishTimeOption = FinishTimeOption, &
               CourantFactorOption = CourantFactorOption, &
               nWriteOption = nWriteOption )
    end select !-- I

  end subroutine Initialize


  subroutine Finalize ( RB )

    type ( RadiationBoxForm ), intent ( inout ) :: &
      RB

    if ( allocated ( RB % Relaxation_RM_BSLL_ASC_CSLD ) ) &
      deallocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
    if ( allocated ( RB % Relaxation_RM_ASC ) ) &
      deallocate ( RB % Relaxation_RM_ASC )

    if ( allocated ( RB % Interactions_BSLL_ASC_CSLD ) ) &
      deallocate ( RB % Interactions_BSLL_ASC_CSLD )
    if ( allocated ( RB % Interactions_ASC ) ) &
      deallocate ( RB % Interactions_ASC )

    call RB % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine AllocateIntegrator ( RB, RadiationName )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName

    integer ( KDI ) :: &
      iC  !-- iCurrent

    select case ( trim ( RB % MomentsType ) )
    case ( 'GREY' )
      allocate ( Integrator_C_1D_PS_C_PS_Form :: RB % Integrator )
    case ( 'SPECTRAL' )
      allocate ( Integrator_C_1D_MS_C_PS_Form :: RB % Integrator )
    case default
      call Show ( 'MomentsType not recognized', CONSOLE % ERROR )
      call Show ( RB % MomentsType, 'MomentsType', CONSOLE % ERROR )
      call Show ( 'RadiationBox_Form', 'module', CONSOLE % ERROR )
      call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- MomentsType

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )
      I % N_CURRENTS_1D  =  size ( RadiationName )
      allocate ( I % TimeStepLabel ( I % N_CURRENTS_1D  +  1 ) )
      do iC = 1, I % N_CURRENTS_1D
        I % TimeStepLabel ( iC )  =  RadiationName ( iC )
      end do !-- iC
      I % TimeStepLabel ( I % N_CURRENTS_1D  +  1 )  =  'Fluid'
    end select !-- I

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )
      allocate ( I % Current_ASC_1D ( I % N_CURRENTS_1D ) )
    class is ( Integrator_C_1D_MS_C_PS_Form )
      allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENTS_1D ) )
    end select !-- I

  end subroutine AllocateIntegrator


  subroutine InitializePositionSpace &
               ( RB, MinCoordinateOption, MaxCoordinateOption, &
                 nCellsPositionOption )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsPositionOption

    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    allocate ( Atlas_SC_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    nCellsPosition = [ 32, 32, 32 ]
    if ( present ( nCellsPositionOption ) ) &
      nCellsPosition = nCellsPositionOption
    call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

    call PS % CreateChart &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsPosition )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize ( PS, GeometryType = 'GALILEAN' )
    call PS % SetGeometry ( GA )

    end select !-- GA
    end select !-- PS
    end select !-- I

  end subroutine InitializePositionSpace


  subroutine InitializeMomentumSpace &
               ( RB, EnergyScaleOption, nCellsEnergyOption )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    real ( KDR ), intent ( in ), optional :: &
      EnergyScaleOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsEnergyOption

    integer ( KDI ) :: &
      nCellsEnergy
    real ( KDR ) :: &
      EnergyScale

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_MS_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Bundle_SLL_ASC_CSLD_Form :: I % MomentumSpace )
    select type ( MS => I % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, 'MomentumSpace' )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    EnergyScale = 10.0_KDR
    if ( present ( EnergyScaleOption ) ) &
      EnergyScale = EnergyScaleOption
    call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

!    CoordinateUnit = UNIT % IDENTITY
!    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

    nCellsEnergy = 16
    if ( present ( nCellsEnergyOption ) ) &
      nCellsEnergy = nCellsEnergyOption
    call PROGRAM_HEADER % GetParameter ( nCellsEnergy, 'nCellsEnergy' )

    call MS % CreateChart &
           ( SpacingOption = [ 'COMPACTIFIED' ], &
             CoordinateSystemOption = 'SPHERICAL', &
!             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = [ EnergyScale ], &
             nCellsOption = [ nCellsEnergy ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    end select !-- MS
    end select !-- PS
    end select !-- I

  end subroutine InitializeMomentumSpace


  subroutine InitializeRadiation ( RB, RadiationName, RadiationType )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType

    integer ( KDI ) :: &
      iC  !-- iCurrent

    select type ( I_1D => RB % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    do iC = 1, I_1D % N_CURRENTS_1D
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )

        select type ( I => I_1D )
        class is ( Integrator_C_1D_PS_C_PS_Form )

          select type ( PS => I % PositionSpace )
          class is ( Atlas_SC_Form )

          allocate &
            ( RadiationMoments_ASC_Form :: &
                I % Current_ASC_1D ( iC ) % Element )
          select type ( RA => I % Current_ASC_1D ( iC ) % Element )
          type is ( RadiationMoments_ASC_Form )
          call RA % Initialize &
                 ( PS, RadiationType ( iC ), &
                   NameShortOption = RadiationName ( iC ) )
                   ! Velocity_U_UnitOption = WHH % VelocityUnit, &
                   ! MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
                   ! MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
                   ! EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
                   ! TemperatureUnitOption = WHH % TemperatureUnit )
          end select !-- RA        
          end select !-- PS

        class is ( Integrator_C_1D_MS_C_PS_Form )

          select type ( MS => I % MomentumSpace )
          class is ( Bundle_SLL_ASC_CSLD_Form )

          allocate &
            ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
                I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
          select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
          class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
          call RMB % Initialize &
                 ( MS, RadiationType ( iC ), &
                   NameShortOption = RadiationName ( iC ) )
                   ! Velocity_U_UnitOption = WHH % VelocityUnit, &
                   ! MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
                   ! MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
                   ! EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
                   ! TemperatureUnitOption = WHH % TemperatureUnit )
          end select !-- RMB
          end select !-- MS

        end select !--I

      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
        call Show ( 'RadiationBox_G_OS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end do !-- iC

    end select !-- I_1D

  end subroutine InitializeRadiation


  subroutine InitializeFluid ( RB )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_C_PS_Template )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( Fluid_ASC_Form :: I % Current_ASC )
    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'IDEAL' )
    end select !-- FA

    end select !-- PS
    end select !-- I

  end subroutine InitializeFluid


  subroutine InitializeSteps &
               ( RB, Name, ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption )

    class ( RadiationBoxForm ), intent ( inout ) :: &
      RB
    character ( * ), intent ( in ) :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption, &
      EvolveFluidOption

    integer ( KDI ) :: &
      iC  !-- iCurrent
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions, &
      EvolveFluid


    !-- Operators

    ApplyStreaming    = .true.
    ApplyInteractions = .true.
    EvolveFluid       = .true.
    if ( present ( ApplyStreamingOption ) ) &
      ApplyStreaming = ApplyStreamingOption
    if ( present ( ApplyInteractionsOption ) ) &
      ApplyInteractions = ApplyInteractionsOption
    if ( present ( EvolveFluidOption ) ) &
      EvolveFluid = EvolveFluidOption


    !-- Relaxation

    if ( ApplyInteractions ) then
      select type ( I => RB % Integrator )
      class is ( Integrator_C_1D_PS_C_PS_Form )

        allocate ( RB % Relaxation_RM_ASC )
        associate ( R => RB % Relaxation_RM_ASC )
        select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
        class is ( RadiationMoments_ASC_Form )
          call R % Initialize ( RMA, Name = RB % Name )
        end select !-- RMA
        end associate !-- R

      class is ( Integrator_C_1D_MS_C_PS_Form )

        allocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
        associate ( R => RB % Relaxation_RM_BSLL_ASC_CSLD )
        select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
          call R % Initialize ( RMB, Name = RB % Name )
        end select !-- RMB
        end associate !-- R

      end select !-- I
    end if !-- ApplyInteractions


    !-- Radiation step

    select type ( I => RB % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )

      allocate ( Step_RK2_C_ASC_1D_Form :: I % Step_1D )
      select type ( S_1D => I % Step_1D )
      class is ( Step_RK2_C_ASC_1D_Form )
      call S_1D % Initialize ( I, I % Current_ASC_1D, Name )

      if ( .not. ApplyStreaming ) then
        do iC = 1, I % N_CURRENTS_1D
          S_1D % ApplyDivergence_1D ( iC ) % Pointer  =>  null ( )  
        end do !-- iC
      end if

      if ( ApplyInteractions ) then
        do iC = 1, I % N_CURRENTS_1D
          S_1D % ApplyRelaxation_1D ( iC ) % Pointer  &
            =>  RB % Relaxation_RM_ASC % Apply 
        end do !-- iC
      end if

      end select !-- S_1D

    class is ( Integrator_C_1D_MS_C_PS_Form )

      allocate ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form :: I % Step_1D )
      select type ( S_1D => I % Step_1D )
      class is ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form )
      call S_1D % Initialize ( I, I % Current_BSLL_ASC_CSLD_1D, Name )

      if ( .not. ApplyStreaming ) then
        do iC = 1, I % N_CURRENTS_1D
          S_1D % ApplyDivergence_S ( iC ) % Pointer  =>  null ( )  
        end do !-- iC
      end if

      if ( ApplyInteractions ) then
        do iC = 1, I % N_CURRENTS_1D
          S_1D % ApplyRelaxation_F ( iC ) % Pointer  &
            =>  RB % Relaxation_RM_BSLL_ASC_CSLD % Apply 
        end do !-- iC
      end if

      end select !-- S_1D

    end select !-- I


  end subroutine InitializeSteps


end module RadiationBox_Form
