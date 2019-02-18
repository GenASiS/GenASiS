module RadiationBox_S_OS__Form

  !-- RadiationBox_Spectral_OperatorSplit_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_1D_MS_C_PS_Template ) :: &
    RadiationBox_S_OS_Form
      class ( Interactions_BSLL_ASC_CSLD_Template ), allocatable :: &
        Interactions_BSLL_ASC_CSLD
      class ( Relaxation_RM_BSLL_ASC_CSLD_Form ), allocatable :: &
        Relaxation_RM_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RadiationBox_S_OS_Form

contains


  subroutine Initialize &
               ( RB, RadiationName, RadiationType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption, MinCoordinateOption, MaxCoordinateOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 EnergyScaleOption, nCellsPositionOption, nCellsEnergyOption, &
                 nWriteOption )

    class ( RadiationBox_S_OS_Form ), intent ( inout ) :: &
      RB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in )  :: &
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

    integer ( KDI ) :: &
      iC, &  !-- iCurrents
      nCellsEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition
    real ( KDR ) :: &
      EnergyScale
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions, &
      EvolveFluid


    if ( RB % Type == '' ) &
      RB % Type = 'a RadiationBox_S_OS'


    !-- PositionSpace

    allocate ( Atlas_SC_Form :: RB % PositionSpace )
    select type ( PS => RB % PositionSpace )
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


    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: RB % MomentumSpace )
    select type ( MS => RB % MomentumSpace )
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


    !-- Prepare for Currents

    RB % N_CURRENTS_MS = size ( RadiationName )
    allocate ( RB % Current_BSLL_ASC_CSLD_1D ( RB % N_CURRENTS_MS ) )
    allocate ( RB % TimeStepLabel ( RB % N_CURRENTS_MS  +  1 ) )
    do iC = 1, RB % N_CURRENTS_MS
      RB % TimeStepLabel ( iC )  =  RadiationName ( iC )
    end do !-- iC
    RB % TimeStepLabel ( RB % N_CURRENTS_MS  +  1 )  =  'Fluid'
    

    !-- Radiation

    do iC = 1, RB % N_CURRENTS_MS
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )
        allocate &
          ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
              RB % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        select type ( RMB => RB % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
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
      case default
        call Show ( 'RadiationType not implemented', CONSOLE % ERROR )
        call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
        call Show ( 'RadiationBox_S_OS__Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end do !-- iC


    !-- Fluid

    allocate ( Fluid_ASC_Form :: RB % Current_ASC )
    select type ( FA => RB % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'IDEAL' )
    end select !-- FA


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
      allocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
      associate ( R => RB % Relaxation_RM_BSLL_ASC_CSLD )
      select type ( RMB => RB % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
      class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
        call R % Initialize ( RMB, Name = RB % Name )
      end select !-- RMA
      end associate !-- R
    end if !-- ApplyInteractions


    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form :: RB % Step_MS ) !-- Radiation
    select type ( S_MS => RB % Step_MS )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form )
    call S_MS % Initialize ( RB, RB % Current_BSLL_ASC_CSLD_1D, Name )

    if ( .not. ApplyStreaming ) then
      do iC = 1, size ( RadiationName )
        S_MS % ApplyDivergence_S ( iC ) % Pointer  =>  null ( )  
      end do !-- iC
    end if

    if ( ApplyInteractions ) then
      do iC = 1, size ( RadiationName )
        S_MS % ApplyRelaxation_F ( iC ) % Pointer  &
          =>  RB % Relaxation_RM_BSLL_ASC_CSLD % Apply 
      end do !-- iC
    end if

    end select !-- S_MS

    allocate ( Step_RK2_C_ASC_Form :: RB % Step_PS ) !-- Fluid
    select type ( S_PS => RB % Step_PS )
    class is ( Step_RK2_C_ASC_Form )
    call S_PS % Initialize ( RB, RB % Current_ASC, Name )
    if ( .not. EvolveFluid ) &
      S_PS % ApplyDivergence % Pointer => null ( )  !-- Disable fluid evolution
    end select !-- S_PS


    !-- Template

    call RB % InitializeTemplate_C_1D_MS_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )


    !-- Cleanup

    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( RB )

    type ( RadiationBox_S_OS_Form ), intent ( inout ) :: &
      RB

    if ( allocated ( RB % Relaxation_RM_BSLL_ASC_CSLD ) ) &
      deallocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
    if ( allocated ( RB % Interactions_BSLL_ASC_CSLD ) ) &
      deallocate ( RB % Interactions_BSLL_ASC_CSLD )

    call RB % FinalizeTemplate_C_1D_MS_C_PS ( )

  end subroutine Finalize

  
end module RadiationBox_S_OS__Form
