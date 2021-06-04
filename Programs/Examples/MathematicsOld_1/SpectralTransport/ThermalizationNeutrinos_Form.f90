module ThermalizationNeutrinos_Form

  use Basics
  use Mathematics

  use RadiationMoments_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( Integrator_C_MS_C_PS_Template ) :: &
    ThermalizationNeutrinosForm
!!$      type ( RadiationMoments_ASC_Form ), allocatable :: &
!!$        Reference_ASC, &
!!$        FractionalDifference_ASC
!!$      type ( Interactions_BSLL_ASC_CSLD_Form ), allocatable :: &
!!$        Interactions_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
!!$    procedure, private, pass :: &
!!$      ComputeTimeStepLocal
  end type ThermalizationNeutrinosForm

    real ( KDR ), private :: &
      MassDensityMin, &
      MassDensityMax, &
      TemperatureMin, &
      TemperatureMax, &
      ElectronFractionMin, &
      ElectronFractionMax, &
      EnergyScale!, &
!!$      Opacity, &
!!$      TimeScale
    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate


contains


  subroutine Initialize ( T, Name )

    class ( ThermalizationNeutrinosForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ), dimension ( 3 ) :: &
      nPositionCells      
    real ( KDR ) :: &
      MaxRadius
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit_PS

    integer ( KDI ) :: &
      nEnergyCells
    real ( KDR ), dimension ( 3 ) :: &
      Scale
    type ( MeasuredValueForm ) :: &
      EnergyDensityUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    MassDensityMin      =  1.0d08  *  UNIT % GRAM / UNIT % CENTIMETER ** 3
    MassDensityMax      =  4.0d14  *  UNIT % GRAM / UNIT % CENTIMETER ** 3
    TemperatureMin      =  5.0d09  *  UNIT % KELVIN
    TemperatureMax      =  2.6d11  *  UNIT % KELVIN
    ElectronFractionMin =  0.30_KDR
    ElectronFractionMax =  0.46_KDR
    EnergyScale         =  10.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT

    call PROGRAM_HEADER % GetParameter ( MassDensityMin, 'MassDensityMin' )
    call PROGRAM_HEADER % GetParameter ( MassDensityMax, 'MassDensityMax' )
    call PROGRAM_HEADER % GetParameter ( TemperatureMin, 'TemperatureMin' )
    call PROGRAM_HEADER % GetParameter ( TemperatureMax, 'TemperatureMax' )
    call PROGRAM_HEADER % GetParameter &
           ( ElectronFractionMin, 'ElectronFractionMin' )
    call PROGRAM_HEADER % GetParameter &
           ( ElectronFractionMax, 'ElectronFractionMax' )
    call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

!!$    Opacity = 1.0_KDR
!!$    call PROGRAM_HEADER % GetParameter ( Opacity, 'Opacity' )
!!$
!!$    associate &
!!$      ( c     => CONSTANT % SPEED_OF_LIGHT, &
!!$        Kappa => Opacity )
!!$    TimeScale  =  1.0 / ( c * Kappa )
!!$    end associate !-- c, etc.
    
    !-- PositionSpace

    allocate ( Atlas_SC_Form :: T % PositionSpace )
    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    MaxRadius = 100.0_KDR * UNIT % KILOMETER

    CoordinateUnit_PS  =  UNIT % KILOMETER

    MinCoordinate  =  - MaxRadius
    MaxCoordinate  =  + MaxRadius

    associate ( nD => PS % nDimensions )
    nPositionCells = 1
    nPositionCells ( : nD ) = 32
    call PROGRAM_HEADER % GetParameter &
           ( nPositionCells ( : nD ), 'nPositionCells' )
    end associate !-- nD

    call PS % CreateChart &
           ( CoordinateUnitOption = CoordinateUnit_PS, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nPositionCells )

    !-- Geometry of PositionSpace

    allocate ( T % Geometry_ASC )
    associate ( GA => T % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- MomentumSpace

    allocate ( Bundle_SLL_ASC_CSLD_Form :: T % MomentumSpace )
    select type ( MS => T % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, Name )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    Scale = 0.0_KDR
    Scale ( 1 ) = EnergyScale

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

    Spacing = ''
    Spacing ( 1 ) = 'COMPACTIFIED'

    nEnergyCells = 20
    call PROGRAM_HEADER % GetParameter ( nEnergyCells, 'nEnergyCells' )

    call MS % CreateChart &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = 'SPHERICAL', &
             CoordinateUnitOption = CoordinateUnit, &
             ScaleOption = Scale, &
             nCellsOption = [ nEnergyCells, 1, 1 ], &
             nGhostLayersOption = [ 0, 0, 0 ] )

    !-- Radiation

    EnergyDensityUnit  =  UNIT % MEGA_ELECTRON_VOLT / UNIT % HBAR_C ** 3

    allocate ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
               T % Current_BSLL_ASC_CSLD )
    select type ( RMB => T % Current_BSLL_ASC_CSLD )
    type is ( RadiationMoments_BSLL_ASC_CSLD_Form )
    call RMB % Initialize &
           ( MS, 'NEUTRINOS_SPECTRAL', &
             EnergyDensityUnitOption = EnergyDensityUnit )

!!$    !-- Fluid
!!$
!!$    allocate ( Fluid_ASC_Form :: T % Current_ASC )
!!$    select type ( FA => T % Current_ASC )
!!$    class is ( Fluid_ASC_Form )
!!$    call FA % Initialize &
!!$           ( PS, 'NON_RELATIVISTIC', &
!!$             TemperatureUnitOption = UNIT % MEGA_ELECTRON_VOLT )
!!$
!!$    !-- Interactions
!!$
!!$    allocate ( T % Interactions_BSLL_ASC_CSLD )
!!$    associate ( IB => T % Interactions_BSLL_ASC_CSLD )
!!$    call IB % Initialize &
!!$           ( MS, 'FIXED', EnergyDensityUnitOption = EnergyDensityUnit )
!!$    call RMB % SetInteractions ( IB )
!!$    end associate !-- IB
!!$
!!$    !-- Step
!!$
!!$    allocate ( Step_RK2_C_BSLL_ASC_CSLD_Form :: T % Step_MS )
!!$    select type ( S_MS => T % Step_MS )
!!$    class is ( Step_RK2_C_BSLL_ASC_CSLD_Form )
!!$    call S_MS % Initialize ( RMB, Name )
!!$    S_MS % ApplyDivergence_S % Pointer => null ( )  !-- Disable spatial 
!!$                                                    !   section evolution
!!$    S_MS % ApplyRelaxation_F % Pointer => ApplyRelaxation_RM_S
!!$    end select !-- S_MS
!!$
!!$    allocate ( Step_RK2_C_ASC_Form :: T % Step_PS )
!!$    select type ( S_PS => T % Step_PS )
!!$    class is ( Step_RK2_C_ASC_Form )
!!$    call S_PS % Initialize ( FA, Name )
!!$    S_PS % ApplyDivergence % Pointer => null ( )  !-- Disable fluid evolution
!!$    end select !-- S_PS
!!$
!!$    !-- Diagnostics
!!$
!!$    EnergyDensityUnit  =  UNIT % MEGA_ELECTRON_VOLT ** 4 / UNIT % HBAR_C ** 3
!!$
!!$    allocate ( T % Reference_ASC )
!!$    allocate ( T % FractionalDifference_ASC )
!!$    call T % Reference_ASC % Initialize &
!!$           ( PS, 'GENERIC', NameShortOption = 'Reference', &
!!$             EnergyDensityUnitOption = EnergyDensityUnit, &
!!$             AllocateSourcesOption = .false., &
!!$             IgnorabilityOption = CONSOLE % INFO_5 )
!!$    call T % FractionalDifference_ASC % Initialize &
!!$           ( PS, 'GENERIC', NameShortOption = 'FractionalDifference', &
!!$             AllocateSourcesOption = .false., &
!!$             IgnorabilityOption = CONSOLE % INFO_5 )
!!$    T % SetReference => SetReference
!!$
!!$    !-- Initial conditions
!!$
!!$    call SetFluid ( T )
!!$    call SetRadiation ( T )
!!$    call SetInteractions ( T )
!!$
!!$    !-- Initialize template
!!$
!!$    call T % InitializeTemplate_C_MS_C_PS &
!!$           ( Name, FinishTimeOption = 10.0_KDR * TimeScale )
!!$
!!$    !-- Cleanup
!!$
!!$    end select !-- FA
    end select !-- RMB
    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( T )

    type ( ThermalizationNeutrinosForm ), intent ( inout ) :: &
      T

!!$    if ( allocated ( T % Interactions_BSLL_ASC_CSLD ) ) &
!!$      deallocate ( T % Interactions_BSLL_ASC_CSLD )
!!$    if ( allocated ( T % FractionalDifference_ASC ) ) &
!!$      deallocate ( T % FractionalDifference_ASC )
!!$    if ( allocated ( T % Reference_ASC ) ) &
!!$      deallocate ( T % Reference_ASC )

!!$    call T % FinalizeTemplate_C_MS_C_PS ( )

  end subroutine Finalize


end module ThermalizationNeutrinos_Form
