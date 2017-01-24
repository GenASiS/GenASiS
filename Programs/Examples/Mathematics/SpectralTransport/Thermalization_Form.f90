module Thermalization_Form

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Matter_Form
  use Matter_ASC__Form
  use SetFermiDiracSpectrum_Command
  use Interactions_BSLL_ASC_CSLD__Form
  use RadiationMoments_BSLL_ASC_CSLD__Form

  implicit none
  private

  type, public, extends ( IntegratorTemplate ) :: ThermalizationForm
    type ( Matter_ASC_Form ), allocatable :: &
      Matter_ASC
    type ( Interactions_BSLL_ASC_CSLD_Form ), allocatable :: &
      Interactions_BSLL_ASC_CSLD
    type ( RadiationMoments_BSLL_ASC_CSLD_Form ), allocatable :: &
      RadiationMoments_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  end type ThermalizationForm

    private :: &
      SetMatter, &
      SetRadiation

    real ( KDR ), private, parameter :: &
      T_Scale_MeV = 3.0_KDR, &
      T_Min_MeV = 0.1_KDR, &
      T_Max_MeV = 10.0_KDR

contains


  subroutine Initialize ( T, Name )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T
    character ( * ), intent ( in )  :: &
      Name

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

    if ( T % Type == '' ) &
      T % Type = 'a Thermalization'

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: T % PositionSpace )
    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call PS % CreateChart ( )

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
    Scale ( 1 ) = T_Scale_MeV * UNIT % MEV

    CoordinateUnit = UNIT % IDENTITY
    CoordinateUnit ( 1 ) = UNIT % MEV

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

    !-- Matter

    allocate ( T % Matter_ASC )
    associate ( MA => T % Matter_ASC )
    call MA % Initialize &
           ( PS, TemperatureUnitOption = UNIT % MEV, &
             ChemicalPotentialUnitOption = UNIT % MEV )
    end associate !-- MA

    !-- Radiation

    EnergyDensityUnit  =  UNIT % MEV / UNIT % HBAR_C ** 3

    allocate ( T % RadiationMoments_BSLL_ASC_CSLD )
    associate ( RMB => T % RadiationMoments_BSLL_ASC_CSLD )
    call RMB % Initialize &
           ( MS, 'GENERIC', EnergyDensityUnitOption = EnergyDensityUnit )

    !-- Interactions

    allocate ( T % Interactions_BSLL_ASC_CSLD )
    associate ( IB => T % Interactions_BSLL_ASC_CSLD )
    call IB % Initialize ( MS, 'FIXED' )
    call RMB % SetInteractions ( IB )
    end associate !-- IB

    !-- Initial conditions

    call SetMatter ( T )
    call SetRadiation ( T )

    !-- Initialize template

    call T % InitializeTemplate &
           ( Name, FinishTimeOption = 0.0_KDR )

    !-- Cleanup

    end associate !-- RMB
    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  subroutine Finalize ( T )

    type ( ThermalizationForm ), intent ( inout ) :: &
      T

    if ( allocated ( T % RadiationMoments_BSLL_ASC_CSLD ) ) &
      deallocate ( T % RadiationMoments_BSLL_ASC_CSLD )
    if ( allocated ( T % Matter_ASC ) ) &
      deallocate ( T % Matter_ASC )

    call T % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetMatter ( T )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T

    real ( KDR ) :: &
      R_Min, R_Max, &
      T_Min, T_Max
    type ( CollectiveOperation_R_Form ) :: &
      CO
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( MatterForm ), pointer :: &
      M

    M => T % Matter_ASC % Matter ( )

    select type ( PS => T % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )

    associate &
      ( R0 => ( C % MaxCoordinate + C % MinCoordinate ) / 2.0_KDR, &
        L  => ( C % MaxCoordinate - C % MinCoordinate ), &
        X  => G % Value ( :, G % CENTER ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER ( 3 ) ), &
        Temperature => M % Value ( :, M % TEMPERATURE ), &
        ChemicalPotential => M % Value ( :, M % CHEMICAL_POTENTIAL ) )
    associate &
      ( R => sqrt ( ( X - R0 ( 1 ) ) ** 2  +  ( Y - R0 ( 2 ) ) ** 2 &
                    +  ( Z - R0 ( 3 ) ) ** 2 ) )

    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = minval ( R )
    CO % Outgoing % Value ( 2 ) = 1.0_KDR / maxval ( R )
    call CO % Reduce ( REDUCTION % MIN )

    R_Min = CO % Incoming % Value ( 1 )
    R_Max = 1.0_KDR / CO % Incoming % Value ( 2 )

    T_Min = T_Min_MeV * UNIT % MEV
    T_Max = T_Max_MeV * UNIT % MEV

!    Temperature &
!      = T_Max  -  ( T_Max - T_Min ) / ( R_Max - R_Min ) * ( R - R_Min )
    Temperature = T_Max * ( ( R_Min / R ) ** ( log10 ( T_Max / T_Min ) &
                                               / log10 ( R_Max / R_Min ) ) )
    ChemicalPotential = 0.0_KDR

    end associate !-- R
    end associate !-- R0, etc.
    end select !-- C
    end select !-- PS

    nullify ( G, M )

  end subroutine SetMatter


  subroutine SetRadiation ( T )

    class ( ThermalizationForm ), intent ( inout ) :: &
      T

    integer ( KDI ) :: &
      iF, &  !-- iFiber
      iE     !-- iEnergy  
    real ( KDR ) :: &
      Amplitude, &
      Perturbation
    class ( RadiationMomentsForm ), pointer :: &
      RM
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( RadiationMomentsForm ), allocatable :: &
      RS  !-- RadiationSection
    class ( MatterForm ), pointer :: &
      M

    associate ( RMB => T % RadiationMoments_BSLL_ASC_CSLD )

    select type ( MS => T % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )

    !-- Initialize primitive spectra in momentum space

    M => T % Matter_ASC % Matter ( )

    call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )

    do iF = 1, MS % nFibers
      associate ( iBC => MS % iaBaseCell ( iF ) )
      RM => RMB % RadiationMomentsFiber ( iF )
      associate &
        ( J   => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
          H_1 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
          H_2 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
          H_3 => RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ), &
          T   => M % Value ( iBC, M % TEMPERATURE ), &
          Mu  => M % Value ( iBC, M % CHEMICAL_POTENTIAL ), &
          E   => RMB % Energy )

      call SetFermiDiracSpectrum ( E, T, Mu, J )

      Amplitude = 0.9_KDR
      do iE = 1, RMB % nEnergyValues

        call random_number ( Perturbation )
        Perturbation = Amplitude * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        J ( iE ) = ( 1.0_KDR + Perturbation )  *  J ( iE )

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_1 ( iE ) = Perturbation

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_2 ( iE ) = Perturbation

        call random_number ( Perturbation )
        Perturbation &
          = 0.01_KDR * J ( iE ) * 2.0_KDR * ( Perturbation - 0.5_KDR ) 
        H_3 ( iE ) = Perturbation

      end do !-- iE

      end associate !-- J, etc.
      end associate !-- iBC
    end do !-- iF

    !-- Switch to position space to compute from primitive

    G => MS % Base_CSLD % Geometry ( )

    allocate ( RS )
    call RS % Initialize &
           ( RMB % Velocity_U_Unit, RMB % MomentumDensity_U_Unit, &
             RMB % MomentumDensity_D_Unit, RMB % EnergyDensityUnit, &
             G % nValues, ClearOption = .true. )

    do iE = 1, RMB % nEnergyValues
      call MS % LoadSection ( RS, RMB, iE )
      call RS % ComputeFromPrimitive ( G )
      call MS % StoreSection ( RMB, RS, iE )
    end do !-- iE

    !-- Cleanup

    end select !-- MS
    end associate !-- RMB
    nullify ( G, M, RM )

  end subroutine SetRadiation


end module Thermalization_Form
