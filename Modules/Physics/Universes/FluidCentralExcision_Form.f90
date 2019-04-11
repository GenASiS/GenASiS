module FluidCentralExcision_Form

  use Basics
  use Mathematics
  use Spaces
  use FluidCentral_Template

  implicit none
  private

  type, public, extends ( FluidCentralTemplate ) :: &
    FluidCentralExcisionForm
  contains
    procedure, private, pass :: &
      Initialize_FCE
    generic, public :: &
      Initialize => Initialize_FCE
    final :: &
      Finalize
    procedure, private, pass :: &
      InitializeAtlas
    procedure, private, pass :: &
      InitializeGeometry
    procedure, public, pass :: &
      SetCoarsening
  end type FluidCentralExcisionForm

    private :: &
      CoarsenSingularities

contains


  subroutine Initialize_FCE &
               ( FCE, FluidType, GeometryType, Name, DimensionlessOption, &
                 FinishTimeOption, CourantFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 RadiusMaxOption, RadiusMinOption, RadialRatioOption, &
                 CentralMassOption, nWriteOption, nCellsPolarOption )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE
    character ( * ), intent ( in )  :: &
      FluidType, &
      GeometryType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      LimiterParameterOption, &
      ShockThresholdOption, &
      RadiusMaxOption, &
      RadiusMinOption, &
      RadialRatioOption, &
      CentralMassOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption, &
      nCellsPolarOption

    if ( FCE % Type == '' ) &
      FCE % Type = 'a FluidCentralExcision'

    call FCE % InitializeTemplate_FC &
               ( FluidType, GeometryType, Name, &
                 DimensionlessOption = DimensionlessOption, &
                 FinishTimeOption = FinishTimeOption, &
                 CourantFactorOption = CourantFactorOption, &
                 LimiterParameterOption = LimiterParameterOption, &
                 ShockThresholdOption = ShockThresholdOption, &
                 RadiusMaxOption = RadiusMaxOption, &
                 RadiusMinOption = RadiusMinOption, &
                 RadialRatioOption = RadialRatioOption, &
                 CentralMassOption = CentralMassOption, &
                 nWriteOption = nWriteOption, &
                 nCellsPolarOption = nCellsPolarOption )

    select type ( I => FCE % Integrator )
    class is ( Integrator_C_PS_Form )
    select type ( S => I % Step )
    class is ( Step_RK_C_ASC_Template )
      S % CoarsenSingularities => CoarsenSingularities
    end select !-- S
    end select !-- I

  end subroutine Initialize_FCE


  impure elemental subroutine Finalize ( FCE )

    type ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE

    call FCE % FinalizeTemplate_FC ( )

  end subroutine Finalize


  subroutine InitializeAtlas &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusExcisionOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusExcision, &
      RadialRatio

    associate ( I => FC % Integrator )

    allocate ( Atlas_SC_CE_Form :: I % PositionSpace )
    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FC % Dimensionless ) then

      call PS % CreateChart_CE &
             ( UseCustomBoundaryInnerOption = FC % UseCustomBoundaryInner, &
               RadiusMaxOption = RadiusMaxOption, &
               RadiusExcisionOption = RadiusExcisionOption, &
               RadialRatioOption = RadialRatioOption )

    else

      RadiusMax       =  1.0e3_KDR  *  UNIT % KILOMETER
      RadiusExcision  =   40.0_KDR  *  UNIT % KILOMETER
      RadialRatio     =  1.0_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax = RadiusMaxOption
      if ( present ( RadiusExcisionOption ) ) &
        RadiusExcision = RadiusExcisionOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio = RadialRatioOption

      call PS % CreateChart_CE &
             ( UseCustomBoundaryInnerOption = FC % UseCustomBoundaryInner, &
               CoordinateUnitOption = FC % Units % Coordinate_PS, &
               RadiusMaxOption = RadiusMax, &
               RadiusExcisionOption = RadiusExcision, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end if !-- Dimensionless

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


  subroutine InitializeGeometry &
               ( FC, GA, PS, GeometryType, CentralMassOption )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FC
    type ( Geometry_ASC_Form ), intent ( inout ) :: &
      GA
    class ( Atlas_SC_Form ), intent ( in ) :: &
      PS
    character ( * ), intent ( in )  :: &
      GeometryType
    real ( KDR ), intent ( in ), optional :: &
      CentralMassOption

    real ( KDR ) :: &
      CentralMass

    if ( FC % Dimensionless ) then
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS', &
               GravitationalConstantOption = 1.0_KDR, &
               CentralMassOption = 1.0_KDR )
    else
      CentralMass  =  1.4_KDR  *  UNIT % SOLAR_MASS
      if ( present ( CentralMassOption ) ) &
        CentralMass = CentralMassOption
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS', &
               CentralMassUnitOption = UNIT % SOLAR_MASS, &
               CentralMassOption = CentralMass )
    end if !-- Dimensionless

  end subroutine InitializeGeometry


  subroutine SetCoarsening ( FC )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FC

    select type ( I => FC % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CE_Form )

    !-- Azimuthal coarsening

    if ( PS % nDimensions < 3 ) &
      return
    call Show ( 'SetCoarsening_3', FC % IGNORABILITY )
    call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
    call FC % SetCoarseningTemplate ( iAngular = 3 )

    end select !-- PS
    end select !-- I

  end subroutine SetCoarsening


  subroutine CoarsenSingularities ( S, Increment )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( StorageForm ), intent ( inout ) :: &
      Increment

    select type ( FCE => S % Integrator % Universe )
    class is ( FluidCentralExcisionForm )

    select type ( I => S % Integrator )
    class is ( IntegratorTemplate )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_CE_Form )

    if ( PS % nDimensions > 2 ) &
      call FCE % CoarsenSingularityTemplate ( Increment, iAngular = 3 )
  
    end select !-- PS
    end select !-- I
    end select !-- FCE

  end subroutine CoarsenSingularities


end module FluidCentralExcision_Form
