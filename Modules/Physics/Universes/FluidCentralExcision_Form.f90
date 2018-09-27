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
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
    procedure, private, pass :: &
      InitializePositionSpace
    procedure, private, pass :: &
      InitializeGeometry
    procedure, private, pass :: &
      SetCoarsening
    procedure, public, nopass :: &
      CoarsenSingularities
  end type FluidCentralExcisionForm

    class ( FluidCentralExcisionForm ), private, pointer :: &
      FluidCentralExcision => null ( )

contains


  subroutine Initialize &
               ( FCE, Name, FluidType, GeometryType, DimensionlessOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 LimiterParameterOption, ShockThresholdOption, &
                 RadiusMaxOption, RadiusMinOption, RadialRatioOption, &
                 CentralMassOption, nWriteOption, nCellsPolarOption )

    class ( FluidCentralExcisionForm ), intent ( inout ), target :: &
      FCE
    character ( * ), intent ( in )  :: &
      Name, &
      FluidType, &
      GeometryType
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
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

    FluidCentralExcision => FCE

    call FCE % InitializeTemplate_FC &
               ( Name, FluidType, GeometryType, &
                 DimensionlessOption = DimensionlessOption, &
                 TimeUnitOption = TimeUnitOption, &
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

    select type ( S => FCE % Step )
    class is ( Step_RK2_C_ASC_Form )
      S % CoarsenSingularities => CoarsenSingularities
    end select !-- S

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCE )

    type ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE

    call FCE % FinalizeTemplate_FC ( )

  end subroutine Finalize


  subroutine InitializePositionSpace &
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

    allocate ( Atlas_SC_CE_Form :: FC % PositionSpace )
    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FC % Dimensionless ) then

      call PS % CreateChart_CE &
             ( RadiusMaxOption = RadiusMaxOption, &
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
             ( CoordinateUnitOption = FC % CoordinateUnit, &
               RadiusMaxOption = RadiusMax, &
               RadiusExcisionOption = RadiusExcision, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end if !-- Dimensionless

    end select !-- PS

  end subroutine InitializePositionSpace


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

    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_CE_Form )

    !-- Azimuthal coarsening

    if ( PS % nDimensions < 3 ) &
      return
    call Show ( 'SetCoarsening_3', FC % IGNORABILITY )
    call Show ( FC % Name, 'Universe', FC % IGNORABILITY )
    call FC % SetCoarseningTemplate ( iAngular = 3 )

    end select !-- PS

  end subroutine SetCoarsening


  subroutine CoarsenSingularities ( S )

    class ( StorageForm ), intent ( inout ) :: &
      S

    associate ( FCE => FluidCentralExcision )
if ( FCE % Time > 1.8e-3_KDR  *  UNIT % SECOND ) then
  call Show ( '>>> Before coarsening' )
  call FluidCentralExcision % Write ( )
end if

    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_CE_Form )

    if ( PS % nDimensions > 2 ) &
      call FCE % CoarsenSingularityTemplate ( S, iAngular = 3 )
  
if ( FCE % Time > 1.8e-3_KDR  *  UNIT % SECOND ) then
  call Show ( '>>> After coarsening' )
  call FluidCentralExcision % Write ( )
end if

    end select !-- PS
    end associate !-- FluidCentralExcision

  end subroutine CoarsenSingularities


end module FluidCentralExcision_Form
