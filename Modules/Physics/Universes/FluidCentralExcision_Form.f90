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
    procedure, public, pass :: &
      CoarsenSingularities
  end type FluidCentralExcisionForm

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

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCE )

    type ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE

    call FCE % FinalizeTemplate_FC ( )

  end subroutine Finalize


  subroutine InitializePositionSpace &
               ( FC, RadiusMaxOption, RadiusCoreOption, RadiusMinOption, &
                 RadialRatioOption, nCellsPolarOption )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FC
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadiusMinOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusMin, &
      RadialRatio

    allocate ( Atlas_SC_CE_Form :: FC % PositionSpace )
    select type ( PS => FC % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    if ( FC % Dimensionless ) then

      call PS % CreateChart_CE &
             ( RadiusMaxOption = RadiusMaxOption, &
               RadiusMinOption = RadiusMinOption, &
               RadialRatioOption = RadialRatioOption )

    else

      RadiusMax    =  1.0e3_KDR  *  UNIT % KILOMETER
      RadiusMin    =   40.0_KDR  *  UNIT % KILOMETER
      RadialRatio  =  1.0_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax = RadiusMaxOption
      if ( present ( RadiusMinOption ) ) &
        RadiusMin = RadiusMinOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio = RadialRatioOption

      call PS % CreateChart_CE &
             ( CoordinateUnitOption = FC % CoordinateUnit, &
               RadiusMaxOption = RadiusMax, &
               RadiusMinOption = RadiusMin, &
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

    !-- FIXME: Fill in along the lines of FluidCentralCore_Form

  end subroutine SetCoarsening


  subroutine CoarsenSingularities ( FC, S, iAngular )

    class ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FC
    class ( StorageForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iAngular

    !-- FIXME: Fill in along the lines of FluidCentralCore_Form

  end subroutine CoarsenSingularities


end module FluidCentralExcision_Form
