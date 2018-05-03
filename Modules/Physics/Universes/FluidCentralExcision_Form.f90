module FluidCentralExcision_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use ApplyGravity_F__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: &
    FluidCentralExcisionForm
      logical ( KDL ) :: &
        Dimensionless = .false.
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FluidCentralExcisionForm

contains


  subroutine Initialize &
               ( FCE, Name, FluidType, GeometryType, DimensionlessOption, &
                 TimeUnitOption, RadiusMaxOption, RadiusMinOption, &
                 RadialRatioOption, FinishTimeOption, nWriteOption )

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
      RadiusMaxOption, &
      RadiusMinOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    real ( KDR ) :: &
      RadiusMin, &
      RadiusMax, &
      RadialRatio!, &
!      FinishTime
!    type ( MeasuredValueForm ) :: &
!      TimeUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit


    if ( FCE % Type == '' ) &
      FCE % Type = 'a FluidCentralExcision'

    FCE % Dimensionless = .false.
    if ( present ( DimensionlessOption ) ) &
      FCE % Dimensionless = DimensionlessOption


    !-- PositionSpace

    allocate ( Atlas_SC_CE_Form :: FCE % PositionSpace )
    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )


    if ( FCE % Dimensionless ) then

      call PS % CreateChart_CE &
             ( RadiusMaxOption = RadiusMaxOption, &
               RadiusMinOption = RadiusMinOption, &
               RadialRatioOption = RadialRatioOption )

    else

      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadiusMin    =   40.0_KDR  *  UNIT % KILOMETER
      RadialRatio  =  6.5_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax = RadiusMaxOption
      if ( present ( RadiusMinOption ) ) &
        RadiusMin = RadiusMinOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio = RadialRatioOption

      CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

      call PS % CreateChart_CE &
             ( CoordinateUnitOption = CoordinateUnit, &
               RadiusMaxOption = RadiusMax, &
               RadiusMinOption = RadiusMin, &
               RadialRatioOption = RadialRatio )

    end if !-- Dimensionless


    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    if ( FCE % Dimensionless ) then
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS', &
               GravitationalConstantOption = 1.0_KDR )
    else
      call GA % Initialize &
             ( PS, GeometryType, GravitySolverTypeOption = 'CENTRAL_MASS' )
    end if !-- Dimensionless

    call PS % SetGeometry ( GA )


    !-- Cleanup

!    end select !-- FA
    end select !-- GA
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( FCE )

    type ( FluidCentralExcisionForm ), intent ( inout ) :: &
      FCE

    call FCE % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


end module FluidCentralExcision_Form
