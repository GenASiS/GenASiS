module GreyRadiationBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ) :: GreyRadiationBoxForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type GreyRadiationBoxForm

contains


  subroutine Initialize &
               ( GRB, RadiationType, Name, MinCoordinateOption, &
                 MaxCoordinateOption, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nCellsOption, nWriteOption )

    class ( GreyRadiationBoxForm ), intent ( inout ) :: &
      GRB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationType
    character ( * ), intent ( in )  :: &
      Name
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iC  !-- iCurrent


    if ( GRB % Type == '' ) &
      GRB % Type = 'a GreyRadiationBox'


    !-- PositionSpace

    allocate ( Atlas_SC_Form :: GRB % PositionSpace )
    select type ( PS => GRB % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    call PS % CreateChart &
           ( MinCoordinateOption = MinCoordinateOption, &
             MaxCoordinateOption = MaxCoordinateOption, &
             nCellsOption = nCellsOption )

    allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    call GA % Initialize ( PS, GeometryType = 'GALILEAN' )
    call PS % SetGeometry ( GA )
    end select !-- GA


    !-- Prepare for Currents

    GRB % N_CURRENTS_PS = size ( RadiationType )
    allocate ( GRB % Current_ASC_1D ( GRB % N_CURRENTS_PS ) )
    allocate ( GRB % TimeStepLabel ( GRB % N_CURRENTS_PS  ) )
    do iC = 1, GRB % N_CURRENTS_PS
      GRB % TimeStepLabel ( iC )  =  RadiationType ( iC )
    end do !-- iC
    

    !-- Radiation

    do iC = 1, GRB % N_CURRENTS_PS
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )
        allocate &
          ( RadiationMoments_ASC_Form :: &
              GRB % Current_ASC_1D ( iC ) % Element )
        select type ( RA => GRB % Current_ASC_1D ( iC ) % Element )
        class is ( RadiationMoments_ASC_Form )
        call RA % Initialize &
               ( PS, RadiationType ( iC ), &
                 NameShortOption = RadiationType ( iC ) )
                 ! Velocity_U_UnitOption = WHH % VelocityUnit, &
                 ! MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
                 ! MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
                 ! EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
                 ! TemperatureUnitOption = WHH % TemperatureUnit )
        end select !-- RA
      case default
        call Show ( 'RadiationType not implemented', CONSOLE % ERROR )
        call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
        call Show ( 'GreyRadiationBox_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end do !-- iC


    !-- Step

    allocate ( Step_RK2_C_ASC_1D_Form :: GRB % Step )
    select type ( S => GRB % Step )
    class is ( Step_RK2_C_ASC_1D_Form )
    call S % Initialize ( GRB, GRB % Current_ASC_1D, Name )
    end select !-- S


    !-- Template

    call GRB % InitializeTemplate_C_1D_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )


    !-- Cleanup

    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( GRB )

    type ( GreyRadiationBoxForm ), intent ( inout ) :: &
      GRB

    call GRB % FinalizeTemplate_C_1D_PS ( )

  end subroutine Finalize


end module GreyRadiationBox_Form
