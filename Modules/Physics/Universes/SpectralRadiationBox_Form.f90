module SpectralRadiationBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_1D_MS_C_PS_Template ) :: &
    SpectralRadiationBoxForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SpectralRadiationBoxForm

contains


  subroutine Initialize &
               ( SRB, RadiationType, Name, MinCoordinateOption, &
                 MaxCoordinateOption, TimeUnitOption, FinishTimeOption, &
                 CourantFactorOption, nCellsPositionOption, nWriteOption )

    class ( SpectralRadiationBoxForm ), intent ( inout ) :: &
      SRB
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
      nCellsPositionOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iC, &  !-- iCurrents
      nCellsEnergy
    integer ( KDI ), dimension ( 3 ) :: &
      nCellsPosition


    if ( SRB % Type == '' ) &
      SRB % Type = 'a SpectralRadiationBox'


    !-- PositionSpace

    allocate ( Atlas_SC_Form :: SRB % PositionSpace )
    select type ( PS => SRB % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

    nCellsPosition = [ 32, 32, 32 ]
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

    allocate ( Bundle_SLL_ASC_CSLD_Form :: SRB % MomentumSpace )
    select type ( MS => SRB % MomentumSpace )
    class is ( Bundle_SLL_ASC_CSLD_Form )
    call MS % Initialize ( PS, 'MomentumSpace' )
    call MS % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

    nCellsEnergy = 4
    call PROGRAM_HEADER % GetParameter ( nCellsEnergy, 'nCellsEnergy' )

    call MS % CreateChart &
           ( CoordinateSystemOption = 'SPHERICAL', &
             nCellsOption = [ nCellsEnergy ], &
             nGhostLayersOption = [ 0, 0, 0 ] )


    !-- Prepare for Currents

    SRB % N_CURRENTS_MS = size ( RadiationType )
    allocate ( SRB % Current_BSLL_ASC_CSLD_1D ( SRB % N_CURRENTS_MS ) )
    allocate ( SRB % TimeStepLabel ( SRB % N_CURRENTS_MS  ) )
    do iC = 1, SRB % N_CURRENTS_MS
      SRB % TimeStepLabel ( iC )  =  RadiationType ( iC )
    end do !-- iC
    

    !-- Radiation

    do iC = 1, SRB % N_CURRENTS_MS
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )
        allocate &
          ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
              SRB % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        select type ( RMB => SRB % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
        class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
        call RMB % Initialize &
               ( MS, RadiationType ( iC ), &
                 NameShortOption = RadiationType ( iC ), &
                 UseLimiterOption = .true., &
                 LimiterParameterOption = 1.4_KDR )
                 ! Velocity_U_UnitOption = WHH % VelocityUnit, &
                 ! MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
                 ! MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
                 ! EnergyDensityUnitOption = WHH % EnergyDensityUnit, &
                 ! TemperatureUnitOption = WHH % TemperatureUnit )
        end select !-- RMB
      case default
        call Show ( 'RadiationType not implemented', CONSOLE % ERROR )
        call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
        call Show ( 'SpectralRadiationBox_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end do !-- iC


    ! !-- Fluid ( This doesn't do anything, just a Step interface placeholder )

    ! allocate ( Fluid_ASC_Form :: PWS % Current_ASC )
    ! select type ( FA => PWS % Current_ASC )
    ! class is ( Fluid_ASC_Form )
    ! call FA % Initialize ( PS, 'NON_RELATIVISTIC' )


    !-- Step

    allocate ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form :: SRB % Step_MS )
    select type ( S_MS => SRB % Step_MS )
    class is ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form )
    call S_MS % Initialize ( SRB, SRB % Current_BSLL_ASC_CSLD_1D, Name )
    end select !-- S_MS

    ! allocate ( Step_RK2_C_ASC_Form :: PWS % Step_PS )
    ! select type ( S_PS => PWS % Step_PS )
    ! class is ( Step_RK2_C_ASC_Form )
    ! call S_PS % Initialize ( FA, Name )
    ! S_PS % ApplyDivergence % Pointer => null ( )  !-- Disable fluid evolution
    ! end select !-- S


    !-- Template

    call SRB % InitializeTemplate_C_1D_MS_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )


    !-- Cleanup

    end select !-- MS
    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( SRB )

    type ( SpectralRadiationBoxForm ), intent ( inout ) :: &
      SRB

    call SRB % FinalizeTemplate_C_1D_MS_C_PS ( )

  end subroutine Finalize

  
end module SpectralRadiationBox_Form
