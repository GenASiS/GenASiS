module GreyRadiationBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_Template ) :: GreyRadiationBoxForm
    integer ( KDI ) :: &
      FLUID = 0
    class ( Interactions_ASC_Template ), allocatable :: &
      Interactions_ASC
    class ( Relaxation_RM_G_Form ), allocatable :: &
      Relaxation_RM_G
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type GreyRadiationBoxForm

contains


  subroutine Initialize &
               ( GRB, RadiationName, RadiationType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 MinCoordinateOption, MaxCoordinateOption, TimeUnitOption, &
                 FinishTimeOption, CourantFactorOption, nCellsOption, &
                 nWriteOption )

    class ( GreyRadiationBoxForm ), intent ( inout ), target :: &
      GRB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in )  :: &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption
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
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions


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

    GRB % N_CURRENTS_PS  =  size ( RadiationName ) + 1  !-- Radiation + Fluid
    allocate ( GRB % Current_ASC_1D ( GRB % N_CURRENTS_PS ) )
    allocate ( GRB % TimeStepLabel ( GRB % N_CURRENTS_PS ) )

    do iC = 1, size ( RadiationName )
      GRB % TimeStepLabel ( iC )  =  RadiationName ( iC )
    end do !-- iC

    GRB % FLUID  =  GRB % N_CURRENTS_PS
    GRB % TimeStepLabel ( GRB % FLUID )  =  'Fluid'


    !-- Radiation

    do iC = 1, size ( RadiationName )
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )
        allocate &
          ( RadiationMoments_ASC_Form :: &
              GRB % Current_ASC_1D ( iC ) % Element )
        select type ( RA => GRB % Current_ASC_1D ( iC ) % Element )
        class is ( RadiationMoments_ASC_Form )
        call RA % Initialize &
               ( PS, RadiationType ( iC ), &
                 NameShortOption = RadiationName ( iC ) )
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


    !-- Fluid

    allocate ( Fluid_ASC_Form :: &
                 GRB % Current_ASC_1D ( GRB % FLUID ) % Element )
    select type ( FA => GRB % Current_ASC_1D ( GRB % FLUID ) % Element )
    class is ( Fluid_ASC_Form )
    call FA % Initialize &
           ( PS, 'IDEAL' )
!             TemperatureUnitOption = UNIT % MEGA_ELECTRON_VOLT )
    end select !-- FA


    !-- Relaxation

    allocate ( GRB % Relaxation_RM_G )
    associate ( R => GRB % Relaxation_RM_G )
    select type ( RMA => GRB % Current_ASC_1D ( 1 ) % Element )
    class is ( RadiationMoments_ASC_Form )
      call R % Initialize ( RMA, Name = GRB % Name )
    end select !-- RMA
    end associate !-- R


    !-- Step

    ApplyStreaming    = .true.
    ApplyInteractions = .true.
    if ( present ( ApplyStreamingOption ) ) &
      ApplyStreaming = ApplyStreamingOption
    if ( present ( ApplyInteractionsOption ) ) &
      ApplyInteractions = ApplyInteractionsOption

    allocate ( Step_RK2_C_ASC_1D_Form :: GRB % Step )
    select type ( S => GRB % Step )
    class is ( Step_RK2_C_ASC_1D_Form )
    call S % Initialize ( GRB, GRB % Current_ASC_1D, Name )

    if ( .not. ApplyStreaming ) then
      do iC = 1, size ( RadiationName )
        S % ApplyDivergence_1D ( iC ) % Pointer  =>  null ( )  
      end do !-- iC
    end if

    if ( ApplyInteractions ) then
      do iC = 1, size ( RadiationName )
        S % ApplyRelaxation_1D ( iC ) % Pointer  &
          =>  GRB % Relaxation_RM_G % Apply 
      end do !-- iC
    end if

    S % ApplyDivergence_1D ( GRB % FLUID ) % Pointer  =>  null ( )  
      !-- Disable fluid evolution

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

    if ( allocated ( GRB % Relaxation_RM_G ) ) &
      deallocate ( GRB % Relaxation_RM_G )
    if ( allocated ( GRB % Interactions_ASC ) ) &
      deallocate ( GRB % Interactions_ASC )

    call GRB % FinalizeTemplate_C_1D_PS ( )

  end subroutine Finalize


end module GreyRadiationBox_Form
