module RadiationBox_G_OS__Form

  !-- RadiationBox_Grey_OperatorSplit_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  type, public, extends ( Integrator_C_1D_PS_C_PS_Template ) :: &
    RadiationBox_G_OS_Form
      class ( Interactions_ASC_Template ), allocatable :: &
        Interactions_ASC
      class ( Relaxation_RM_ASC_Form ), allocatable :: &
        Relaxation_RM_ASC
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RadiationBox_G_OS_Form

contains


  subroutine Initialize &
               ( RB, RadiationName, RadiationType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption, MinCoordinateOption, MaxCoordinateOption, &
                 TimeUnitOption, FinishTimeOption, CourantFactorOption, &
                 nCellsOption, nWriteOption )

    class ( RadiationBox_G_OS_Form ), intent ( inout ), target :: &
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
      CourantFactorOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsOption
    integer ( KDI ), intent ( in ), optional :: &
      nWriteOption

    integer ( KDI ) :: &
      iC  !-- iCurrent
    logical ( KDL ) :: &
      ApplyStreaming, &
      ApplyInteractions, &
      EvolveFluid


    if ( RB % Type == '' ) &
      RB % Type = 'a RadiationBox_G_OS'


    !-- PositionSpace

    allocate ( Atlas_SC_Form :: RB % PositionSpace )
    select type ( PS => RB % PositionSpace )
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

    RB % N_CURRENTS_1D  =  size ( RadiationName )
    allocate ( RB % Current_ASC_1D ( RB % N_CURRENTS_1D ) )
    allocate ( RB % TimeStepLabel ( RB % N_CURRENTS_1D  +  1 ) )
    do iC = 1, RB % N_CURRENTS_1D
      RB % TimeStepLabel ( iC )  =  RadiationName ( iC )
    end do !-- iC
    RB % TimeStepLabel ( RB % N_CURRENTS_1D  +  1 )  =  'Fluid'


    !-- Radiation

    do iC = 1, RB % N_CURRENTS_1D
      select case ( trim ( RadiationType ( iC ) ) )
      case ( 'GENERIC' )
        allocate &
          ( RadiationMoments_ASC_Form :: &
              RB % Current_ASC_1D ( iC ) % Element )
        select type ( RA => RB % Current_ASC_1D ( iC ) % Element )
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
        call Show ( 'RadiationBox_G_OS__Form', 'module', CONSOLE % ERROR )
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
      allocate ( RB % Relaxation_RM_ASC )
      associate ( R => RB % Relaxation_RM_ASC )
      select type ( RMA => RB % Current_ASC_1D ( 1 ) % Element )
      class is ( RadiationMoments_ASC_Form )
        call R % Initialize ( RMA, Name = RB % Name )
      end select !-- RMA
      end associate !-- R
    end if !-- ApplyInteractions


    !-- Step

    allocate ( Step_RK2_C_ASC_1D_Form :: RB % Step_1D ) !-- Radiation
    select type ( S_1D => RB % Step_1D )
    class is ( Step_RK2_C_ASC_1D_Form )
    call S_1D % Initialize ( RB, RB % Current_ASC_1D, Name )

    if ( .not. ApplyStreaming ) then
      do iC = 1, RB % N_CURRENTS_1D
        S_1D % ApplyDivergence_1D ( iC ) % Pointer  =>  null ( )  
      end do !-- iC
    end if

    if ( ApplyInteractions ) then
      do iC = 1, RB % N_CURRENTS_1D
        S_1D % ApplyRelaxation_1D ( iC ) % Pointer  &
          =>  RB % Relaxation_RM_ASC % Apply 
      end do !-- iC
    end if

    end select !-- S_1D

    allocate ( Step_RK2_C_ASC_Form :: RB % Step ) !-- Fluid
    select type ( S => RB % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( RB, RB % Current_ASC, Name )
    if ( .not. EvolveFluid ) &
      S % ApplyDivergence % Pointer  =>  null ( )  

    end select !-- S


    !-- Template

    call RB % InitializeTemplate_C_1D_PS_C_PS &
           ( Name, TimeUnitOption = TimeUnitOption, &
             FinishTimeOption = FinishTimeOption, &
             CourantFactorOption = CourantFactorOption, &
             nWriteOption = nWriteOption )


    !-- Cleanup

    end select !-- PS

  end subroutine Initialize


  impure elemental subroutine Finalize ( RB )

    type ( RadiationBox_G_OS_Form ), intent ( inout ) :: &
      RB

    if ( allocated ( RB % Relaxation_RM_ASC ) ) &
      deallocate ( RB % Relaxation_RM_ASC )
    if ( allocated ( RB % Interactions_ASC ) ) &
      deallocate ( RB % Interactions_ASC )

    call RB % FinalizeTemplate_C_1D_PS_C_PS ( )

  end subroutine Finalize


end module RadiationBox_G_OS__Form
