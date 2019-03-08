module RadiationBox_Form

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies
  use Universe_Template

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RadiationBoxForm
!     real ( KDR ) :: &
!       InteractionFactor
!     logical ( KDL ) :: &
!       ApplyStreaming, &
!       ApplyInteractions, &
!       EvolveFluid
!     character ( LDL ) :: &
!       MomentsType = ''
!     class ( Interactions_ASC_Template ), allocatable :: &
!       Interactions_ASC
!     class ( Interactions_BSLL_ASC_CSLD_Template ), allocatable :: &
!       Interactions_BSLL_ASC_CSLD
!     class ( Relaxation_RM_ASC_Form ), allocatable :: &
!       Relaxation_RM_ASC
!     class ( Relaxation_RM_BSLL_ASC_CSLD_Form ), allocatable :: &
!       Relaxation_RM_BSLL_ASC_CSLD
  contains
    procedure, public, pass :: &
      Initialize_RB
    generic, public :: &
      Initialize => Initialize_RB
!     final :: &
!       Finalize
  end type RadiationBoxForm

!     private :: &
!       AllocateIntegrator, &
!       InitializePositionSpace, &
!       InitializeMomentumSpace, &
!       InitializeRadiation, &
!       InitializeFluid, &
!       InitializeSteps, &
!       ComputeTimeStepLocal, &
!       PrepareStep, &
!       IntegrateSources, &
!       ApplySources_Fluid

!     private :: &
!         ComputeTimeStepInteractions, &
!         ComputeFluidSource_G_S_Radiation_Kernel, &
!         ApplySources_Fluid_Kernel

!     class ( RadiationBoxForm ), pointer :: &
!       RadiationBox => null ( )

contains

  
  subroutine Initialize_RB &
               ( RB, RadiationName, RadiationType, MomentsType, Name, &
                 ApplyStreamingOption, ApplyInteractionsOption, &
                 EvolveFluidOption, MinCoordinateOption, MaxCoordinateOption, &
                 FinishTimeOption, CourantFactorOption, EnergyScaleOption, &
                 BaryonMassReferenceOption, nCellsPositionOption, &
                 nCellsEnergyOption, nWriteOption )

    class ( RadiationBoxForm ), intent ( inout ), target :: &
      RB
    character ( * ), dimension ( : ), intent ( in )  :: &
      RadiationName, &
      RadiationType
    character ( * ), intent ( in ) :: &
      MomentsType, &
      Name
    logical ( KDL ), intent ( in ), optional :: &
      ApplyStreamingOption, &
      ApplyInteractionsOption, &
      EvolveFluidOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      CourantFactorOption, &
      EnergyScaleOption, &
      BaryonMassReferenceOption
    integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
      nCellsPositionOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsEnergyOption, &
      nWriteOption

!     type ( MeasuredValueForm ) :: &
!       CoordinateUnit_PS, &
!       VelocityUnit, &
!       BaryonMassUnit, &
!       NumberDensityUnit, &
!       EnergyDensityUnit, &
!       NumberUnit, &
!       EnergyUnit, &
!       MomentumUnit, &
!       AngularMomentumUnit

!     if ( RB % Type == '' ) &
!       RB % Type = 'a RadiationBox'
    
!     call RB % InitializeTemplate ( Name )

!     RadiationBox => RB

!     RB % MomentsType  =  MomentsType

!     call AllocateIntegrator &
!            ( RB, RadiationName )
!     call InitializePositionSpace &
!            ( RB, BoundaryConditionsFaceOption, CoordinateUnit_PS_Option, &
!              MinCoordinateOption, MaxCoordinateOption, nCellsPositionOption )
!     call InitializeMomentumSpace &
!            ( RB, CoordinateUnit_MS_Option, EnergyScaleOption, &
!              nCellsEnergyOption )

!     if ( present ( CoordinateUnit_PS_Option ) ) &
!       CoordinateUnit_PS = CoordinateUnit_PS_Option ( 1 )
!     if ( present ( Velocity_U_UnitOption ) ) &
!       VelocityUnit = Velocity_U_UnitOption ( 1 )
!     if ( present ( BaryonMassUnitOption ) ) &
!       BaryonMassUnit = BaryonMassUnitOption
!     if ( present ( NumberDensityUnitOption ) ) &
!       NumberDensityUnit = NumberDensityUnitOption
!     if ( present ( EnergyDensityUnitOption ) ) &
!       EnergyDensityUnit = EnergyDensityUnitOption

!     associate ( nD => RB % Integrator % PositionSpace % nDimensions )
!     NumberUnit    =  NumberDensityUnit  *  CoordinateUnit_PS ** nD
!     EnergyUnit    =  EnergyDensityUnit  *  CoordinateUnit_PS ** nD
!     MomentumUnit  =  BaryonMassUnit  *  CoordinateUnit_PS ** (-3) &
!                      *  VelocityUnit  *  CoordinateUnit_PS ** nD
!     end associate !-- nD
!     AngularMomentumUnit  =  MomentumUnit  *  CoordinateUnit_PS

!     call InitializeRadiation &
!            ( RB, RadiationName, RadiationType, Velocity_U_UnitOption, &
!              MomentumDensity_U_UnitOption, MomentumDensity_D_UnitOption, &
!              BaryonMassUnitOption, NumberDensityUnitOption, &
!              EnergyDensityUnitOption, TemperatureUnitOption, &
!              NumberUnit, EnergyUnit, MomentumUnit, AngularMomentumUnit, &
!              TimeUnitOption )
!     call InitializeFluid &
!            ( RB, Velocity_U_UnitOption, MomentumDensity_D_UnitOption, &
!              BaryonMassUnitOption, NumberDensityUnitOption, &
!              EnergyDensityUnitOption, TemperatureUnitOption, NumberUnit, &
!              EnergyUnit, MomentumUnit, AngularMomentumUnit, &
!              TimeUnitOption, BaryonMassReferenceOption )
!     call InitializeSteps &
!            ( RB, Name, ApplyStreamingOption, ApplyInteractionsOption, &
!              EvolveFluidOption )

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_C_PS_Template )
!       call I % Initialize &
!              ( Name, TimeUnitOption = TimeUnitOption, &
!                FinishTimeOption = FinishTimeOption, &
!                CourantFactorOption = CourantFactorOption, &
!                nWriteOption = nWriteOption )
!       I % ComputeTimeStepLocal => ComputeTimeStepLocal
!     end select !-- I

!     RB % InteractionFactor  =  1.0e-2_KDR
!     call PROGRAM_HEADER % GetParameter &
!            ( RB % InteractionFactor, 'InteractionFactor' )
!     call Show ( RB % InteractionFactor, 'InteractionFactor', &
!                 RB % IGNORABILITY )

  end subroutine Initialize_RB


!   subroutine Finalize ( RB )

!     type ( RadiationBoxForm ), intent ( inout ) :: &
!       RB

!     if ( allocated ( RB % Relaxation_RM_BSLL_ASC_CSLD ) ) &
!       deallocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
!     if ( allocated ( RB % Relaxation_RM_ASC ) ) &
!       deallocate ( RB % Relaxation_RM_ASC )

!     if ( allocated ( RB % Interactions_BSLL_ASC_CSLD ) ) &
!       deallocate ( RB % Interactions_BSLL_ASC_CSLD )
!     if ( allocated ( RB % Interactions_ASC ) ) &
!       deallocate ( RB % Interactions_ASC )

!     call RB % FinalizeTemplate ( )

!   end subroutine Finalize


!   subroutine AllocateIntegrator ( RB, RadiationName )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     character ( * ), dimension ( : ), intent ( in )  :: &
!       RadiationName

!     integer ( KDI ) :: &
!       iC  !-- iCurrent

!     select case ( trim ( RB % MomentsType ) )
!     case ( 'GREY' )
!       allocate ( Integrator_C_1D_PS_C_PS_Form :: RB % Integrator )
!     case ( 'SPECTRAL' )
!       allocate ( Integrator_C_1D_MS_C_PS_Form :: RB % Integrator )
!     case default
!       call Show ( 'MomentsType not recognized', CONSOLE % ERROR )
!       call Show ( RB % MomentsType, 'MomentsType', CONSOLE % ERROR )
!       call Show ( 'RadiationBox_Form', 'module', CONSOLE % ERROR )
!       call Show ( 'AllocateIntegrator', 'subroutine', CONSOLE % ERROR )
!       call PROGRAM_HEADER % Abort ( )
!     end select !-- MomentsType

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_C_PS_Template )

!       I % N_CURRENTS_1D  =  size ( RadiationName )
!       allocate ( I % TimeStepLabel &
!                    ( I % N_CURRENTS_1D  +  1  +  I % N_CURRENTS_1D ) )

!       do iC = 1, I % N_CURRENTS_1D
!         I % TimeStepLabel ( iC )  &
!           =  trim ( RadiationName ( iC ) ) // ' Streaming'
!       end do !-- iC

!       I % TimeStepLabel ( I % N_CURRENTS_1D  +  1 )  &
!           =  'Fluid Advection'

!       do iC = 1, I % N_CURRENTS_1D
!         I % TimeStepLabel ( I % N_CURRENTS_1D  +  1  +  iC )  &
!           =  trim ( RadiationName ( iC ) ) // ' Interactions'
!       end do !-- iC

!     end select !-- I

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_PS_C_PS_Form )
!       allocate ( I % Current_ASC_1D ( I % N_CURRENTS_1D ) )
!     class is ( Integrator_C_1D_MS_C_PS_Form )
!       allocate ( I % Current_BSLL_ASC_CSLD_1D ( I % N_CURRENTS_1D ) )
!     end select !-- I

!   end subroutine AllocateIntegrator


!   subroutine InitializePositionSpace &
!                ( RB, BoundaryConditionsFaceOption, CoordinateUnit_PS_Option, &
!                  MinCoordinateOption, MaxCoordinateOption, &
!                  nCellsPositionOption )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     type ( Character_1D_Form ), dimension ( : ), intent ( in ), optional :: &
!       BoundaryConditionsFaceOption
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       CoordinateUnit_PS_Option
!     real ( KDR ), dimension ( : ), intent ( in ), optional :: &
!       MinCoordinateOption, &
!       MaxCoordinateOption
!     integer ( KDI ), dimension ( 3 ), intent ( in ), optional :: &
!       nCellsPositionOption

!     integer ( KDI ) :: &
!       iD  !-- iDimension
!     integer ( KDI ), dimension ( 3 ) :: &
!       nCellsPosition

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_C_PS_Template )

!     allocate ( Atlas_SC_Form :: I % PositionSpace )
!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )
!     call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )

!     if ( present ( BoundaryConditionsFaceOption ) ) then
!       do iD = 1, PS % nDimensions
!         call PS % SetBoundaryConditionsFace &
!                ( BoundaryConditionsFaceOption ( iD ) % Value, &
!                  iDimension = iD )
!       end do !-- iD
!     end if

!     nCellsPosition = [ 32, 32, 32 ]
!     if ( present ( nCellsPositionOption ) ) &
!       nCellsPosition = nCellsPositionOption
!     call PROGRAM_HEADER % GetParameter ( nCellsPosition, 'nCellsPosition' )

!     call PS % CreateChart &
!            ( CoordinateUnitOption = CoordinateUnit_PS_Option, &
!              MinCoordinateOption = MinCoordinateOption, &
!              MaxCoordinateOption = MaxCoordinateOption, &
!              nCellsOption = nCellsPosition )

!     allocate ( Geometry_ASC_Form :: PS % Geometry_ASC )
!     select type ( GA => PS % Geometry_ASC )
!     class is ( Geometry_ASC_Form )
!     call GA % Initialize ( PS, GeometryType = 'GALILEAN' )
!     call PS % SetGeometry ( GA )

!     end select !-- GA
!     end select !-- PS
!     end select !-- I

!   end subroutine InitializePositionSpace


!   subroutine InitializeMomentumSpace &
!                ( RB, CoordinateUnit_MS_Option, EnergyScaleOption, &
!                  nCellsEnergyOption )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       CoordinateUnit_MS_Option
!     real ( KDR ), intent ( in ), optional :: &
!       EnergyScaleOption
!     integer ( KDI ), intent ( in ), optional :: &
!       nCellsEnergyOption

!     integer ( KDI ) :: &
!       nCellsEnergy
!     real ( KDR ) :: &
!       EnergyScale

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_MS_C_PS_Form )

!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )

!     allocate ( Bundle_SLL_ASC_CSLD_Form :: I % MomentumSpace )
!     select type ( MS => I % MomentumSpace )
!     class is ( Bundle_SLL_ASC_CSLD_Form )
!     call MS % Initialize ( PS, 'MomentumSpace' )
!     call MS % SetBoundaryConditionsFace &
!            ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 1 )

!     EnergyScale = 10.0_KDR
!     if ( present ( EnergyScaleOption ) ) &
!       EnergyScale = EnergyScaleOption
!     call PROGRAM_HEADER % GetParameter ( EnergyScale, 'EnergyScale' )

! !    CoordinateUnit = UNIT % IDENTITY
! !    CoordinateUnit ( 1 ) = UNIT % MEGA_ELECTRON_VOLT

!     nCellsEnergy = 16
!     if ( present ( nCellsEnergyOption ) ) &
!       nCellsEnergy = nCellsEnergyOption
!     call PROGRAM_HEADER % GetParameter ( nCellsEnergy, 'nCellsEnergy' )

!     call MS % CreateChart &
!            ( SpacingOption = [ 'COMPACTIFIED' ], &
!              CoordinateSystemOption = 'SPHERICAL', &
!              CoordinateUnitOption = CoordinateUnit_MS_Option, &
!              ScaleOption = [ EnergyScale ], &
!              nCellsOption = [ nCellsEnergy ], &
!              nGhostLayersOption = [ 0, 0, 0 ] )

!     end select !-- MS
!     end select !-- PS
!     end select !-- I

!   end subroutine InitializeMomentumSpace


!   subroutine InitializeRadiation &
!                ( RB, RadiationName, RadiationType, Velocity_U_UnitOption, &
!                  MomentumDensity_U_UnitOption, MomentumDensity_D_UnitOption, &
!                  BaryonMassUnitOption, NumberDensityUnitOption, &
!                  EnergyDensityUnitOption, TemperatureUnitOption, &
!                  NumberUnitOption, EnergyUnitOption, MomentumUnitOption, &
!                  AngularMomentumUnitOption, TimeUnitOption )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     character ( * ), dimension ( : ), intent ( in )  :: &
!       RadiationName, &
!       RadiationType
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       Velocity_U_UnitOption, &
!       MomentumDensity_U_UnitOption, &
!       MomentumDensity_D_UnitOption
!     type ( MeasuredValueForm ), intent ( in ), optional :: &
!       BaryonMassUnitOption, &
!       NumberDensityUnitOption, &
!       EnergyDensityUnitOption, &
!       TemperatureUnitOption, &
!       NumberUnitOption, &
!       EnergyUnitOption, &
!       MomentumUnitOption, &
!       AngularMomentumUnitOption, &
!       TimeUnitOption

!     integer ( KDI ) :: &
!       iC  !-- iCurrent
!     character ( LDL ) :: &
!       RadiationTypeLocal

!     select type ( I_1D => RB % Integrator )
!     class is ( Integrator_C_1D_C_PS_Template )

!     do iC = 1, I_1D % N_CURRENTS_1D

!       select type ( I => I_1D )
!       class is ( Integrator_C_1D_PS_C_PS_Form )

!         select type ( PS => I % PositionSpace )
!         class is ( Atlas_SC_Form )

!         select case ( trim ( RadiationType ( iC ) ) )
!         case ( 'GENERIC' )
!           allocate &
!             ( RadiationMoments_ASC_Form :: I % Current_ASC_1D ( iC ) % Element )
!           RadiationTypeLocal = 'GENERIC'
!         case ( 'PHOTONS' )
!           allocate &
!             ( PhotonMoments_ASC_Form :: I % Current_ASC_1D ( iC ) % Element )
!           RadiationTypeLocal = 'PHOTONS_GREY'
!         case default
!           call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
!           call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
!           call Show ( 'RadiationBox_Form', 'module', CONSOLE % ERROR )
!           call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
!           call PROGRAM_HEADER % Abort ( )
!         end select !-- RadiationType
       
!         select type ( RA => I % Current_ASC_1D ( iC ) % Element )
!         class is ( RadiationMoments_ASC_Form )
!           call RA % Initialize &
!                  ( PS, RadiationTypeLocal, &
!                    NameShortOption = RadiationName ( iC ), &
!                    Velocity_U_UnitOption = Velocity_U_UnitOption, &
!                    MomentumDensity_U_UnitOption &
!                      = MomentumDensity_U_UnitOption, &
!                    MomentumDensity_D_UnitOption &
!                      = MomentumDensity_D_UnitOption, &
!                    EnergyDensityUnitOption = EnergyDensityUnitOption, &
!                    TemperatureUnitOption = TemperatureUnitOption, &
!                    EnergyUnitOption = EnergyUnitOption, &
!                    MomentumUnitOption = MomentumUnitOption, &
!                    AngularMomentumUnitOption = AngularMomentumUnitOption, &
!                    TimeUnitOption = TimeUnitOption )
!         end select !-- RA        
!         end select !-- PS

!       class is ( Integrator_C_1D_MS_C_PS_Form )

!         select type ( MS => I % MomentumSpace )
!         class is ( Bundle_SLL_ASC_CSLD_Form )

!         select case ( trim ( RadiationType ( iC ) ) )
!         case ( 'GENERIC' )
!           allocate &
!             ( RadiationMoments_BSLL_ASC_CSLD_Form :: &
!                 I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
!           RadiationTypeLocal = 'GENERIC'
!         case ( 'PHOTONS' )
!           allocate &
!             ( PhotonMoments_BSLL_ASC_CSLD_Form :: &
!                 I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
!           RadiationTypeLocal = 'PHOTONS_SPECTRAL'
!         case default
!           call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
!           call Show ( RadiationType ( iC ), 'RadiationType', CONSOLE % ERROR )
!           call Show ( 'RadiationBox_Form', 'module', CONSOLE % ERROR )
!           call Show ( 'InitializeRadiation', 'subroutine', CONSOLE % ERROR )
!           call PROGRAM_HEADER % Abort ( )
!         end select !-- RadiationType

!         select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
!         class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
!           call RMB % Initialize &
!                  ( MS, RadiationTypeLocal, &
!                    NameShortOption = RadiationName ( iC ), &
!                    Velocity_U_UnitOption = Velocity_U_UnitOption, &
!                    MomentumDensity_U_UnitOption &
!                      = MomentumDensity_U_UnitOption, &
!                    MomentumDensity_D_UnitOption &
!                      = MomentumDensity_D_UnitOption, &
!                    EnergyDensityUnitOption = EnergyDensityUnitOption, &
!                    TemperatureUnitOption = TemperatureUnitOption, &
!                    EnergyUnitOption = EnergyUnitOption, &
!                    MomentumUnitOption = MomentumUnitOption, &
!                    AngularMomentumUnitOption = AngularMomentumUnitOption, &
!                    TimeUnitOption = TimeUnitOption )
!         end select !-- RMB
!         end select !-- MS

!         end select !--I

!     end do !-- iC

!     end select !-- I_1D

!   end subroutine InitializeRadiation


!   subroutine InitializeFluid &
!                ( RB, Velocity_U_UnitOption, MomentumDensity_D_UnitOption, &
!                  BaryonMassUnitOption, NumberDensityUnitOption, &
!                  EnergyDensityUnitOption, TemperatureUnitOption, &
!                  NumberUnitOption, EnergyUnitOption, MomentumUnitOption, &
!                  AngularMomentumUnitOption, TimeUnitOption, &
!                  BaryonMassReferenceOption )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
!       Velocity_U_UnitOption, &
!       MomentumDensity_D_UnitOption
!     type ( MeasuredValueForm ), intent ( in ), optional :: &
!       BaryonMassUnitOption, &
!       NumberDensityUnitOption, &
!       EnergyDensityUnitOption, &
!       TemperatureUnitOption, &
!       NumberUnitOption, &
!       EnergyUnitOption, &
!       MomentumUnitOption, &
!       AngularMomentumUnitOption, &
!       TimeUnitOption
!     real ( KDR ), intent ( in ), optional :: &
!       BaryonMassReferenceOption

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_C_PS_Template )

!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )

!     allocate ( Fluid_ASC_Form :: I % Current_ASC )
!     select type ( FA => I % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     call FA % Initialize &
!            ( PS, 'IDEAL', &
!              Velocity_U_UnitOption = Velocity_U_UnitOption, &
!              MomentumDensity_D_UnitOption = MomentumDensity_D_UnitOption, &
!              BaryonMassUnitOption = BaryonMassUnitOption, &
!              NumberDensityUnitOption = NumberDensityUnitOption, &
!              EnergyDensityUnitOption = EnergyDensityUnitOption, &
!              TemperatureUnitOption = TemperatureUnitOption, &
!              NumberUnitOption = NumberUnitOption, &
!              EnergyUnitOption = EnergyUnitOption, &
!              MomentumUnitOption = MomentumUnitOption, &
!              AngularMomentumUnitOption = AngularMomentumUnitOption, &
!              TimeUnitOption = TimeUnitOption, &
!              BaryonMassReferenceOption = BaryonMassReferenceOption )
!     end select !-- FA

!     end select !-- PS
!     end select !-- I

!   end subroutine InitializeFluid


!   subroutine InitializeSteps &
!                ( RB, Name, ApplyStreamingOption, ApplyInteractionsOption, &
!                  EvolveFluidOption )

!     class ( RadiationBoxForm ), intent ( inout ) :: &
!       RB
!     character ( * ), intent ( in ) :: &
!       Name
!     logical ( KDL ), intent ( in ), optional :: &
!       ApplyStreamingOption, &
!       ApplyInteractionsOption, &
!       EvolveFluidOption

!     integer ( KDI ) :: &
!       iC  !-- iCurrent


!     !-- Operators

!     RB % ApplyStreaming    = .true.
!     RB % ApplyInteractions = .true.
!     RB % EvolveFluid       = .true.
!     if ( present ( ApplyStreamingOption ) ) &
!       RB % ApplyStreaming = ApplyStreamingOption
!     if ( present ( ApplyInteractionsOption ) ) &
!       RB % ApplyInteractions = ApplyInteractionsOption
!     if ( present ( EvolveFluidOption ) ) &
!       RB % EvolveFluid = EvolveFluidOption


!     !-- Relaxation

!     if ( RB % ApplyInteractions ) then
!       select type ( I => RB % Integrator )
!       class is ( Integrator_C_1D_PS_C_PS_Form )

!         allocate ( RB % Relaxation_RM_ASC )
!         associate ( R => RB % Relaxation_RM_ASC )
!         select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
!         class is ( RadiationMoments_ASC_Form )
!           call R % Initialize ( RMA, Name = RB % Name )
!         end select !-- RMA
!         end associate !-- R

!       class is ( Integrator_C_1D_MS_C_PS_Form )

!         allocate ( RB % Relaxation_RM_BSLL_ASC_CSLD )
!         associate ( R => RB % Relaxation_RM_BSLL_ASC_CSLD )
!         select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
!         class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
!           call R % Initialize ( RMB, Name = RB % Name )
!         end select !-- RMB
!         end associate !-- R

!       end select !-- I
!     end if !-- ApplyInteractions


!     !-- Radiation step

!     select type ( I => RB % Integrator )
!     class is ( Integrator_C_1D_PS_C_PS_Form )

!       allocate ( Step_RK2_C_ASC_1D_Form :: I % Step_1D )
!       select type ( S_1D => I % Step_1D )
!       class is ( Step_RK2_C_ASC_1D_Form )
!       call S_1D % Initialize ( I, I % Current_ASC_1D, Name )

!       if ( .not. RB % ApplyStreaming ) then
!         do iC = 1, I % N_CURRENTS_1D
!           S_1D % ApplyDivergence_1D ( iC ) % Pointer  =>  null ( )  
!         end do !-- iC
!       end if

!       if ( RB % ApplyInteractions ) then
!         do iC = 1, I % N_CURRENTS_1D
!           S_1D % ApplyRelaxation_1D ( iC ) % Pointer  &
!             =>  RB % Relaxation_RM_ASC % Apply 
!         end do !-- iC
!       end if

!       end select !-- S_1D

!     class is ( Integrator_C_1D_MS_C_PS_Form )

!       allocate ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form :: I % Step_1D )
!       select type ( S_1D => I % Step_1D )
!       class is ( Step_RK2_C_BSLL_ASC_CSLD_1D_Form )
!       call S_1D % Initialize ( I, I % Current_BSLL_ASC_CSLD_1D, Name )

!       if ( .not. RB % ApplyStreaming ) then
!         do iC = 1, I % N_CURRENTS_1D
!           S_1D % ApplyDivergence_S ( iC ) % Pointer  =>  null ( )  
!         end do !-- iC
!       end if

!       if ( RB % ApplyInteractions ) then
!         do iC = 1, I % N_CURRENTS_1D
!           S_1D % ApplyRelaxation_F ( iC ) % Pointer  &
!             =>  RB % Relaxation_RM_BSLL_ASC_CSLD % Apply 
!         end do !-- iC
!       end if

!       end select !-- S_1D

!     end select !-- I


!     !-- Fluid

!     if ( RB % EvolveFluid ) then
!       select type ( I => RB % Integrator )
!       class is ( Integrator_C_1D_C_PS_Template )

!       allocate ( Step_RK2_C_ASC_Form :: I % Step )
!       select type ( S => I % Step )
!       class is ( Step_RK2_C_ASC_Form )
!         call S % Initialize ( I, I % Current_ASC, 'Fluid' )
!         if ( RB % ApplyInteractions ) &
!           S % ApplySources % Pointer  =>  ApplySources_Fluid
!       end select !-- S

!       I % PrepareStep => PrepareStep

!       end select !-- I
!     end if !-- EvolveFluid

!   end subroutine InitializeSteps


!   subroutine ComputeTimeStepLocal ( I, TimeStepCandidate )

!     class ( IntegratorTemplate ), intent ( inout ), target :: &
!       I
!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       TimeStepCandidate

!     integer ( KDI ) :: &
!       oC  !-- oCandidate

!     associate ( RB => RadiationBox )

!     select type ( I )
!     class is ( Integrator_C_1D_C_PS_Template )

!     !-- Fluid advection and Radiation streaming

!     call I % ComputeTimeStepLocalTemplate ( TimeStepCandidate )

!     if ( .not. RB % ApplyStreaming ) &
!       TimeStepCandidate ( 1 : I % N_CURRENTS_1D )  =  huge ( 1.0_KDR )

!     if ( .not. RB % EvolveFluid ) &
!       TimeStepCandidate ( I % N_CURRENTS_1D  +  1 )  =  huge ( 1.0_KDR )

!     !-- Interactions

!     oC  =  I % N_CURRENTS_1D  +  1
!     if ( RB % ApplyInteractions .and. RB % EvolveFluid ) then
!       call ComputeTimeStepInteractions &
!              ( I, TimeStepCandidate ( oC + 1 : oC + I % N_CURRENTS_1D ) )
!     else
!       TimeStepCandidate ( oC + 1 : oC + I % N_CURRENTS_1D )  &
!         =  huge ( 1.0_KDR )
!     end if

!     end select !-- I
!     end associate !-- RB

!   end subroutine ComputeTimeStepLocal


!   subroutine PrepareStep ( I )

!     class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
!       I

!     class ( Fluid_P_I_Form ), pointer :: &
!       F
!     class ( RadiationMomentsForm ), pointer :: &
!       RM
!     class ( Sources_RM_Form ), pointer :: &
!       SR

!     ! associate &
!     !   ( R => Radiation, &
!     !     F => Fluid )

!     select type ( I )
!     class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

!       select type ( RMA => I % Current_ASC_1D ( 1 ) % Element )
!       class is ( RadiationMoments_ASC_Form )
!         RM => RMA % RadiationMoments ( )
!       end select !-- RMA

!     class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

!       call IntegrateSources ( I )

!       select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
!       class is ( RadiationMoments_BSLL_ASC_CSLD_Form )
!       select type ( RMA => RMB % EnergyIntegral )
!       class is ( RadiationMoments_ASC_Form )
!         RM => RMA % RadiationMoments ( )
!       end select !-- RMA
!       end select !-- RMB

!     end select !-- I

!     select type ( FA => I % Current_ASC )
!     class is ( Fluid_ASC_Form )
!       F => FA % Fluid_P_I ( )
!     end select !-- FA

!     select type ( SR => RM % Sources )
!     class is ( Sources_RM_Form )

!     select type ( SF => F % Sources )
!     class is ( Sources_F_Form )

!     select type ( PS => I % PositionSpace )
!     class is ( Atlas_SC_Form )

!     select type ( PSC => PS % Chart )
!     class is ( Chart_SL_Template )

!     call Clear ( SF % Value )

!     call ComputeFluidSource_G_S_Radiation_Kernel &
!            ( SF % Value ( :, SF % RADIATION_G ), & 
!              SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
!              SF % Value ( :, SF % RADIATION_S_D ( 2 ) ), &
!              SF % Value ( :, SF % RADIATION_S_D ( 3 ) ), &
!              PSC % IsProperCell, &
!              SR % Value ( :, SR % INTERACTIONS_J ), &
!              SR % Value ( :, SR % INTERACTIONS_H_D ( 1 ) ), &
!              SR % Value ( :, SR % INTERACTIONS_H_D ( 2 ) ), &
!              SR % Value ( :, SR % INTERACTIONS_H_D ( 3 ) ) )

!     end select !-- PSC
!     end select !-- PS
!     end select !-- SF
!     end select !-- SR
!     nullify ( F, RM, SR )

!   end subroutine PrepareStep


!   subroutine IntegrateSources ( I )

!     type ( Integrator_C_1D_MS_C_PS_Form ), intent ( inout ), target :: &
!       I

!     integer ( KDI ) :: &
!       iI, &  !-- iIntegral
!       iF, &  !-- iFiber
!       nIntegrals
!     real ( KDR ), dimension ( : ), allocatable :: &
!       Integral
!     type ( Real_1D_Form ), dimension ( : ), allocatable :: &
!       Integrand
!     type ( VolumeIntegralForm ) :: &
!       VI
!     class ( RadiationMomentsForm ), pointer :: &
!       RMEI, &
!       RMF

!     associate ( IB => RadiationBox % Interactions_BSLL_ASC_CSLD ) 

!     select type ( RMB => I % Current_BSLL_ASC_CSLD_1D ( 1 ) % Element )
!     class is ( RadiationMoments_BSLL_ASC_CSLD_Form )

!     select type ( RMA => RMB % EnergyIntegral )
!     class is ( RadiationMoments_ASC_Form )

!     RMEI => RMA % RadiationMoments ( )

!     nIntegrals = 4

!     allocate ( Integral  ( nIntegrals ) )
!     allocate ( Integrand ( nIntegrals ) )
!     do iI = 1, nIntegrals
!       call Integrand ( iI ) % Initialize ( RMB % nEnergyValues )
!     end do !-- iC

!     associate ( MS => RMB % Bundle_SLL_ASC_CSLD )
!     do iF = 1, MS % nFibers
!       associate &
!         ( iBC => MS % iaBaseCell ( iF ), &
!           CF  => MS % Fiber_CSLL )

!       RMF => RMB % RadiationMoments ( iF )
!       select type ( SF => RMF % Sources )
!       class is ( Sources_RM_Form )
!         Integrand ( 1 ) % Value  &
!           =  SF % Value ( :, SF % INTERACTIONS_J )
!         Integrand ( 2 ) % Value  &
!           =  SF % Value ( :, SF % INTERACTIONS_H_D ( 1 ) )
!         Integrand ( 3 ) % Value  &
!           =  SF % Value ( :, SF % INTERACTIONS_H_D ( 2 ) )
!         Integrand ( 4 ) % Value  &
!           =  SF % Value ( :, SF % INTERACTIONS_H_D ( 3 ) )
!       end select !-- SF

!       call VI % Compute ( CF, Integrand, Integral )

!       select type ( SEI => RMEI % Sources )
!       class is ( Sources_RM_Form )
!         SEI % Value ( iBC, SEI % INTERACTIONS_J )         = Integral ( 1 ) 
!         SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 1 ) ) = Integral ( 2 ) 
!         SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 2 ) ) = Integral ( 3 ) 
!         SEI % Value ( iBC, SEI % INTERACTIONS_H_D ( 3 ) ) = Integral ( 4 ) 
!       end select !-- SEI
      
!       end associate !-- iBC, etc.
!     end do !-- iF
!     end associate !-- MS

!     end select !-- RMA
!     end select !-- RMB
!     end associate !-- IB

!     nullify ( RMEI, RMF )
 
!   end subroutine IntegrateSources


!   subroutine ApplySources_Fluid &
!                ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

!     class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
!       S
!     class ( Sources_C_Form ), intent ( inout ) :: &
!       Sources_F
!     type ( StorageForm ), intent ( inout ), target :: &
!       Increment
!     class ( CurrentTemplate ), intent ( in ) :: &
!       Fluid
!     real ( KDR ), intent ( in ) :: &
!       TimeStep
!     integer ( KDI ), intent ( in ) :: &
!       iStage

!     integer ( KDI ) :: &
!       iEnergy, &
!       iMomentum_1, &
!       iMomentum_2, &
!       iMomentum_3

!     select type ( SF => Sources_F )
!     class is ( Sources_F_Form )

!     select type ( F => Fluid )
!     class is ( Fluid_P_I_Form )

!     select type ( Chart => S % Chart )
!     class is ( Chart_SL_Template )

!     call Search ( F % iaConserved, F % CONSERVED_ENERGY, iEnergy )
!     call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
!     call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
!     call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 3 ), iMomentum_3 )

!     call ApplySources_Fluid_Kernel &
!            ( Increment % Value ( :, iEnergy ), &
!              Increment % Value ( :, iMomentum_1 ), &
!              Increment % Value ( :, iMomentum_2 ), &
!              Increment % Value ( :, iMomentum_3 ), &
!              Chart % IsProperCell, &
!              SF % Value ( :, SF % RADIATION_G ), &
!              SF % Value ( :, SF % RADIATION_S_D ( 1 ) ), &
!              SF % Value ( :, SF % RADIATION_S_D ( 2 ) ), &
!              SF % Value ( :, SF % RADIATION_S_D ( 3 ) ), &
!              TimeStep )

!     end select !-- Chart
!     end select !-- F
!     end select !-- SF

!   end subroutine ApplySources_Fluid


!   subroutine ComputeTimeStepInteractions ( I, TimeStepCandidate )

!     class ( Integrator_C_1D_C_PS_Template ), intent ( inout ) :: &
!       I
!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       TimeStepCandidate

!     integer ( KDI ) :: &
!       iC     !-- iCurrent

!     associate ( RB => RadiationBox )

!     select type ( I )
!     class is ( Integrator_C_1D_PS_C_PS_Form )  !-- Grey

!       associate ( IA => RB % Interactions_ASC )  
!       do iC = 1, I % N_CURRENTS_1D  
!         associate ( RMA => I % Current_ASC_1D ( iC ) % Element )
!         call IA % ComputeTimeScale ( RMA, TimeStepCandidate ( iC ) )
!         TimeStepCandidate ( iC )  &
!           =  RB % InteractionFactor  *  TimeStepCandidate ( iC )
!         end associate !-- RA
!       end do !-- iC
!       end associate !-- IA

!     class is ( Integrator_C_1D_MS_C_PS_Form )  !-- Spectral

!       associate ( IB => RB % Interactions_BSLL_ASC_CSLD )
!       do iC = 1, I % N_CURRENTS_1D  
!         associate ( RMB => I % Current_BSLL_ASC_CSLD_1D ( iC ) % Element )
!         call IB % ComputeTimeScale ( RMB, TimeStepCandidate ( iC ) )
!         TimeStepCandidate ( iC )  &
!           =  RB % InteractionFactor  *  TimeStepCandidate ( iC )
!         end associate !-- RB
!       end do !-- iC
!       end associate !-- IB

!     end select !-- I
!     end associate !-- RB

!   end subroutine ComputeTimeStepInteractions


!   subroutine ComputeFluidSource_G_S_Radiation_Kernel &
!                ( FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, IsProperCell, &
!                  Q, A_1, A_2, A_3 )

!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       FS_R_G, &
!       FS_R_S_1, FS_R_S_2, FS_R_S_3
!     logical ( KDL ), dimension ( : ), intent ( in ) :: &
!       IsProperCell
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       Q, &
!       A_1, A_2, A_3

!     integer ( KDI ) :: &
!       iV, &
!       nV

!     nV = size ( FS_R_G )

!     !$OMP parallel do private ( iV )
!     do iV = 1, nV
!       if ( .not. IsProperCell ( iV ) ) &
!         cycle
!       FS_R_G ( iV )    =    FS_R_G ( iV )  -    Q ( iV ) 
!       FS_R_S_1 ( iV )  =  FS_R_S_1 ( iV )  -  A_1 ( iV )
!       FS_R_S_2 ( iV )  =  FS_R_S_2 ( iV )  -  A_2 ( iV )
!       FS_R_S_3 ( iV )  =  FS_R_S_3 ( iV )  -  A_3 ( iV )
!     end do
!     !$OMP end parallel do

!   end subroutine ComputeFluidSource_G_S_Radiation_Kernel


!   subroutine ApplySources_Fluid_Kernel &
!                ( K_G, K_S_1, K_S_2, K_S_3, IsProperCell, &
!                  FS_R_G, FS_R_S_1, FS_R_S_2, FS_R_S_3, dt )

!     real ( KDR ), dimension ( : ), intent ( inout ) :: &
!       K_G, &
!       K_S_1, K_S_2, K_S_3
!     logical ( KDL ), dimension ( : ), intent ( in ) :: &
!       IsProperCell
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       FS_R_G, &
!       FS_R_S_1, FS_R_S_2, FS_R_S_3
!     real ( KDR ), intent ( in ) :: &
!       dt

!     integer ( KDI ) :: &
!       iV, &
!       nV

!     nV = size ( K_G )

!     !$OMP parallel do private ( iV )
!     do iV = 1, nV
!       if ( .not. IsProperCell ( iV ) ) &
!         cycle
!       K_G ( iV )    =  K_G ( iV )    +  FS_R_G ( iV )    *  dt
!       K_S_1 ( iV )  =  K_S_1 ( iV )  +  FS_R_S_1 ( iV )  *  dt
!       K_S_2 ( iV )  =  K_S_2 ( iV )  +  FS_R_S_2 ( iV )  *  dt
!       K_S_3 ( iV )  =  K_S_3 ( iV )  +  FS_R_S_3 ( iV )  *  dt
!     end do
!     !$OMP end parallel do

!   end subroutine ApplySources_Fluid_Kernel


end module RadiationBox_Form
