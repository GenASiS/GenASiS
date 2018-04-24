program FluidFeatures_CSL__Form_Test

  use Basics
  use Mathematics
  use Spaces
  use FluidFeatures_CSL__Form
  use Fluid_CSL__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'FluidFeatures_CSL__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Atlas_SC_Form ), allocatable :: &
    A
  type ( Geometry_ASC_Form ), allocatable :: &
    GA
  type ( Fluid_CSL_Form ), allocatable :: &
    FC
  type ( FluidFeatures_CSL_Form ), allocatable :: &
    FFC

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  allocate ( A )
  call A % Initialize ( 'Atlas_SC', PROGRAM_HEADER % Communicator )
  call A % CreateChart ( )

  allocate ( GA )
  call GA % Initialize ( A, 'GALILEAN' )
  call A % SetGeometry ( GA )

  select type ( C => A % Chart )
  class is ( Chart_SL_Template )
  associate ( nValues => C % nProperCells + C % nGhostCells )

  allocate ( FC )
  call FC % Initialize &
         ( C, 'Fluid', 'IDEAL', RiemannSolverType = 'HLLC', &
           UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           BaryonMassUnit = UNIT % IDENTITY, &
           NumberDensityUnit = UNIT % IDENTITY, &
           EnergyDensityUnit = UNIT % IDENTITY, &
           TemperatureUnit = UNIT % IDENTITY, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, nValues = nValues )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  allocate ( FFC )

  call FFC % Initialize &
         ( FC, 'FluidFeatures', FC % FluidType, ShockThreshold = 0.1_KDR, &
           nValues = nValues )
  call FC % SetFeatures ( FFC )

  deallocate ( FFC )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( FC )
  end associate !-- nValues
  end select !-- C
  deallocate ( GA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program FluidFeatures_CSL__Form_Test
