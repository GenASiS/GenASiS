program FluidFeatures_P__Form_Test

  use Basics
  use Mathematics
  use Fluid_P_I__Form
  use FluidFeatures_P__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'FluidFeatures_P__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Chart_SLD_Form ) :: &
    C
  type ( Fluid_P_I_Form ), allocatable :: &
    F
  type ( FluidFeatures_P_Form ), allocatable :: &
    FF

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )

  call CONSOLE % SetVerbosity ( 'INFO_4' )

  allocate ( F, FF )
  call F % Initialize &
         ( RiemannSolverType = 'HLLC', &
           UseEntropy = .false., UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           BaryonMassUnit = UNIT % IDENTITY, &
           NumberDensityUnit = UNIT % IDENTITY, &
           EnergyDensityUnit = UNIT % IDENTITY, &
           TemperatureUnit = UNIT % IDENTITY, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, nValues = 10 )
  call FF % Initialize &
         ( F, C, ShockThreshold = 0.1_KDR, nValues = 10 )

  deallocate ( FF, F )
  deallocate ( PROGRAM_HEADER )

end program FluidFeatures_P__Form_Test
