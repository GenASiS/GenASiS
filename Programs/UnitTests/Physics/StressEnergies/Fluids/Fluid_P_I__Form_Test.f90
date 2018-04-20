program Fluid_P_I__Form_Test

  use Basics
  use Fluid_P_I__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_P_I__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Fluid_P_I_Form ) :: &
    F

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )

  call F % Initialize &
         ( RiemannSolverType = 'HLLC', UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           BaryonMassUnit = UNIT % IDENTITY, &
           NumberDensityUnit = UNIT % IDENTITY, &
           EnergyDensityUnit = UNIT % IDENTITY, &
           TemperatureUnit = UNIT % IDENTITY, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, nValues = 10 )

  deallocate ( PROGRAM_HEADER )

end program Fluid_P_I__Form_Test
