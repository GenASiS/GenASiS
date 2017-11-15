program Fluid_D__Form_Test

  use Basics
  use Fluid_D__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_D__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Fluid_D_Form ) :: &
    F

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize ( ProgramName )

  call F % Initialize &
         ( RiemannSolverType = 'HLL', UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           NumberDensityUnit = UNIT % IDENTITY, &
           LimiterParameter = 1.4_KDR, nValues = 10 )

  deallocate ( PROGRAM_HEADER )

end program Fluid_D__Form_Test