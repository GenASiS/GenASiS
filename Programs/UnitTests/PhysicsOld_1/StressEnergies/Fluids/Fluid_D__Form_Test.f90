program Fluid_D__Form_Test

  use Basics
  use Fluid_D__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_D__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Fluid_D_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )

  call CONSOLE % SetVerbosity ( 'INFO_4' )

  allocate ( F )
  call F % Initialize &
         ( RiemannSolverType = 'HLL', UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           BaryonMassUnit = UNIT % IDENTITY, &
           NumberDensityUnit = UNIT % IDENTITY, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, nValues = 10 )

  call Show ( 'Fluid_D Variables', F % IGNORABILITY )
  call Show ( F % Variable, 'Variable', F % IGNORABILITY )

  deallocate ( F )
  deallocate ( PROGRAM_HEADER )

end program Fluid_D__Form_Test
