program Sources_F__Form_Test
  
  use Basics
  use Fluid_D__Form
  use Sources_F__Form

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Sources_F__Form_Test'

  integer ( KDI ) :: &
    iD
  type ( Fluid_D_Form ), allocatable :: &
    F
  type ( Sources_F_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )

  call CONSOLE % SetVerbosity ( 'INFO_4' )

  allocate ( F, S )
  call F % Initialize &
         ( RiemannSolverType = 'HLL', UseLimiter = .true., &
           Velocity_U_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           MomentumDensity_D_Unit = [ ( UNIT % IDENTITY, iD = 1, 3 ) ], &
           BaryonMassUnit = UNIT % IDENTITY, &
           NumberDensityUnit = UNIT % IDENTITY, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, nValues = 10 )
  call F % SetPrimitiveConserved ( )
  call S % Initialize ( F, TimeUnit = UNIT % IDENTITY )

  call Show ( 'Sources_F Variables', S % IGNORABILITY )
  call Show ( S % Variable, 'Variable', S % IGNORABILITY )

  deallocate ( S, F )
  deallocate ( PROGRAM_HEADER )

end program Sources_F__Form_Test
