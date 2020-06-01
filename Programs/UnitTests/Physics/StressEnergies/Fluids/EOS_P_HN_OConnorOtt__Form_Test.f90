program EOS_P_HN_OConnorOtt__Form_Test

  use Basics
  use StressEnergyUnits_Form
  use Fluid_P_HN__Form
  use EOS_P_HN_OConnorOtt__Form
  
  implicit none
  
  character ( LDF ) :: &
    ProgramName = 'EOS_P_HN_OConnorOtt_Form_Test'

  integer ( KDI ) :: &
    iD, &
    iV, &
    oV, &
    nV
  integer ( KDI ), dimension ( 3 ) :: &
    nCells
  integer ( KDI ), dimension ( : ), allocatable :: &
    iaFluidOutput, &
    iaSelected_EOS
  type ( TimerForm ) :: &
    Timer
  type ( StressEnergyUnitsForm ) :: &
    SU
  type ( Fluid_P_HN_Form ), allocatable :: &
    F
  type ( EOS_P_HN_OConnorOtt_Form ), allocatable :: &
    EOS 

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )
  
  nCells = 1
  call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  
  call Timer % Initialize ( 'Timer', Level = 1 )
  
  oV = 0
  nV = product ( nCells )
  
  allocate ( F )
  call F % Initialize_P_HN &
         ( FluidType = 'HEAVY_NUCLEUS', &
           RiemannSolverType = 'HLLC', &
           ReconstructedType = 'PRIMITIVE', &
           UseEntropy = .false., UseLimiter = .true., &
           Units = SU, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, &
           nValues = nV )
  
  call F % AllocateDevice ( AssociateVariablesOption = .false. )
  
  associate ( FV => F % Value )
  associate &
    ( FEP_1 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 1 ) ), & 
      FEP_2 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 2 ) ), & 
      FEP_3 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_PLUS ( 3 ) ), & 
      FEM_1 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 1 ) ), &
      FEM_2 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 2 ) ), &
      FEM_3 => FV ( oV + 1 : oV + nV, F % FAST_EIGENSPEED_MINUS ( 3 ) ), &
      M     => FV ( oV + 1 : oV + nV, F % BARYON_MASS ), &
      N     => FV ( oV + 1 : oV + nV, F % COMOVING_BARYON_DENSITY ), &
      V_1   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 1 ) ), &
      V_2   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 2 ) ), &
      V_3   => FV ( oV + 1 : oV + nV, F % VELOCITY_U ( 3 ) ), &
      D     => FV ( oV + 1 : oV + nV, F % CONSERVED_BARYON_DENSITY ), &
      S_1   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 1 ) ), &
      S_2   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 2 ) ), &
      S_3   => FV ( oV + 1 : oV + nV, F % MOMENTUM_DENSITY_D ( 3 ) ), &
      E     => FV ( oV + 1 : oV + nV, F % INTERNAL_ENERGY ), & 
      G     => FV ( oV + 1 : oV + nV, F % CONSERVED_ENERGY ), &
      P     => FV ( oV + 1 : oV + nV, F % PRESSURE ), &   
      T     => FV ( oV + 1 : oV + nV, F % TEMPERATURE ), &
      SB    => FV ( oV + 1 : oV + nV, F % ENTROPY_PER_BARYON ), &
      DS    => FV ( oV + 1 : oV + nV, F % CONSERVED_ENTROPY ), &
  !    Gamma => FV ( oV + 1 : oV + nV, F % ADIABATIC_INDEX ), &
      CS    => FV ( oV + 1 : oV + nV, F % SOUND_SPEED ), &
      MN    => FV ( oV + 1 : oV + nV, F % MACH_NUMBER ), &
      YE    => FV ( oV + 1 : oV + nV, F % ELECTRON_FRACTION ), &
      DE    => FV ( oV + 1 : oV + nV, F % CONSERVED_ELECTRON_DENSITY ), &
      X_P   => FV ( oV + 1 : oV + nV, F % MASS_FRACTION_PROTON ), & 
      X_N   => FV ( oV + 1 : oV + nV, F % MASS_FRACTION_NEUTRON ), &
      X_He  => FV ( oV + 1 : oV + nV, F % MASS_FRACTION_ALPHA ), &
      X_A   => FV ( oV + 1 : oV + nV, F % MASS_FRACTION_HEAVY ), &
      Z     => FV ( oV + 1 : oV + nV, F % ATOMIC_NUMBER_HEAVY ), &
      A     => FV ( oV + 1 : oV + nV, F % MASS_NUMBER_HEAVY ), &
      Mu_NP => FV ( oV + 1 : oV + nV, F % CHEMICAL_POTENTIAL_N_P ), &
      Mu_E  => FV ( oV + 1 : oV + nV, F % CHEMICAL_POTENTIAL_E ), &
      U_V   => FV ( oV + 1 : oV + nV, F % UNUSED_VARIABLE ) )
  
  allocate ( EOS )
  
  call EOS % Initialize ( )
  
  allocate ( iaFluidOutput ( 12 ) )
  allocate ( iaSelected_EOS ( 12 ) )
  
  iaFluidOutput &
    = [ F % INTERNAL_ENERGY, &
        F % PRESSURE, &
        F % ENTROPY_PER_BARYON, &
        F % SOUND_SPEED, &
        F % MASS_FRACTION_ALPHA, &
        F % MASS_FRACTION_HEAVY, &
        F % MASS_FRACTION_NEUTRON, &
        F % MASS_FRACTION_PROTON, &
        F % MASS_NUMBER_HEAVY, &
        F % ATOMIC_NUMBER_HEAVY, &
        F % CHEMICAL_POTENTIAL_E, &
        F % CHEMICAL_POTENTIAL_N_P ]
  
  iaSelected_EOS &
    = [ EOS % LOG_ENERGY, &
        EOS % LOG_PRESSURE, &
        EOS % ENTROPY, &
        EOS % SOUND_SPEED_SQUARE, &
        EOS % MASS_FRACTION_A, &
        EOS % MASS_FRACTION_H, &
        EOS % MASS_FRACTION_N, &
        EOS % MASS_FRACTION_P, &
        EOS % MASS_NUMBER_BAR, &
        EOS % ATOMIC_NUMBER_BAR, &
        EOS % CHEMICAL_POTENTIAL_E, &
        EOS % CHEMICAL_POTENTIAL_HAT ]
  
  call Show ( iaFluidOutput, 'iaFluidOutput' )
  call Show ( iaSelected_EOS, 'iaSelected_EOS' )
  call EOS % SelectVariables ( iaFluidOutput, iaSelected_EOS )
  
  call EOS % AllocateDevice ( )
   
  associate &
    ( T_EOS   => EOS % Table, &
      T_L_D   => EOS % LogDensity, &
      T_L_T   => EOS % LogTemperature, &
      T_YE    => EOS % ElectronFraction, &
      Error   => EOS % Error )
  
  call Timer % Start ( )
  call Show ( 'Setting Initial Values' )
  !E  = -8.5871946287513969E+018_KDR
  E  = 0.0_KDR
  N  = 12144578686.985090_KDR
  T  = 0.61009347651875323_KDR
  P  = 0.0_KDR
  YE = 0.43110487829424210
  call Timer % Stop ( )
  
  call Timer % ShowInterval (  )
  
  call Show ( 'Input' )
  call Show ( [ N ( 1 ), T ( 1 ), YE ( 1 ) ], 'N, T, YE' )
  
  call Show ( 'ComputeFromTemperature' )
  call Timer % Start ( )
  call EOS % ComputeFromTemperature &
        ( F, iaFluidInput = [ F % COMOVING_BARYON_DENSITY, &
                              F % TEMPERATURE, F % ELECTRON_FRACTION ] )
  call Timer % Stop ( )
  call Timer % ShowInterval (  )
  
  call Show ( 'Output' )
  call Show ( [ N ( 1 ), T ( 1 ), YE ( 1 ) ], 'N, T, YE' )
  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
  
  !-- Change temperature by 10%
  T = 1.10 * T 
  
  call Show ( 'ComputeFromEnergy' )
  
  call Show ( 'Input' )
  call Show ( [ N ( 1 ), T ( 1 ), YE ( 1 ) ], 'N, T, YE' )
  
  call Timer % Start ( )
  call EOS % ComputeFromEnergy &
        ( F, iaFluidInput = [ F % COMOVING_BARYON_DENSITY, &
                              F % TEMPERATURE, F % ELECTRON_FRACTION ], &
          iSolve = F % INTERNAL_ENERGY )
  call Timer % Stop ( )
  call Timer % ShowInterval (  )
  
  call Show ( 'Output' )
  call Show ( [ N ( 1 ), T ( 1 ), YE ( 1 ) ], 'N, T, YE' )
  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
  
  call CONSOLE % Mute ( )
  call Show ( 'Fluid_P_HN Variables', F % IGNORABILITY )
  call Show ( F % Variable, 'Variable', F % IGNORABILITY )
  
  end associate

  end associate
  end associate

  deallocate ( F )
  deallocate ( PROGRAM_HEADER )


end program EOS_P_HN_OConnorOtt__Form_Test
