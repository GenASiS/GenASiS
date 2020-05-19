#include "Preprocessor"

program Fluid_P_HN__Form_Test

  use Basics
  use StressEnergyUnits_Form
  use Fluid_P_HN__Form
  
  use NUC_EOS

  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_P_HN__Form_Test'

  integer ( KDI ) :: &
    iD, &
    iV, &
    oV, &
    nV
  integer ( KDI ), dimension ( 3 ) :: &
    nCells
  type ( TimerForm ) :: &
    Timer
  type ( StressEnergyUnitsForm ) :: &
    SU
  type ( Fluid_P_HN_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )
  
  nCells = 32
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
  
  call F % AllocateDevice ( )
  
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
  
  associate &
    ( T_EOS   => F % EOS % Table, &
      T_L_D   => F % EOS % LogDensity, &
      T_L_T   => F % EOS % LogTemperature, &
      T_YE    => F % EOS % ElectronFraction, &
      Error   => F % EOS % Error )
  
  call Timer % Start ( )
  call Show ( 'Setting Initial Values' )
  E  = -8.5871946287513969E+018_KDR
  N  = 12144578686.985090_KDR
  T  = 0.61009347651875323_KDR
  P  = 0.0_KDR
  YE = 0.43110487829424210
  call Timer % Stop ( )
  
  call Timer % ShowInterval (  )

  call Timer % Start ( )
  call F % UpdateDevice ( )
  call Timer % Stop ( )
  call Show ( 'Update Device' )
  call Timer % ShowInterval (  )

  
  call Timer % Start ( )
  call Show ( 'Offload OpenMP NUC_EOS' )
  !$OMP  OMP_TARGET_DIRECTIVE parallel do &
  !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
  !$OMP& private ( iV )
  do iV = 1, nV
    call NUC_EOS_FULL &
         ( N ( iV ), T ( iV ), YE ( iV ), E ( iV ), &
           P ( iV ), SB ( iV ), CS ( iV ), U_V ( iV ), &
           U_V ( iV ), U_V ( iV ), X_He ( iV ), X_A ( iV ), &
           X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
           U_V ( iV ), U_V ( iV ), Mu_NP ( iV ), EOS_Apply_EOS_HN_T, &
           Error ( iV ), EOS_RF_Accuracy, T_L_D, T_L_T, T_YE, T_EOS )
  end do
  call Timer % Stop ( )
  call Timer % ShowInterval (  )

  print*, E ( 1 ), N ( 1 ), T ( 1 ), P ( 1 ), YE ( 1 ), &
          CS ( 1 ), X_He ( 1 )
  
  call Timer % Start ( )
  call Show ( 'CPU OpenMP NUC_EOS' )
  !$OMP  parallel do &
  !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
  !$OMP& private ( iV )
  do iV = 1, nV
    call NUC_EOS_FULL &
         ( N ( iV ), T ( iV ), YE ( iV ), E ( iV ), &
           P ( iV ), SB ( iV ), CS ( iV ), U_V ( iV ), &
           U_V ( iV ), U_V ( iV ), X_He ( iV ), X_A ( iV ), &
           X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
           U_V ( iV ), U_V ( iV ), Mu_NP ( iV ), EOS_Apply_EOS_HN_T, &
           Error ( iV ), EOS_RF_Accuracy, T_L_D, T_L_T, T_YE, T_EOS )
  end do
  call Timer % Stop ( )
  call Timer % ShowInterval (  )
  
  call CONSOLE % Mute ( )

  call Show ( 'Fluid_P_HN Variables', F % IGNORABILITY )
  call Show ( F % Variable, 'Variable', F % IGNORABILITY )
  
  end associate
  end associate
  end associate

  deallocate ( F )
  deallocate ( PROGRAM_HEADER )

end program Fluid_P_HN__Form_Test
