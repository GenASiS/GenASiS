program Fluid_P_HN__Form_Test

  use Basics
  use Mathematics
  use StressEnergyUnits_Form
  use FluidFeatures_P__Form
  use Fluid_P_HN__Form
  
  implicit none

  character ( LDF ) :: &
    ProgramName = 'Fluid_P_HN__Form_Test'
    
  integer ( KDI ), parameter :: &
    iRADIUS_TS            = 2, &  !-- must match the profile file columns
    iRADIAL_VELOCITY_TS   = 3, &
    iDENSITY_TS           = 4, &
    iTEMPERATURE_TS       = 5, &
    iSPECIFIC_ENERGY_TS   = 10, &
    iELECTRON_FRACTION_TS = 11

  integer ( KDI ) :: &
    iD, &
    iV, &
    iP
  integer ( KDI ), dimension ( 3 ) :: &
    nCells
  real ( KDR ), dimension ( : ), allocatable :: &
    T_Copy
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Profile
  character ( LDF ) :: &
    Filename, &
    Path
  type ( TimerForm ) :: &
    Timer
  type ( MeasuredValueForm ), dimension ( 3 ) :: &
    CoordinateUnit
  type ( TableStreamForm ) :: &
    TS  
  type ( GeometryFlatForm ) :: &
    Dummy_G
  type ( Chart_SLD_Form ) :: &
    C
  type ( StressEnergyUnitsForm ) :: &
    SU
  type ( FluidFeatures_P_Form ), allocatable :: &
    FF
  type ( Fluid_P_HN_Form ), allocatable :: &
    F
    
  allocate ( PROGRAM_HEADER )  
  call PROGRAM_HEADER % Initialize &
         ( ProgramName, AppendDimensionalityOption = .false. )
  
  nCells = 1
  call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

  call CONSOLE % SetVerbosity ( 'INFO_4' )
  
  call Timer % Initialize ( 'Timer', Level = 1 )
  
  Path = '../../../../../Applications/WoosleyHeger_07/Parameters/'
  Filename = 'WH07_S12_08.d.stripped'
  call TS % Initialize &
         ( Filename, PROGRAM_HEADER % Communicator % Rank, &
           PathOption = Path )
  call TS % Read ( Profile, nRowsOption = 300, oRowOption = 2 )
  
  allocate ( F, FF )
  call F % Initialize_P_HN &
         ( FluidType = 'HEAVY_NUCLEUS', &
           RiemannSolverType = 'HLLC', &
           ReconstructedType = 'PRIMITIVE', &
           UseEntropy = .true., UseLimiter = .true., &
           Units = SU, &
           BaryonMassReference = 1.0_KDR, &
           LimiterParameter = 1.4_KDR, &
           UseInitialTemperature = .true., &
           nValues = product ( nCells ) )
  
  call FF % Initialize &
         ( F, C, ShockThreshold = 0.1_KDR, nValues = product ( nCells ) )
  FF % Value = 2.0_KDR  !-- Set to no shock

  call F % SetFeatures ( FF )
           
  call Dummy_G % Initialize &
         ( 'RECTANGULAR', CoordinateUnit, product ( nCells ) )
  
  call F % AllocateDevice ( )
  
  associate ( FV => F % Value )
  associate &
    ( FEP_1 => FV ( :, F % FAST_EIGENSPEED_PLUS ( 1 ) ), & 
      FEP_2 => FV ( :, F % FAST_EIGENSPEED_PLUS ( 2 ) ), & 
      FEP_3 => FV ( :, F % FAST_EIGENSPEED_PLUS ( 3 ) ), & 
      FEM_1 => FV ( :, F % FAST_EIGENSPEED_MINUS ( 1 ) ), &
      FEM_2 => FV ( :, F % FAST_EIGENSPEED_MINUS ( 2 ) ), &
      FEM_3 => FV ( :, F % FAST_EIGENSPEED_MINUS ( 3 ) ), &
      M     => FV ( :, F % BARYON_MASS ), &
      N     => FV ( :, F % COMOVING_BARYON_DENSITY ), &
      V_1   => FV ( :, F % VELOCITY_U ( 1 ) ), &
      V_2   => FV ( :, F % VELOCITY_U ( 2 ) ), &
      V_3   => FV ( :, F % VELOCITY_U ( 3 ) ), &
      D     => FV ( :, F % CONSERVED_BARYON_DENSITY ), &
      S_1   => FV ( :, F % MOMENTUM_DENSITY_D ( 1 ) ), &
      S_2   => FV ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), &
      S_3   => FV ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), &
      E     => FV ( :, F % INTERNAL_ENERGY ), & 
      G     => FV ( :, F % CONSERVED_ENERGY ), &
      P     => FV ( :, F % PRESSURE ), &   
      T     => FV ( :, F % TEMPERATURE ), &
      SB    => FV ( :, F % ENTROPY_PER_BARYON ), &
      DS    => FV ( :, F % CONSERVED_ENTROPY ), &
      CS    => FV ( :, F % SOUND_SPEED ), &
      MN    => FV ( :, F % MACH_NUMBER ), &
      Y_E   => FV ( :, F % ELECTRON_FRACTION ), &
      DE    => FV ( :, F % CONSERVED_ELECTRON_DENSITY ), &
      X_P   => FV ( :, F % MASS_FRACTION_PROTON ), & 
      X_N   => FV ( :, F % MASS_FRACTION_NEUTRON ), &
      X_He  => FV ( :, F % MASS_FRACTION_ALPHA ), &
      X_A   => FV ( :, F % MASS_FRACTION_HEAVY ), &
      Z     => FV ( :, F % ATOMIC_NUMBER_HEAVY ), &
      A     => FV ( :, F % MASS_NUMBER_HEAVY ), &
      Mu_NP => FV ( :, F % CHEMICAL_POTENTIAL_N_P ), &
      Mu_E  => FV ( :, F % CHEMICAL_POTENTIAL_E ), &
      U_V   => FV ( :, F % UNUSED_VARIABLE ) )
  
  call Timer % Start ( )
  call Show ( 'Setting Initial Values' )

  iP = 1  
  do iV = 1, size ( N )
    N ( iV )    = Profile ( iP, iDENSITY_TS ) / CONSTANT % ATOMIC_MASS_UNIT
    E ( iV )    = Profile ( iP, iSPECIFIC_ENERGY_TS ) &
                   * Profile ( iP, iDENSITY_TS )
    T ( iV )    = Profile ( iP, iTEMPERATURE_TS )
    Y_E ( iV )  = Profile ( iP, iELECTRON_FRACTION_TS )
    
    iP = iP + 1
    if ( iP > size ( Profile, dim = 1 ) ) &
      iP = 1
  end do
  
  !-- Set unit
  N = N * UNIT % MASS_DENSITY_CGS
  T = T * UNIT % KELVIN
  E = E * ( UNIT % ERG / UNIT % GRAM ) * UNIT % MASS_DENSITY_CGS
  
  V_1 = 0.0_KDR
  V_2 = 0.0_KDR
  V_3 = 0.0_KDR
  
  allocate ( T_Copy, source = T )
  
  call Timer % Stop ( )
  call Timer % ShowInterval (  )
  
  call Timer % Start ( )
  call Show ( 'Updating Device' )
  call F % UpdateDevice ( )
  call Timer % Stop ( )
  call Timer % ShowInterval (  )
  
  call Show ( 'ComputeFromTemperature' )
  call Show ( 'Input' )
  call Show ( [ N ( 1  ), T ( 1 ), E ( 1 ), Y_E ( 1 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromTemperature ( F % StorageForm, Dummy_G, Dummy_G )  
  call Timer % Stop ( )
  
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
  call Timer % ShowInterval (  )
  
  !-- Perturbed T by 20%
  T = 1.20 * T
  call F % UpdateDevice ( )
  
  call Show ( 'ComputeFromPrimitive' )
  call Show ( 'Input' )
  call Show ( [ N ( 1  ), T ( 1 ), E ( 1 ), Y_E ( 1 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromPrimitive ( F % StorageForm, Dummy_G, Dummy_G )
  call Timer % Stop ( )
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
  call Show ( sum ( abs ( T - T_Copy ) ) / sum ( abs ( T_Copy ) ), &
              'L1 Error Temperature' )
  call Timer % ShowInterval (  )
  
  !-- Perturbed T by 20%
  T = 1.15 * T
  call F % UpdateDevice ( )
  
  call Show ( 'ComputeFromConserved' )
  call Show ( 'Input' )
  call Show ( [ N ( 1  ), T ( 1 ), E ( 1 ), Y_E ( 1 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromConserved ( F % StorageForm, Dummy_G, Dummy_G )
  call Timer % Stop ( )
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
  call Show ( sum ( abs ( T - T_Copy ) ) / sum ( abs ( T_Copy ) ), &
              'L1 Error Temperature' )
  call Timer % ShowInterval (  )
  

!  call Timer % Start ( )
!  call Show ( 'Offload OpenMP NUC_EOS' )
!  !$OMP  OMP_TARGET_DIRECTIVE parallel do &
!  !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
!  !$OMP& private ( iV )
!  do iV = 1, nV
!    call NUC_EOS_FULL &
!         ( N ( iV ), T ( iV ), Y_E ( iV ), E ( iV ), &
!           P ( iV ), SB ( iV ), CS ( iV ), U_V ( iV ), &
!           U_V ( iV ), U_V ( iV ), X_He ( iV ), X_A ( iV ), &
!           X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
!           U_V ( iV ), U_V ( iV ), Mu_NP ( iV ), EOS_Apply_EOS_HN_T, &
!           Error ( iV ), EOS_RF_Accuracy, T_L_D, T_L_T, T_Y_E, T_EOS )
!  end do
!  call Timer % Stop ( )
!  call Timer % ShowInterval (  )
!
!  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
!  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
!  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
!  
!  call Timer % Start ( )
!  call Show ( 'CPU OpenMP NUC_EOS' )
!  !$OMP  parallel do &
!  !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
!  !$OMP& private ( iV )
!  do iV = 1, nV
!    call NUC_EOS_FULL &
!         ( N ( iV ), T ( iV ), Y_E ( iV ), E ( iV ), &
!           P ( iV ), SB ( iV ), CS ( iV ), U_V ( iV ), &
!           U_V ( iV ), U_V ( iV ), X_He ( iV ), X_A ( iV ), &
!           X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
!           U_V ( iV ), U_V ( iV ), Mu_NP ( iV ), EOS_Apply_EOS_HN_T, &
!           Error ( iV ), EOS_RF_Accuracy, T_L_D, T_L_T, T_Y_E, T_EOS )
!  end do
!  call Timer % Stop ( )
!  call Timer % ShowInterval (  )
!  
!  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
!  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
!  call Show ( [ X_P ( 1 ), X_N ( 1 ), X_He ( 1 ) ] , 'X_P, X_N, X_He' )
!  
!  
!  !-- Change temperature by 10%
!  T = 1.10 * T 
!  
!  call Show ( 'Input' )
!  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )

!  
!  do iV = 1, size ( N )
!    call NUC_EOS_FULL &
!         ( N ( iV ), T ( iV ), Y_E ( iV ), E ( iV ), &
!           P ( iV ), SB ( iV ), CS ( iV ), U_V ( iV ), &
!           U_V ( iV ), U_V ( iV ), X_He ( iV ), X_A ( iV ), &
!           X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
!           U_V ( iV ), U_V ( iV ), Mu_NP ( iV ), EOS_Apply_EOS_HN_E, &
!           Error ( iV ), EOS_RF_Accuracy, T_L_D, T_L_T, T_Y_E, T_EOS )
!  end do
!  call Show ( [ N ( 1 ), T ( 1 ), Y_E ( 1 ) ], 'N, T, Y_E' )
!  call Show ( [ P ( 1 ), E ( 1 ), CS ( 1 ) ], 'P, E, CS' )
  
  call CONSOLE % Mute ( )

  call Show ( 'Fluid_P_HN Variables', F % IGNORABILITY )
  call Show ( F % Variable, 'Variable', F % IGNORABILITY )
  
  end associate
  end associate

  deallocate ( F )
  deallocate ( PROGRAM_HEADER )

end program Fluid_P_HN__Form_Test
