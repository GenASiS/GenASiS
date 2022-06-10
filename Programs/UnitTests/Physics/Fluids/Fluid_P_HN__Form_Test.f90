program Fluid_P_HN__Form_Test

  !-- Fluid_Perfect_HeavyNucleus__Form_Test

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none

  integer ( KDI ), parameter :: &
    iRADIUS_TS            = 2, &  !-- must match the profile file columns
    iRADIAL_VELOCITY_TS   = 3, &
    iDENSITY_TS           = 4, &
    iTEMPERATURE_TS       = 5, &
    iSPECIFIC_ENERGY_TS   = 10, &
    iELECTRON_FRACTION_TS = 11

  integer ( KDI ) :: &
!    iD, &
    iV, &
    iP
  real ( KDR ) :: &
    RadiusMax, &
    RadiusCore, &
    RadialRatio
  real ( KDR ), dimension ( : ), allocatable :: &
    T_Copy
  real ( KDR ), dimension ( :, : ), allocatable :: &
    Profile
  character ( LDF ) :: &
    Filename, &
    Path
  type ( TimerForm ) :: &
    Timer
  type ( TableStreamForm ) :: &
    TS  
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Gravitation_G_Form ), allocatable :: &
    G
  type ( Units_F_Form ), dimension ( : ), allocatable :: &
    U
  type ( Fluid_P_HN_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Fluid_P_HN__Form_Test', DimensionalityOption = '1D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( U ( 1 ) )
  call U ( 1 ) % Initialize ( TypeOption = 'ASTROPHYSICS' )

  RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
  RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
  RadialRatio  =  2.4_KDR

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = RadiusMax, &
           RadiusCore = RadiusCore, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace', &
           CoordinateUnitOption = U ( 1 ) % Coordinate_PS, &
           RadialRatioOption = RadialRatio )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( F )
  call F % Initialize ( G, U )
  call F % SetStream ( S )

  call A       % Show ( )
  call U ( 1 ) % Show ( )
  call F       % Show ( )
  call S       % Show ( )

  call Timer % Initialize ( 'Timer', Level = 1 )
  
  Path = '../../../../Applications/WoosleyHeger_07/Parameters/'
  Filename = 'WH07_S12_08.d.stripped'
  call TS % Initialize &
         ( Filename, PROGRAM_HEADER % Communicator % Rank, &
           PathOption = Path )
!  call TS % Read ( Profile, nRowsOption = A % Chart_GS % nCells ( 1 ), &
!                   oRowOption = 2 )
  call TS % Read ( Profile, nRowsOption = 25, &
                   oRowOption = 2 )
call Show (  A % Chart_GS % nCells ( 1 ), '>>> nCells ( 1 ) ' )
call Show ( shape ( Profile ), '>>> shape ( Profile )' )
  
  associate ( FV => F % Storage_GS % Value )
  associate &
    ( N    => FV ( :, F % BARYON_DENSITY_C ), &
      E    => FV ( :, F % ENERGY_DENSITY_C ), & 
      P    => FV ( :, F % PRESSURE ), &   
      T    => FV ( :, F % TEMPERATURE ), &
      CS   => FV ( :, F % SOUND_SPEED ), &
      Y_E  => FV ( :, F % ELECTRON_FRACTION ), &
      X_P  => FV ( :, F % MASS_FRACTION_PROTON ), & 
      X_N  => FV ( :, F % MASS_FRACTION_NEUTRON ), &
      X_He => FV ( :, F % MASS_FRACTION_ALPHA ), &
      X_A  => FV ( :, F % MASS_FRACTION_HEAVY ) &
    )

  call Timer % Start ( )
  call Show ( 'Setting Initial Values' )

  iP = 1  
  do iV = 3, size ( N ) - 2
    N ( iV )    =  Profile ( iP, iDENSITY_TS ) / CONSTANT % ATOMIC_MASS_UNIT
    E ( iV )    =  Profile ( iP, iSPECIFIC_ENERGY_TS ) &
                   *  Profile ( iP, iDENSITY_TS )
    T ( iV )    =  Profile ( iP, iTEMPERATURE_TS )
    Y_E ( iV )  =  Profile ( iP, iELECTRON_FRACTION_TS )
    
    iP = iP + 1
    if ( iP > size ( Profile, dim = 1 ) ) &
      iP = 1
  end do
  
  !-- Set unit
  N  =  N  *  UNIT % MASS_DENSITY_CGS
  T  =  T  *  UNIT % KELVIN
  E  =  E  *  ( UNIT % ERG / UNIT % GRAM ) * UNIT % MASS_DENSITY_CGS
  
!   V_1 = 0.0_KDR
!   V_2 = 0.0_KDR
!   V_3 = 0.0_KDR
  
  allocate ( T_Copy, source = T )
  
  call Timer % Stop ( )
  call Timer % ShowInterval ( CONSOLE % INFO_1 )
  
  call Timer % Start ( )
  call Show ( 'Updating Device' )
  call F % UpdateDevice ( )
  call Timer % Stop ( )
  call Timer % ShowInterval ( CONSOLE % INFO_1 )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  call Show ( 'ComputeFromTemperature' )
  call Show ( 'Input' )
  call Show ( [ N ( 3 ), T ( 3 ), E ( 3 ), Y_E ( 3 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromTemperature ( )  
  call Timer % Stop ( )
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 3 ), T ( 3 ), Y_E ( 3 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 3 ), E ( 3 ), CS ( 3 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 3 ), X_N ( 3 ), X_He ( 3 ), X_A ( 3 ) ] , &
              'X_P, X_N, X_He, X_A' )
  call Timer % ShowInterval ( CONSOLE % INFO_1 )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  !-- Perturbed T by 20%
  T = 1.20 * T
  call F % UpdateDevice ( )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  call Show ( 'ComputeFromPrimitive' )
  call Show ( 'Input' )
  call Show ( [ N ( 3 ), T ( 3 ), E ( 3 ), Y_E ( 3 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromPrimitive ( F )
  call Timer % Stop ( )
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 3 ), T ( 3 ), Y_E ( 3 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 3 ), E ( 3 ), CS ( 3 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 3 ), X_N ( 3 ), X_He ( 3 ), X_A ( 3 ) ] , &
              'X_P, X_N, X_He, X_A' )
  call Show ( sum ( abs ( T - T_Copy ) ) / sum ( abs ( T_Copy ) ), &
              'L1 Error Temperature' )
  call Timer % ShowInterval ( CONSOLE % INFO_1 )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  !-- Perturbed T by 15%
  T = 1.15 * T
  call F % UpdateDevice ( )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

  call Show ( 'ComputeFromBalanced' )
  call Show ( 'Input' )
  call Show ( [ N ( 3 ), T ( 3 ), E ( 3 ), Y_E ( 3 ) ], 'Input N, T, E, Y_E' )
  call Timer % Start ( )
  call F % ComputeFromBalanced ( )
  call Timer % Stop ( )
  call F % UpdateHost ( )
  call Show ( 'Output' )
  call Show ( [ N ( 3 ), T ( 3 ), Y_E ( 3 ) ], 'N, T, Y_E' )
  call Show ( [ P ( 3 ), E ( 3 ), CS ( 3 ) ], 'P, E, CS' )
  call Show ( [ X_P ( 3 ), X_N ( 3 ), X_He ( 3 ), X_A ( 3 ) ] , &
              'X_P, X_N, X_He, X_A' )
  call Show ( sum ( abs ( T - T_Copy ) ) / sum ( abs ( T_Copy ) ), &
              'L1 Error Temperature' )
  call Timer % ShowInterval ( CONSOLE % INFO_1 )
  
  call GIS % Open ( GIS % ACCESS_CREATE )
  call S % Write ( )
  call GIS % Close ( )

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
  
!   call CONSOLE % Mute ( )

  end associate !-- FV
  end associate !-- N, etc. 

  deallocate ( F )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( U )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Fluid_P_HN__Form_Test
