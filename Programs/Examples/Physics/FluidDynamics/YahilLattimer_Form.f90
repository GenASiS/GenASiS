module YahilLattimer_Form

  !-- Yahil and Lattimer 1982, in Supernovae: A Survey of Current Research,
  !   ed. M. J. Rees and R. J. Stoneham, 53-70

  !-- Yahil 1983, ApJ 265, 1047-1055

  !-- (-t) = Minus_t_YL = CollapseTime - Time.
  !   Here (-t) is the time expression appearing in the above papers, in 
  !   terms of which infinite density ("catastrophe") occurs at t = 0. 
  !   "Time" ( >= 0 ) is the time variable in the code ("code time"),
  !   which begins at Time = 0. "CollapseTime" ( > 0 ) is the code time at 
  !   which catastrophe would be reached.

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: YahilLattimerForm
    real ( KDR ) :: &
      AdiabaticIndex, &
      PolytropicConstant, &
      DensityDimensionless, &
      DensityInitial, &
      PressureInitial, &
      CollapseTime
    type ( SplineInterpolationForm ), dimension ( : ), allocatable :: &
      SplineInterpolation
    type ( Fluid_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, private, pass :: &
      Initialize_YL
    generic, public :: &
      Initialize => Initialize_YL
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type YahilLattimerForm

    private :: &
      SetReference, &
      InitializeFluidCentralCore, &
      InitializeDiagnostics, &
      SetProblem

      private :: &
        PrepareInterpolation, &
        SetFluid

        private :: &
          SetFluidKernel

    integer ( KDI ), private, parameter :: &
      iProfile_X = 1, &  !-- must match the Profile file columns
      iProfile_D = 2, &
      iProfile_V = 3
    integer ( KDI ), private, parameter :: &
      iSpline_D = 1, & !-- spline interpolation
      iSpline_V = 2

contains


  subroutine Initialize_YL ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ), target :: &
      YL
    character ( * ), intent ( in ) :: &
      Name

    if ( YL % Type == '' ) &
      YL % Type = 'a YahilLattimer'

    call InitializeFluidCentralCore ( YL, Name )
    call InitializeDiagnostics ( YL )
    call SetProblem ( YL )

  end subroutine Initialize_YL


  subroutine ComputeError ( YL )

    class ( YahilLattimerForm ) , intent ( in ) :: &
      YL
    
    real ( KDR ) :: &
      L1_Rho, &
      L1_V, &
      L1_P
    class ( Fluid_P_I_Form ), pointer :: &
      F_D, &
      F_R
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( PS => YL % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    
    select type ( PSC => PS % Chart ) 
    class is ( Chart_SLD_Form )
    
    F_D => YL % Difference % Fluid_P_I ( )
    F_R => YL % Reference % Fluid_P_I ( )
    
    call CO % Initialize ( PS % Communicator, [ 6 ], [ 6 ] )

    associate &
      ( D_Rho => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
        R_Rho => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ), &
          D_V => F_D % Value ( :, F_D % VELOCITY_U ( 1 ) ), &
          R_V => F_R % Value ( :, F_R % VELOCITY_U ( 1 ) ), &
          D_P => F_D % Value ( :, F_D % PRESSURE ), &
          R_P => F_R % Value ( :, F_R % PRESSURE ), &
        Norm_D_Rho => CO % Incoming % Value ( 1 ), &
        Norm_R_Rho => CO % Incoming % Value ( 2 ), &
          Norm_D_V => CO % Incoming % Value ( 3 ), &
          Norm_R_V => CO % Incoming % Value ( 4 ), &
          Norm_D_P => CO % Incoming % Value ( 5 ), &
          Norm_R_P => CO % Incoming % Value ( 6 ) )

    CO % Outgoing % Value ( 1 ) &
      = sum ( abs ( D_Rho ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 2 ) &
      = sum ( abs ( R_Rho ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 3 ) &
      = sum ( abs ( D_V ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 4 ) &
      = sum ( abs ( R_V ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 5 ) &
      = sum ( abs ( D_P ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 6 ) &
      = sum ( abs ( R_P ), mask = PSC % IsProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1_Rho = Norm_D_Rho / Norm_R_Rho
      L1_V = Norm_D_V / Norm_R_V
      L1_P = Norm_D_P / Norm_R_P

    call Show ( L1_Rho, '*** L1_Rho error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_V, '*** L1_V error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_P, '*** L1_P error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D_Rho, etc.
    end select !-- PSC
    end select !-- PS
    nullify ( F_D, F_R )

  end subroutine ComputeError


  impure elemental subroutine Finalize ( YL )
    
    type ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    if ( allocated ( YL % Difference ) ) &
      deallocate ( YL % Difference )
    if ( allocated ( YL % Reference ) ) &
      deallocate ( YL % Reference )
    if ( allocated ( YL % SplineInterpolation ) ) &
      deallocate ( YL % SplineInterpolation )

  end subroutine Finalize


  subroutine SetReference ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference

    select type ( I )
    class is ( Integrator_C_PS_Form )

    select type ( YL => I % Universe )
    class is ( YahilLattimerForm )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( PSC => PS % Chart )
    class is ( Chart_SLD_Form )

    F_R => YL % Reference % Fluid_P_I ( )
    call SetFluid ( YL, F_R, G, PSC )

    F_D => YL % Difference % Fluid_P_I ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    end select !-- PSC
    end select !-- PS
    end select !-- FA
    end select !-- OS
    end select !-- I
    nullify ( G, F, F_R, F_D )

  end subroutine SetReference


  subroutine InitializeFluidCentralCore ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    character ( * ), intent ( in )  :: &
      Name
    
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = ( OffloadEnabled ( ) .and. GetNumberOfDevices ( ) >= 1 )
    call PROGRAM_HEADER % GetParameter ( UseDevice, 'UseDevice' )

    call YL % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', &
             Name = Name, FluidUseDeviceOption = UseDevice, &
             GeometryUseDeviceOption = UseDevice )
             
    YL % Integrator % SetReference => SetReference

  end subroutine InitializeFluidCentralCore


  subroutine InitializeDiagnostics ( YL )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    select type ( PS => YL % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( YL % Reference )
    allocate ( YL % Difference )
    call YL % Reference % Initialize &
           ( PS, 'IDEAL', YL % Units, NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             AllocateFeaturesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call YL % Difference % Initialize &
           ( PS, 'IDEAL', YL % Units, NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             AllocateFeaturesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine SetProblem ( YL )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    real ( KDR ) :: &
      GC, &
      DensityFinal, &
      CollapseTime
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => YL % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    call PrepareInterpolation ( YL )

    GC  =  CONSTANT % GRAVITATIONAL

    associate &
      ( Gamma => YL % AdiabaticIndex, &
        Kappa => YL % PolytropicConstant, &
          T_C => YL % CollapseTime, &
          D_0 => YL % DensityDimensionless, &
        Rho_I => YL % DensityInitial, &
          P_I => YL % PressureInitial, &
        Rho_F => DensityFinal, &
          T_F => I % FinishTime )

    Gamma  =  1.30_KDR
    Rho_I  =  7.0e9_KDR * UNIT % MASS_DENSITY_CGS
      P_I  =  6.0e27_KDR * UNIT % BARYE
    Rho_F  =  1.0e14_KDR * UNIT % MASS_DENSITY_CGS
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( Rho_I, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( P_I, 'PressureInitial' )
    call PROGRAM_HEADER % GetParameter ( Rho_F, 'DensityFinal' )

    Kappa  =  P_I  /  Rho_I ** Gamma
      T_C  =  sqrt ( D_0 / ( GC * Rho_I ) ) 
      T_F  =  T_C  -  sqrt ( D_0 / ( GC * Rho_F ) ) 

    call Show ( Gamma, 'AdiabaticIndex' )
    call Show ( Rho_I, UNIT % MASS_DENSITY_CGS, 'DensityInitial' )
    call Show ( P_I, UNIT % BARYE, 'PressureInitial' )
    call Show ( Rho_F, UNIT % MASS_DENSITY_CGS, 'DensityFinal' )
    call Show ( T_C, UNIT % SECOND, 'CollapseTime' )
    call Show ( T_F, UNIT % SECOND, 'Reset FinishTime' )

    select type ( PSC => PS % Chart )
    class is ( Chart_SLD_Form )

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( YL, F, G, PSC )

    end select !-- PSC
    end associate !-- Gamma, etc.
    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine PrepareInterpolation ( YL )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    real ( KDR ), dimension ( : ), allocatable :: &
      X, &
      D, &
      V
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Profile
    character ( LDF ) :: &
      Path, &
      Filename
    type ( TableStreamForm ) :: &
      TS

    call Show ( 'Preparing Interpolation' )

    Path = '../Parameters/'
    FileName = 'YahilLattimerCollapse_Gm_130.dat'
    
    call PROGRAM_HEADER % GetParameter ( Filename, 'Filename' )
    call Show ( Filename, 'Filename' )

    call TS % Initialize &
           ( Filename, PROGRAM_HEADER % Communicator % Rank, &
             PathOption = Path )
    call TS % Read ( Profile, oRowOption = 1 )

    YL % DensityDimensionless  =  Profile ( 1, iProfile_D )

    allocate ( YL % SplineInterpolation ( 2 ) )
    associate &
      ( SI => YL % SplineInterpolation, &
        nProfile => size ( Profile, dim = 1 ) )

    allocate ( X ( nProfile + 1 ) )
    allocate ( D ( nProfile + 1 ) ) 
    allocate ( V ( nProfile + 1 ) )

    X ( 2 : )  =  Profile ( :, iProfile_X )
    D ( 2 : )  =  Profile ( :, iProfile_D )
    V ( 2 : )  =  Profile ( :, iProfile_V )

    X ( 1 )  =  0.0_KDR
    D ( 1 )  =  D ( 2 )
    V ( 1 )  =  0.0_KDR

    call SI ( iSpline_D ) % Initialize &
           ( X, D, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iSpline_V ) % Initialize &
           ( X, V, VerbosityOption = CONSOLE % INFO_3 )

    end associate !-- SI, etc.

  end subroutine PrepareInterpolation


  subroutine SetFluid ( YL, F, G, PSC )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( Chart_SLD_Form ), intent ( inout ) :: &
      PSC
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    call F % SetAdiabaticIndex ( YL % AdiabaticIndex )
    call F % SetFiducialParameters &
           ( FiducialBaryonDensity = YL % DensityInitial, &
             FiducialPressure = YL % PressureInitial )

    call SetFluidKernel &
           (    N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                P = F % Value ( :, F % PRESSURE ), &
                E = F % Value ( :, F % INTERNAL_ENERGY ), &
              V_1 = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
              V_2 = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
              V_3 = F % Value ( :, F % VELOCITY_U ( 3 ) ), &
             SI_D = YL % SplineInterpolation ( iSpline_D ), &
             SI_V = YL % SplineInterpolation ( iSpline_V ), &
             IsProperCell = PSC % IsProperCell, &
                        R = G % Value ( :, G % CENTER_U ( 1 ) ), &
               Minus_t_YL = YL % CollapseTime - YL % Integrator % Time, &
                    Gamma = YL % AdiabaticIndex, &
                    Kappa = YL % PolytropicConstant, &
                        G = CONSTANT % GRAVITATIONAL, &
                      amu = CONSTANT % ATOMIC_MASS_UNIT )
    
  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( N, P, E, V_1, V_2, V_3, SI_D, SI_V, IsProperCell, R, &
                 Minus_t_YL, Gamma, Kappa, G, amu )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      P, &
      E, &
      V_1, V_2, V_3
    type ( SplineInterpolationForm ), intent ( in ) :: &
      SI_D, &
      SI_V
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R
    real ( KDR ) :: &
      Minus_t_YL, &
      Gamma, &
      Kappa, &
      G, &
      amu

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      X, &
      D, &
      V

    nValues = size ( N )
      
    !$OMP parallel do &
    !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV )
    do iV = 1, nValues

      if ( .not. IsProperCell ( iV ) ) &
        cycle

      X  =  Kappa ** ( -0.5_KDR )  *  G ** ( ( Gamma - 1.0_KDR ) / 2.0_KDR )  &
            *  R ( iV )  *  Minus_t_YL ** ( Gamma - 2.0_KDR )
 
      call SI_D % Evaluate ( X, D )
      call SI_V % Evaluate ( X, V )

      N ( iV )  =  D / ( amu  *  G  *  Minus_t_YL ** 2 )

      P ( iV )  =  Kappa  *  ( amu * N ( iV ) ) ** Gamma

      E ( iV )  =  P ( iV )  /  ( Gamma - 1.0_KDR )

      V_1 ( iV )  =  V  *  Kappa ** ( 0.5_KDR )  &
                     *  G ** ( ( 1.0_KDR - Gamma ) / 2.0_KDR )  &
                     *  Minus_t_YL ** ( 1.0_KDR - Gamma )

      V_2 ( iV )  =  0.0_KDR
      V_3 ( iV )  =  0.0_KDR

    end do
    !$OMP end parallel do

  end subroutine SetFluidKernel


end module YahilLattimer_Form
