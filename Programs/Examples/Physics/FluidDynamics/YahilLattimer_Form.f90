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

!     class ( Fluid_P_I_Form ), pointer :: &
!       F, &
!       F_R, &  !-- F_Reference
!       F_D     !-- F_Difference
!     integer ( KDI ) :: &
!       iV
!     real ( KDR ) :: &
!       R_Max

!     select type ( FCC )
!     class is ( FluidCentralCoreForm )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_P_I ( )

!     F_R => YahilLattimer % Reference % Fluid_P_I ( )
!     call SetFluid ( YahilLattimer, F_R, FCC % Time )

!     F_D => YahilLattimer % Difference % Fluid_P_I ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

!     F_D % Value = abs ( F_D % Value &
!                           / max ( abs ( F_R % Value ), &
!                                     sqrt ( tiny ( 0.0_KDR ) ) ) )

!     end select !-- FA
!     end select !-- FCC
!     nullify ( F, F_R, F_D )

  end subroutine SetReference


  subroutine InitializeFluidCentralCore ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    character ( * ), intent ( in )  :: &
      Name

    call YL % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', Name = Name )
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
             IgnorabilityOption = CONSOLE % INFO_2 )
    call YL % Difference % Initialize &
           ( PS, 'IDEAL', YL % Units, NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
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
    class is ( Chart_SL_Template )

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( YL, F, G, PSC % IsProperCell )

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


  subroutine SetFluid ( YL, F, G, IsProperCell )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell

    call F % SetAdiabaticIndex ( YL % AdiabaticIndex )
    call F % SetFiducialParameters &
           ( FiducialBaryonDensity = YL % DensityInitial, &
             FiducialPressure = YL % PressureInitial )

    call SetFluidKernel &
           (    N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                E = F % Value ( :, F % INTERNAL_ENERGY ), &
              V_1 = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
              V_2 = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
              V_3 = F % Value ( :, F % VELOCITY_U ( 3 ) ), &
             SI_D = YL % SplineInterpolation ( iSpline_D ), &
             SI_V = YL % SplineInterpolation ( iSpline_V ), &
             IsProperCell = IsProperCell, &
                        R = G % Value ( :, G % CENTER_U ( 1 ) ), &
               Minus_t_YL = YL % CollapseTime - YL % Integrator % Time, &
                    Gamma = YL % AdiabaticIndex, &
                    Kappa = YL % PolytropicConstant, &
                        G = CONSTANT % GRAVITATIONAL, &
                      amu = CONSTANT % ATOMIC_MASS_UNIT )

    call F % ComputeFromPrimitive ( G )

  end subroutine SetFluid


!   subroutine SetFluid ( YL, F, Time )

!     class ( YahilLattimerForm ), intent ( inout ) :: &
!       YL
!     class ( Fluid_P_I_Form ), intent ( inout ) :: &
!       F
!     real ( KDR ), intent ( in ) :: &
!       Time

!     class ( GeometryFlatForm ), pointer :: &
!       G
!     integer ( KDI ) :: &
!       iV
!     real ( KDR ) :: &
!       T, &
!       R_TM
!     real ( KDR ), dimension ( : ), allocatable :: &
!       R, &
!       Rho, &
!       V
!     type ( SplineInterpolationForm ), dimension ( 2 ) :: &
!       SI

!     !-- Interpolate self similar solution

!     T = YL % t_collapse - Time

!     call PrepareInterpolation &
!            ( SI, YL % AnalyticProfile, YL % Kappa, &
!              YL % AdiabaticIndex, T )

!     select type ( PS => YL % Integrator % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )

!     select type ( Grid => PS % Chart )
!     class is ( Chart_SL_Template )

!     select type ( C => PS % Chart )
!     class is ( Chart_SLD_Form )

!     associate &
!       ( R    => G % Value ( :, G % CENTER_U ( 1 ) ), &
!         Rho  => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
!         V_1  => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
!         V_2  => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
!         V_3  => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
!         E    => F % Value ( :, F % INTERNAL_ENERGY ), &
!         P    => F % Value ( :, F % PRESSURE ) )

!     V_2 = 0.0_KDR
!     V_3 = 0.0_KDR

!       do iV = 1, size ( R )
!         if ( R ( iV ) < 0.0_KDR ) &
!              cycle
!         call SI ( iRHO_SI ) &
!                 % Evaluate ( R ( iV ), Rho ( iV ) )
!         call SI ( iV_SI ) &
!                 % Evaluate ( R ( iV ), V_1 ( iV ) )
!         P ( iV ) = YL % Kappa * Rho ( iV ) ** ( YL % AdiabaticIndex )
!         E ( iV ) = P ( iV ) / ( YL % AdiabaticIndex - 1.0_KDR )
!       end do
      
!       Rho = Rho / CONSTANT % ATOMIC_MASS_UNIT

!       call F % ComputeFromPrimitive ( G )

!     end associate !-- R, etc.
!     end select    !-- C
!     end select    !-- Grid
!     end select    !-- PS

!     nullify ( G )

!   end subroutine SetFluid


  subroutine SetFluidKernel &
               ( N, E, V_1, V_2, V_3, SI_D, SI_V, IsProperCell, R, &
                 Minus_t_YL, Gamma, Kappa, G, amu )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
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
      
    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      if ( .not. IsProperCell ( iV ) ) &
        cycle

      X  =  Kappa ** ( -0.5_KDR )  *  G ** ( ( Gamma - 1.0_KDR ) / 2.0_KDR )  &
            *  R ( iV )  *  Minus_t_YL ** ( Gamma - 2.0_KDR )
 
      call SI_D % Evaluate ( X, D )
      call SI_V % Evaluate ( X, V )

      N ( iV )  =  D / ( amu  *  G  *  Minus_t_YL ** 2 )

      E ( iV )  =  ( Kappa  *  ( amu * N ( iV ) ) ** Gamma )  &
                   /  ( Gamma - 1.0_KDR )

      V_1 ( iV )  =  V  *  Kappa ** ( 0.5_KDR )  &
                     *  G ** ( ( 1.0_KDR - Gamma ) / 2.0_KDR )  &
                     *  Minus_t_YL ** ( 1.0_KDR - Gamma )

      V_2 ( iV )  =  0.0_KDR
      V_3 ( iV )  =  0.0_KDR

    end do
    !$OMP end parallel do

  end subroutine SetFluidKernel


  subroutine ComputeError ( YL )

    class ( YahilLattimerForm ) , intent ( in ) :: &
      YL
    
!     real ( KDR ) :: &
!       L1_Rho, &
!       L1_V, &
!       L1_P
!     class ( Fluid_P_I_Form ), pointer :: &
!       F, &
!       F_D, &
!       F_R
!     type ( CollectiveOperation_R_Form ) :: &
!       CO

!     select type ( FCC => YL % Integrator )
!     type is ( FluidCentralCoreForm )
!     select type ( PS => FCC % PositionSpace )
!     class is ( Atlas_SC_Form ) 
!     select type ( C => PS % Chart )
!     class is ( Chart_SL_Template )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_P_I ( )

!     F_R => YL % Reference % Fluid_P_I ( )
!     call SetFluid ( YahilLattimer, F_R, FCC % Time )

!     F_D => YL % Difference % Fluid_P_I ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )
       
    
!     associate &
!       ( Difference_Rho &
!           => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
!         Reference_Rho &
!           => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ), &
!         Difference_V   &
!           => F_D % Value ( :, F_D % VELOCITY_U ( 1 ) ), &
!         Reference_V   &
!           => F_R % Value ( :, F_R % VELOCITY_U ( 1 ) ), &
!         Difference_P   &
!           => F_D % Value ( :, F_D % PRESSURE ), &
!         Reference_P   &
!           => F_R % Value ( :, F_R % PRESSURE ) )

!     call CO % Initialize ( PS % Communicator, [ 6 ], [ 6 ] )

!     CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference_Rho ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 2 ) = sum ( abs ( Reference_Rho ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 3 ) = sum ( abs ( Difference_V ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 4 ) = sum ( abs ( Reference_V ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 5 ) = sum ( abs ( Difference_P ), &
!                                         mask = C % IsProperCell )
!     CO % Outgoing % Value ( 6 ) = sum ( abs ( Reference_P ), &
!                                         mask = C % IsProperCell )

!     call CO % Reduce ( REDUCTION % SUM )

!     end associate !-- Difference_Rho, etc.
    
!     associate &
!       ( DifferenceSum_Rho  => CO % Incoming % Value ( 1 ), &
!         ReferenceSum_Rho   => CO % Incoming % Value ( 2 ), &
!         DifferenceSum_V    => CO % Incoming % Value ( 3 ), &
!         ReferenceSum_V     => CO % Incoming % Value ( 4 ), &
!         DifferenceSum_P    => CO % Incoming % Value ( 5 ), &
!         ReferenceSum_P     => CO % Incoming % Value ( 6 ) )

!     L1_Rho = DifferenceSum_Rho / ReferenceSum_Rho
!     L1_V   = DifferenceSum_V   / ReferenceSum_V
!     L1_P   = DifferenceSum_P   / ReferenceSum_P
!     end associate

!     call Show ( L1_Rho, '*** L1_Rho error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )
!     call Show ( L1_V, '*** L1_V error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )
!     call Show ( L1_P, '*** L1_P error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )

!     end select !-- FA
!     end select !-- C
!     end select !-- PS
!     end select !--FCC

!     nullify ( F, F_R, F_D )

  end subroutine ComputeError


end module YahilLattimer_Form
