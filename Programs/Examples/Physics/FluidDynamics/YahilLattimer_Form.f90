#include "Preprocessor"

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

  type, public, extends ( Universe_F_CC_Form ) :: YahilLattimerForm
    real ( KDR ) :: &
      AdiabaticIndex, &
      PolytropicConstant, &
      DensityDimensionless, &
      DensityInitial, &
      PressureInitial, &
      DensityFinal, &
      CollapseTime
    type ( InterpolationForm ), allocatable :: &
      Interpolation_D, &
      Interpolation_V
    type ( Fluid_P_I_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, private, pass :: &
      Initialize_H
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
    procedure, public, pass :: &
      ShowDiagnostics
  end type YahilLattimerForm

    private :: &
      InitializeUniverse, &
      InitializeDiagnostics, &
      SetInitial, &
      ResetInitial, &
      SetReference

      private :: &
        PrepareInterpolation, &
        SetFluid

        private :: &
          SetFluidKernel

    integer ( KDI ), private, parameter :: &
      iProfile_X = 1, &  !-- must match the Profile file columns
      iProfile_D = 2, &
      iProfile_V = 3

contains


  subroutine Initialize_H ( U, NameOption )

    class ( YahilLattimerForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a YahilLattimer'

    Name  =  'YahilLattimer'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )
    call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  subroutine ComputeError ( YL )

    class ( YahilLattimerForm ) , intent ( in ) :: &
      YL
    
    real ( KDR ) :: &
      L1_Rho, &
      L1_V, &
      L1_P
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( A  =>  YL % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  YL % Reference, &
        F_D  =>  YL % Difference )
    associate &
      ( FV_R  =>  F_R % Storage_GS % Value, &
        FV_D  =>  F_D % Storage_GS % Value )

    call CO % Initialize ( C % Communicator, [ 6 ], [ 6 ] )

    associate &
      ( D_Rho  =>  FV_D ( :, F_D % BARYON_DENSITY_C ), &
        R_Rho  =>  FV_R ( :, F_R % BARYON_DENSITY_C ), &
        D_V    =>  FV_D ( :, F_D % VELOCITY_U_1 ), &
        R_V    =>  FV_R ( :, F_R % VELOCITY_U_1 ), &
        D_P    =>  FV_D ( :, F_D % PRESSURE ), &
        R_P    =>  FV_R ( :, F_R % PRESSURE ), &
        Norm_D_Rho  =>  CO % Incoming % Value ( 1 ), &
        Norm_R_Rho  =>  CO % Incoming % Value ( 2 ), &
        Norm_D_V    =>  CO % Incoming % Value ( 3 ), &
        Norm_R_V    =>  CO % Incoming % Value ( 4 ), &
        Norm_D_P    =>  CO % Incoming % Value ( 5 ), &
        Norm_R_P    =>  CO % Incoming % Value ( 6 ) )

    CO % Outgoing % Value ( 1 ) &
      =  sum ( abs ( D_Rho ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 ) &
      =  sum ( abs ( R_Rho ), mask = C % ProperCell )
    CO % Outgoing % Value ( 3 ) &
      =  sum ( abs ( D_V ), mask = C % ProperCell )
    CO % Outgoing % Value ( 4 ) &
      =  sum ( abs ( R_V ), mask = C % ProperCell )
    CO % Outgoing % Value ( 5 ) &
      =  sum ( abs ( D_P ), mask = C % ProperCell )
    CO % Outgoing % Value ( 6 ) &
      =  sum ( abs ( R_P ), mask = C % ProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1_Rho  =  Norm_D_Rho / Norm_R_Rho
    L1_V    =  Norm_D_V   / Norm_R_V
    L1_P    =  Norm_D_P   / Norm_R_P

    call Show ( L1_Rho, '*** L1_Rho error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_V, '*** L1_V error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_P, '*** L1_P error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D_Rho, etc.
    end associate !-- FV_R, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine ComputeError


  impure elemental subroutine Finalize ( YL )
    
    type ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    if ( allocated ( YL % Difference ) ) &
      deallocate ( YL % Difference )
    if ( allocated ( YL % Reference ) ) &
      deallocate ( YL % Reference )
    if ( allocated ( YL % Interpolation_V ) ) &
      deallocate ( YL % Interpolation_V )
    if ( allocated ( YL % Interpolation_D ) ) &
      deallocate ( YL % Interpolation_D )

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( YahilLattimerForm ), intent ( in ) :: &
      U

    call U % Universe_F_CC_Form % ShowParameters ( )

    call Show ( U % AdiabaticIndex, 'AdiabaticIndex' )
    call Show ( U % DensityInitial, UNIT % MASS_DENSITY_CGS, &
                'MassDensityInitial' )
    call Show ( U % DensityInitial, UNIT % ENERGY_DENSITY_NUCLEAR, &
                'MassDensityInitial' )
    call Show ( U % DensityInitial  /  CONSTANT % ATOMIC_MASS_UNIT, &
                UNIT % NUMBER_DENSITY_NUCLEAR, 'BaryonDensityInitial' )
    call Show ( U % PressureInitial, UNIT % BARYE, 'PressureInitial' )
    call Show ( U % PressureInitial, UNIT % ENERGY_DENSITY_NUCLEAR, &
                'PressureInitial' )
    call Show ( U % DensityFinal, UNIT % MASS_DENSITY_CGS, &
                'MassDensityFinal' )
    call Show ( U % DensityFinal, UNIT % ENERGY_DENSITY_NUCLEAR, &
                'MassDensityFinal' )
    call Show ( U % DensityFinal  /  CONSTANT % ATOMIC_MASS_UNIT, &
                UNIT % NUMBER_DENSITY_NUCLEAR, 'BaryonDensityFinal' )
    call Show ( U % CollapseTime, UNIT % SECOND, 'CollapseTime' )

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

    class ( YahilLattimerForm ), intent ( in ) :: &
      U

    call U % Reference % Show ( )
    call U % Universe_F_CC_Form % ShowDiagnostics ( )

  end subroutine ShowDiagnostics


  subroutine InitializeUniverse ( YL, Name )

    class ( YahilLattimerForm ), intent ( inout ), target :: &
      YL
    character ( * ), intent ( in )  :: &
      Name
    
    call YL % Initialize &
           ( FluidType = 'IDEAL', &
             GravitationType = 'NEWTON_SG', &
             NameOption = Name, &
             nCellsPolarOption = 128 )

    YL % Integrator % SetInitial    =>  SetInitial
    YL % Integrator % ResetInitial  =>  ResetInitial
    YL % Integrator % SetReference  =>  SetReference
    YL % Integrator % System        =>  YL

  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( YL )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL

    allocate ( YL % Reference )
    allocate ( YL % Difference )
    associate &
      ( F_R  =>  YL % Reference, &
        F_D  =>  YL % Difference, &
        G    =>  YL % Integrator % Geometry_X, &
        S    =>  YL % Integrator % Checkpoint_X )

    call F_R % Initialize ( G, YL % Units_F, NameOption = 'Reference' )
    call F_D % Initialize ( G, YL % Units_F, NameOption = 'Difference' )
    call F_R % SetStream ( S )
    call F_D % SetStream ( S )

    end associate !-- FA_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    real ( KDR ) :: &
      GC

    select type ( YL  =>  I % System )
      class is ( YahilLattimerForm )
    select type ( I => YL % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    call PrepareInterpolation ( YL )

    GC  =  CONSTANT % GRAVITATIONAL

    associate &
      ( Gamma => YL % AdiabaticIndex, &
        Kappa => YL % PolytropicConstant, &
          T_C => YL % CollapseTime, &
          D_0 => YL % DensityDimensionless, &
        Rho_I => YL % DensityInitial, &
          P_I => YL % PressureInitial, &
        Rho_F => YL % DensityFinal, &
          T_F => I % T_Finish )

    Gamma  =  1.30_KDR
    Rho_I  =  7.0e9_KDR   *  UNIT % MASS_DENSITY_CGS
      P_I  =  6.0e27_KDR  *  UNIT % BARYE
    Rho_F  =  1.0e14_KDR  *  UNIT % MASS_DENSITY_CGS
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( Rho_I, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( P_I,   'PressureInitial' )
    call PROGRAM_HEADER % GetParameter ( Rho_F, 'DensityFinal' )

    Kappa  =  P_I  /  Rho_I ** Gamma
      T_C  =  sqrt ( D_0 / ( GC * Rho_I ) ) 
      T_F  =  T_C  -  sqrt ( D_0 / ( GC * Rho_F ) ) 

    call F % SetAdiabaticIndex &
           ( Gamma )
    call F % SetFiducialParameters &
           ( FiducialBaryonDensity = Rho_I, &
             FiducialPressure = P_I )

    call SetFluid ( YL, F )

    associate ( F_R  =>  YL % Reference )
      call F_R % SetAdiabaticIndex &
             ( Gamma )
      call F_R % SetFiducialParameters &
             ( FiducialBaryonDensity = Rho_I, &
               FiducialPressure = P_I )
    end associate !-- F_R

    if ( allocated ( YL % SA_Fluid ) ) then
      select type ( F_SA  =>  YL % SA_Fluid % FieldSet_SA )
      class is ( Fluid_P_I_Form )
        call F_SA % SetAdiabaticIndex &
               ( Gamma )
        call F_SA % SetFiducialParameters &
               ( FiducialBaryonDensity = Rho_I, &
                 FiducialPressure = P_I )
      end select !-- F_SA
    end if

    if ( allocated ( YL % AA_Fluid ) ) then
      select type ( F_AA  =>  YL % AA_Fluid % FieldSet_AA )
      class is ( Fluid_P_I_Form )
        call F_AA % SetAdiabaticIndex &
               ( Gamma )
        call F_AA % SetFiducialParameters &
               ( FiducialBaryonDensity = Rho_I, &
                 FiducialPressure = P_I )
      end select !-- F_AA
    end if

    end associate !-- Gamma, etc.

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- YL

  end subroutine SetInitial
  
  
  subroutine ResetInitial ( I, RestartFrom, T_Restart )
  
  
    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      RestartFrom
    type ( MeasuredValueForm ), intent ( out ) :: &
      T_Restart
    
    call SetInitial ( I )
    call I % ResetInitial_H ( RestartFrom, T_Restart )

  end subroutine ResetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( YL  =>  I % System )
      class is ( YahilLattimerForm )
    select type ( I  =>  YL % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    associate &
      ( F_R  =>  YL % Reference, &
        F_D  =>  YL % Difference )

    call SetFluid ( YL, F_R )
    call F_R % ComputeFromPrimitive ( F_R )

    call F_D % MultiplyAdd ( F, F_R, -1.0_KDR, UseDeviceOption = .false. )

    end associate !-- F_R, etc.
    end select !-- F
    end select !-- I
    end select !-- YL

  end subroutine SetReference


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

    allocate ( YL % Interpolation_D )
    allocate ( YL % Interpolation_V )
    associate &
      ( I_D => YL % Interpolation_D, &
        I_V => YL % Interpolation_V, &
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

    call I_D % Initialize &
           ( X, D, VerbosityOption = CONSOLE % INFO_3 )
    call I_V % Initialize &
           ( X, V, VerbosityOption = CONSOLE % INFO_3 )

    end associate !-- I_D, etc.

  end subroutine PrepareInterpolation


  subroutine SetFluid ( YL, F )

    class ( YahilLattimerForm ), intent ( inout ) :: &
      YL
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    call SetFluidKernel &
           (    N = FV ( :, F % BARYON_DENSITY_C ), &
                P = FV ( :, F % PRESSURE ), &
                E = FV ( :, F % ENERGY_DENSITY_C ), &
              V_1 = FV ( :, F % VELOCITY_U_1 ), &
              V_2 = FV ( :, F % VELOCITY_U_2 ), &
              V_3 = FV ( :, F % VELOCITY_U_3 ), &
              I_D = YL % Interpolation_D, &
              I_V = YL % Interpolation_V, &
             ProperCell = C % ProperCell, &
                R = GV ( :, G % CENTER_U_1 ), &
             Minus_t_YL = YL % CollapseTime - YL % Integrator % T, &
                  Gamma = YL % AdiabaticIndex, &
                  Kappa = YL % PolytropicConstant, &
                      G = CONSTANT % GRAVITATIONAL, &
                    amu = CONSTANT % ATOMIC_MASS_UNIT )
    
    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( N, P, E, V_1, V_2, V_3, I_D, I_V, ProperCell, R, &
                 Minus_t_YL, Gamma, Kappa, G, amu )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      N, &
      P, &
      E, &
      V_1, V_2, V_3
    type ( InterpolationForm ), intent ( in ) :: &
      I_D, &
      I_V
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R
    real ( KDR ), intent ( in ) :: &
      Minus_t_YL, &
      Gamma, &
      Kappa, &
      G, &
      amu

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    real ( KDR ) :: &
      X, &
      D, &
      V

    nV = size ( N )
      
    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( X, D, V )
    do iV  =  1,  nV

      if ( .not. ProperCell ( iV ) ) &
        cycle

      X  =  Kappa ** ( -0.5_KDR )  *  G ** ( ( Gamma - 1.0_KDR ) / 2.0_KDR )  &
            *  R ( iV )  *  Minus_t_YL ** ( Gamma - 2.0_KDR )
      
      call I_D % Evaluate ( X, D )
      call I_V % Evaluate ( X, V )

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
