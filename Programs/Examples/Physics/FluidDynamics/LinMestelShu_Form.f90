module LinMestelShu_Form

  !-- From Lin, Mestel, and Shu 1965,
  !   "The Gravitational Collapse of a Uniform Spheroid"

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CC_Form ) :: LinMestelShuForm
    real ( KDR ) :: &
      Mass, &
      DensityInitial, &
      Eccentricity, &
      SemiMajor, &
      SemiMinor, &
      DensityFactor, &
      Density_OS, &        !-- OppenheimerSnyder, initial
      Radius_OS, &         !-- OppenheimerSnyder, initial
      RadiusFactor_OS, &   !-- OppenheimerSnyder
      TimeScale_OS, &      !-- OppenheimerSnyder
      TimeFinal_OS, &      !-- OppenheimerSnyder
      AtmosphereParameter
    type ( DifferentialEquationForm ), allocatable :: &
      DifferentialEquation
    type ( Fluid_D_Form ), allocatable :: &
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
  end type LinMestelShuForm

    private :: &
      InitializeUniverse, &
      InitializeDiagnostics, &
      SetInitial, &
      SetReference

      private :: &
        ComputeSlope_LMS, &
        SetFinishTime, &
        SetFluid

        private :: &
          Zero_T, &
          SetFluidKernel


contains


  subroutine Initialize_H ( U, NameOption )

    class ( LinMestelShuForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional  :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a LinMestelShu'

    Name  =  'LinMestelShu'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )
    call InitializeDiagnostics ( U )

  end subroutine Initialize_H


  subroutine ComputeError ( LMS )

    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS

    real ( KDR ) :: &
      L1
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( A  =>  LMS % Reference % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C    =>  A % Chart_GS, &
        F_R  =>  LMS % Reference, &
        F_D  =>  LMS % Difference )
    associate &
      ( FV_R  =>  F_R % Storage_GS % Value, &
        FV_D  =>  F_D % Storage_GS % Value )

    call CO % Initialize ( C % Communicator, [ 2 ], [ 2 ] )

    associate &
      ( D  =>  FV_D ( :, F_D % BARYON_DENSITY_C ), &
        R  =>  FV_R ( :, F_R % BARYON_DENSITY_C ), &
        Norm_D  =>  CO % Incoming % Value ( 1 ), &
        Norm_R  =>  CO % Incoming % Value ( 2 ) )

    CO % Outgoing % Value ( 1 ) &
      =  sum ( abs ( D ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 ) &
      =  sum ( abs ( R ), mask = C % ProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1  =  Norm_D / Norm_R
    call Show ( L1, '*** L1 error', &
                nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D, etc.
    end associate !-- FV_R, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine ComputeError


  impure elemental subroutine Finalize ( LMS )

    type ( LinMestelShuForm ), intent ( inout ) :: &
      LMS 
    
    if ( allocated ( LMS % Difference ) ) &
      deallocate ( LMS % Difference )
    if ( allocated ( LMS % Reference ) ) &
      deallocate ( LMS % Reference )
    if ( allocated ( LMS % DifferentialEquation ) ) &
      deallocate ( LMS % DifferentialEquation )

  end subroutine Finalize

  
  subroutine ShowParameters ( U )

    class ( LinMestelShuForm ), intent ( in ) :: &
      U

    real ( KDR ) :: &
      Pi

    call U % Universe_F_CC_Form % ShowParameters ( )

    Pi  =  CONSTANT % PI

    call Show ( U % Mass, 'Mass' )
    call Show ( U % DensityInitial, 'DensityInitial' )
    call Show ( U % Eccentricity, 'Eccentricity' )
    call Show ( U % SemiMajor, 'SemiMajor' )
    call Show ( U % SemiMinor, 'SemiMinor' )
    call Show ( U % DensityFactor, 'DensityFactor' )
    call Show ( U % Density_OS, 'Density_OS' )
    call Show ( U % Radius_OS, 'Radius_OS' )
    call Show ( U % RadiusFactor_OS, 'RadiusFactor_OS' )
    call Show ( Pi / 2  *  U % TimeScale_OS, 'TimeSingularity_OS' )
    call Show ( U % TimeFinal_OS, 'TimeFinal_OS' )
    call Show ( U % AtmosphereParameter, 'AtmosphereParameter' )

  end subroutine ShowParameters


  subroutine ShowDiagnostics ( U )

    class ( LinMestelShuForm ), intent ( in ) :: &
      U

    call U % Reference % Show ( )
    call U % Universe_F_CC_Form % ShowDiagnostics ( )

  end subroutine ShowDiagnostics


  subroutine InitializeUniverse ( LMS, Name )

    class ( LinMestelShuForm ), intent ( inout ), target :: &
      LMS
    character ( * ), intent ( in )  :: &
      Name

    call LMS % Initialize &
           ( FluidType = 'DUST', &
             GravitationType = 'NEWTON_SG', &
             NameOption = Name, &
             DimensionlessOption = .true., &
             GravityFactorOption = 0.01_KDR, &
             nCellsPolarOption = 128 )

    !-- Modify from default OUTFLOW to INFLOW outer radial boundary condition
    select type ( I  =>  LMS % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'INFLOW    ' ], iC = 1, iD = 1 )
    end associate !-- F
    end select !-- I

    LMS % Integrator % SetInitial    =>  SetInitial
    LMS % Integrator % SetReference  =>  SetReference
    LMS % Integrator % System        =>  LMS

  end subroutine InitializeUniverse


  subroutine InitializeDiagnostics ( LMS )

    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS

    allocate &
      ( LMS % Reference, &
        LMS % Difference )
    associate &
      ( F_R  =>  LMS % Reference, &
        F_D  =>  LMS % Difference, &
        G    =>  LMS % Integrator % Geometry_X, &
        S    =>  LMS % Integrator % Checkpoint_X )

    call F_R % Initialize ( G, LMS % Units_F, NameOption = 'Reference' )
    call F_D % Initialize ( G, LMS % Units_F, NameOption = 'Difference' )
    call F_R % SetStream ( S )
    call F_D % SetStream ( S )

    allocate ( LMS % DifferentialEquation )
    associate ( DE  =>  LMS % DifferentialEquation )
    call DE % Initialize ( LMS, nEquations = 4 )
    DE % ComputeSlope  =>  ComputeSlope_LMS
    end associate !-- R

    end associate !-- FA_R, etc.

  end subroutine InitializeDiagnostics


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    real ( KDR ) :: &
      Pi, &
      Eta_OS

    select type ( LMS  =>  I % System )
      class is ( LinMestelShuForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    associate &
      (   M      =>  LMS % Mass, &
          D_0    =>  LMS % DensityInitial, &
          E_0    =>  LMS % Eccentricity, &
          R_0    =>  LMS % SemiMajor, &
          Z_0    =>  LMS % SemiMinor, &
         DF      =>  LMS % DensityFactor, &
          D_OS   =>  LMS % Density_OS, &
          R_OS   =>  LMS % Radius_OS, &
         RF_OS   =>  LMS % RadiusFactor_OS, &
        Tau_OS   =>  LMS % TimeScale_OS, &
         TF_OS   =>  LMS % TimeFinal_OS, &
         AP      =>  LMS % AtmosphereParameter, &
          R_Max  =>  C % MaxCoordinate ( 1 ) )

      Pi    =  CONSTANT % PI
      M     =  1.0_KDR
      D_0   =  1.0e-3_KDR
      E_0   =  0.6_KDR
     DF     =  1.0e1_KDR
     AP     =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter (   M, 'Mass' )
    call PROGRAM_HEADER % GetParameter ( D_0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( E_0, 'Eccentricity' )
    call PROGRAM_HEADER % GetParameter (  DF, 'DensityFactor' )
    call PROGRAM_HEADER % GetParameter (  AP, 'AtmosphereParameter' )

      D_OS  =  D_0
      R_OS  =  ( 3.0 * M / ( 4.0 * Pi * D_0 ) ) ** ( 1.0_KDR / 3.0_KDR )
     RF_OS  =  DF ** ( - 1.0_KDR / 3.0_KDR ) 
    Tau_OS  =  sqrt ( 3.0 / ( 8.0 * Pi * D_0 ) )
    Eta_OS  =  acos ( 2.0 * RF_OS  -  1.0 )
     TF_OS  =  0.5 * Tau_OS * ( Eta_OS  +  sin ( Eta_OS ) )

    call SetFinishTime ( LMS, TF_OS )

    if ( R_0  >  R_Max  ) then
      call Show ( 'SemiMajor axis too large', CONSOLE % ERROR )
      call Show ( R_0, 'SemiMajor', CONSOLE % ERROR )
      call Show ( R_Max, 'RadiusMax', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    !-- Demand that the spheroid have the same volume as a sphere
    !   (OS = OppenheimerSnyder case) with the same mass and density
    R_0  =  R_OS  *  ( 1.0_KDR  -  E_0 ** 2 ) ** ( - 1.0_KDR / 6.0_KDR )
    Z_0  =  R_0  *  sqrt ( 1.0_KDR  -  E_0 ** 2 )

    call SetFluid ( LMS, F )
    call F % SetBaryonDensityMin ( )

    end associate !-- e0, etc.

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- LMS

  end subroutine SetInitial


  subroutine SetReference ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( LMS  =>  I % System )
      class is ( LinMestelShuForm )
    select type ( I  =>  LMS % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_D_Form )
    associate &
      ( F_R  =>  LMS % Reference, &
        F_D  =>  LMS % Difference )

    call SetFluid ( LMS, F_R )

    call F_D % MultiplyAdd ( F, F_R, -1.0_KDR, UseDeviceOption = .false. )

    end associate !-- F_R, etc.
    end select !-- F
    end select !-- I
    end select !-- LMS

  end subroutine SetReference


  ! subroutine ComputeSlope_LMS ( LMS, X, Y, dYdX )

  !   class ( * ), intent ( in ) :: &
  !     LMS
  !   real ( KDR ), intent ( in ) :: &
  !     X
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     Y
  !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  !     dYdX

  !   real ( KDR ) :: &
  !     SqrtTiny, &
  !     TwoPi, &
  !     E, &
  !     A, &
  !     C
    
  !   SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
  !   TwoPi     =  2.0_KDR  *  CONSTANT % PI

  !   select type ( LMS )
  !     class is ( LinMestelShuForm )
  !   associate &
  !     ( E_0  =>  LMS % Eccentricity, &
  !       D_0  =>  LMS % DensityInitial )

  !   E  =  max ( sqrt ( max ( 1.0_KDR &
  !                            -  ( Y ( 3 )  /  Y ( 1 ) ) ** 2  &
  !                               *  ( 1.0_KDR  -  E_0 ** 2 ),  &
  !                            SqrtTiny ) ), &
  !               SqrtTiny )

  !   A  =  TwoPi  *  sqrt ( ( 1.0_KDR  -  max ( E ** 2, SqrtTiny ) ) )  &
  !         /  max ( E ** 3, SqrtTiny )  &
  !         *  ( asin ( E )  &
  !              -  E  *  sqrt ( ( 1.0_KDR - max ( E ** 2, SqrtTiny ) ) ) )

  !   C  =  2.0_KDR  *  TwoPi  /  max ( E ** 2, SqrtTiny )  &
  !         *  ( 1.0_KDR  -  sqrt ( 1.0_KDR  -  max ( E ** 2, SqrtTiny ) ) &
  !                          *  asin ( E )  /  E )

  !   dYdX ( 1 ) = Y ( 2 )
  !   dYdX ( 3 ) = Y ( 4 )

  !   if ( ( Y ( 1 )  *  Y ( 3 ) )  >  0.0_KDR ) then
  !     dYdX ( 2 )  =  - D_0  *  A  &
  !                      /  max ( ( Y ( 1 )  *  Y ( 3 ) ), SqrtTiny )
  !   else
  !     dYdX ( 2 )  =  - D_0  *  A  &
  !                      /  min ( ( Y ( 1 ) * Y ( 3 ) ), - SqrtTiny )
  !   end if

  !   dYdX ( 4 )  =  - D_0  *  C &
  !                    /  max ( ( Y ( 1 ) ** 2 ), SqrtTiny )

  !   end associate !-- E_0 
  !   end select !-- LMS

  ! end subroutine ComputeSlope_LMS


  subroutine ComputeSlope_LMS ( LMS, X, Y, dYdX )

    class ( * ), intent ( in ) :: &
      LMS
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      dYdX

    real ( KDR ) :: &
      SqrtTiny, &
      Pi, &
      E, &
      A, &
      C
    
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    Pi  =  CONSTANT % PI

    select type ( LMS )
      class is ( LinMestelShuForm )
    associate &
      ( E_0  =>  LMS % Eccentricity, &
        D_0  =>  LMS % DensityInitial )

    E  =  sqrt ( 1.0_KDR &
                 -  ( Y ( 3 )  /  Y ( 1 ) ) ** 2  &
                      *  ( 1.0_KDR  -  E_0**2 ) )

    if ( E > 1.0e-2_KDR ) then

      A  =  2 * Pi  *  sqrt ( 1.0_KDR  -  E**2 )  &
            /  E**3  &
            *  ( asin ( E )  &
                 -  E  *  sqrt ( 1.0_KDR - E**2 ) )

      C  =  4 * Pi  /  E**2  &
            *  ( 1.0_KDR  -  sqrt ( 1.0_KDR  -  E**2 ) &
                             *  asin ( E )  /  E )
      
    else

      A  =  4 * Pi / 3  -  14 * Pi * E**2 / 15  &
            -  13 * Pi * E**4 / 70  -  19 * Pi * E**6 / 252  &
            -  125 * Pi * E**8 / 3168

      C  =  4 * Pi / 3  +  8 * Pi * E**6 / 15  &
            +  32 * Pi * E**4 / 105  +  64 * Pi * E**6 / 315 &
            +  512 * Pi * E**8 / 3465

    end if

    dYdX ( 1 ) = Y ( 2 )
    dYdX ( 3 ) = Y ( 4 )

    dYdX ( 2 )  =  - D_0  *  A  /  ( Y ( 1 )  *  Y ( 3 ) )
    dYdX ( 4 )  =  - D_0  *  C  /  Y ( 1 ) ** 2

    end associate !-- E_0 
    end select !-- LMS

  end subroutine ComputeSlope_LMS


  subroutine SetFinishTime ( LMS, T_Finish_OS )

    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS
    real ( KDR ), intent ( in ) :: &
      T_Finish_OS

    type ( RootForm ) :: &
      R

    call R % Initialize ( LMS )
    R % Zero  =>  Zero_T

    call Show ( 'Solving for T_Finish' )
    call R % Solve ( [ 0.9_KDR  *  T_Finish_OS,  &
                       1.05_KDR  *  T_Finish_OS ], &
                     LMS % Integrator % T_Finish )
    if ( R % Success ) then
      call Show ( 'Solve for T_Finish succeeded' )
      call Show ( R % nIterations, 'nIterations' )
      call Show ( R % Accuracy, 'Accuracy' )
      call Show ( R % RequestedAccuracy, 'RequestedAccuracy' )
    else
      call Show ( 'Solve for T_Finish failed', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

  end subroutine SetFinishTime


  subroutine SetFluid ( LMS, F )

    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F

    real ( KDR ) :: &
      X_Start, &
      X_Finish, &
      H_Start, &
      A_1, &
      A_3, &
      D

    associate &
      ( DE  =>  LMS % DifferentialEquation, &
        T   =>  LMS % Integrator % T )
    associate &
      ( Y  =>  DE % Solution )

    if ( T  >  0.0_KDR ) then

      X_Start   =  0.0_KDR
      X_Finish  =  T
      H_Start   =  LMS % TimeScale_OS  *  1.0e-3

      Y ( 1 )  =  1.0_KDR
      Y ( 2 )  =  0.0_KDR
      Y ( 3 )  =  1.0_KDR
      Y ( 4 )  =  0.0_KDR

      call DE % Integrate ( X_Start, X_Finish, H_Start )

      A_1  =  Y ( 1 )  *  LMS % SemiMajor
      A_3  =  Y ( 3 )  *  LMS % SemiMinor
      D    =  LMS % DensityInitial  /  ( Y ( 1 ) ** 2  *  Y ( 3 ) )

    else

      A_1  =  LMS % SemiMajor
      A_3  =  LMS % SemiMinor
      D    =  LMS % DensityInitial

    end if

    end associate !-- Y
    end associate !-- DE, etc.

    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    call SetFluidKernel &
           ( ProperCell = C % ProperCell, &
              R_E  = GV ( :, G % EDGE_I_U_1 ), &
             Th_E  = GV ( :, G % EDGE_I_U_2 ), &
             Ph_E  = GV ( :, G % EDGE_I_U_3 ), &
              R_W  = GV ( :, G % WIDTH_U_1 ), &
             Th_W  = GV ( :, G % WIDTH_U_2 ), &
             Ph_W  = GV ( :, G % WIDTH_U_3 ), &
              R_C  = GV ( :, G % CENTER_U_1 ), &
              A_1  = A_1, &
              A_3  = A_3, &
              D    = D, &
              M    = LMS % Mass, &
              D_OS = LMS % Density_OS, &
              R_OS = LMS % Radius_OS, &
              AP   = LMS % AtmosphereParameter, &
             nD    = C % nDimensions, &
              N    = FV ( :, F % BARYON_DENSITY_C ), &
              V_1  = FV ( :, F % VELOCITY_U_1 ), &
              V_2  = FV ( :, F % VELOCITY_U_2 ), &
              V_3  = FV ( :, F % VELOCITY_U_3 ) )

    end associate !-- FV, etc.
    end associate !-- C, etc.
    end select !-- A

  end subroutine SetFluid


  function Zero_T ( LMS, T ) result ( F )

    class ( * ), intent ( inout ) :: &
      LMS
    real ( KDR ), intent ( in ) :: &
      T
    real ( KDR ) :: &
      F

    real ( KDR ) :: &
      X_Start, &
      X_Finish, &
      H_Start

    select type ( LMS )
      class is ( LinMestelShuForm )
    associate &
      ( DE  =>  LMS % DifferentialEquation, &
        DF  =>  LMS % DensityFactor )
    associate &
      ( Y  =>  DE % Solution )

    X_Start   =  0.0_KDR
    X_Finish  =  T
    H_Start   =  LMS % TimeScale_OS  *  1.0e-3

    Y ( 1 )  =  1.0_KDR
    Y ( 2 )  =  0.0_KDR
    Y ( 3 )  =  1.0_KDR
    Y ( 4 )  =  0.0_KDR

    call DE % Integrate ( X_Start, X_Finish, H_Start )

    F  =  DF  -  1.0_KDR  /  ( Y ( 1 ) ** 2  *  Y ( 3 ) )

    end associate !-- Y
    end associate !-- DE, etc.
    end select !-- LMS
    
  end function Zero_T


  subroutine SetFluidKernel &
               ( ProperCell, R_E, Th_E, Ph_E, R_W, Th_W, Ph_W, R_C, &
                 A_1, A_3, D, M, D_OS, R_OS, AP, nD, N, V_1, V_2, V_3 )

    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E, Th_E, Ph_E, &
      R_W, Th_W, Ph_W, &
      R_C
    real ( KDR ), intent ( in ) :: &
      A_1, &
      A_3, &
      D, &
      M, &
      D_OS, &
      R_OS, &
      AP
    integer ( KDI ), intent ( in ) :: &
      nD  !-- nDimensions
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      V_1, V_2, V_3

    integer ( KDI ) :: &
      iV, &          !-- iValue
      iS, jS, kS, &  !-- iSubcell
      nV
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF  !-- VolumeFraction
    real ( KDR ) :: &
       R_I,  R_O, &
      Th_I, Th_O, &
      Ph_I, Ph_O, &
       RS,  ThS,  PhS, &  !-- subcell
      dRS, dThS, dPhS, &  !-- subcell
       VS, &  !-- VolumeSubcell
      dVS, &  !-- dVolumeSubcell
      Rho_sq_in_in, Rho_sq_in_out, Rho_sq_out_in, Rho_sq_out_out, &
        Z_sq_in_in,   Z_sq_in_out,   Z_sq_out_in,   Z_sq_out_out, &
      Rho_S_sq, &
        Z_S_sq, &
      Pi

    Pi  =  CONSTANT % PI

    nV  =  size ( N )

    allocate ( BVF ( nV ) )
    call Clear ( BVF )

    nSubCells           =   1
    nSubcells ( : nD )  =  20

    do iV  =  1,  nV

!-- Establish atmosphere at outer radial boundary   
!      if ( .not. ProperCell ( iV ) ) &
!        cycle

       R_I  =   R_E ( iV )
       R_O  =   R_E ( iV )  +   R_W ( iV )

      Th_I  =  Th_E ( iV )
      Th_O  =  Th_E ( iV )  +  Th_W ( iV )

      Ph_I  =  Ph_E ( iV )
      Ph_O  =  Ph_E ( iV )  +  Ph_W ( iV )

      Rho_sq_in_in    =  ( R_I  *  sin ( Th_I ) ) ** 2
      Rho_sq_in_out   =  ( R_O  *  sin ( Th_I ) ) ** 2
      Rho_sq_out_in   =  ( R_I  *  sin ( Th_O ) ) ** 2
      Rho_sq_out_out  =  ( R_O  *  sin ( Th_O ) ) ** 2

        Z_sq_in_in    =  ( R_I  *  cos ( Th_I ) ) ** 2
        Z_sq_in_out   =  ( R_O  *  cos ( Th_I ) ) ** 2
        Z_sq_out_in   =  ( R_I  *  cos ( Th_O ) ) ** 2
        Z_sq_out_out  =  ( R_O  *  cos ( Th_O ) ) ** 2

       if (       Rho_sq_in_in / A_1 ** 2    +  Z_sq_in_in / A_3 ** 2  &
                    <=  1.0_KDR &
            .and. Rho_sq_in_out / A_1 ** 2   +  Z_sq_in_out / A_3 ** 2  &
                    <=  1.0_KDR &
            .and. Rho_sq_out_in / A_1 ** 2   +  Z_sq_out_in / A_3 ** 2  &
                    <=  1.0_KDR &
            .and. Rho_sq_out_out / A_1 ** 2  +  Z_sq_out_out / A_3 ** 2  &
                    <=  1.0_KDR ) &
       then 
         BVF ( iV )  =  1.0_KDR
         cycle
       end if

      if (       Rho_sq_in_in / A_1 ** 2    +  Z_sq_in_in / A_3 ** 2  &
                   >  1.0_KDR &
           .and. Rho_sq_in_out / A_1 ** 2   +  Z_sq_in_out / A_3 ** 2  &
                   >  1.0_KDR &
           .and. Rho_sq_out_in / A_1 ** 2   +  Z_sq_out_in / A_3 ** 2  &
                   >  1.0_KDR &
           .and. Rho_sq_out_out / A_1 ** 2  +  Z_sq_out_out / A_3 ** 2  &
                   >  1.0_KDR ) &
      then 
        BVF ( iV )  =  0.0_KDR
        cycle
      end if
      
       dRS  =  (  R_O  -   R_I )  /  nSubcells ( 1 )
      dThS  =  ( Th_O  -  Th_I )  /  nSubcells ( 2 )
      dPhS  =  ( Ph_O  -  Ph_I )  /  nSubcells ( 2 )

      VS  =  0.0_KDR
      do kS  =  1,  nSubcells ( 3 )
        do jS  =  1,  nSubcells ( 2 )
          do iS  =  1,  nSubcells ( 1 )
             RS  =   R_I  +  ( iS - 0.5_KDR ) *  dRS
            ThS  =  Th_I  +  ( jS - 0.5_KDR ) * dThS
            PhS  =  Ph_I  +  ( kS - 0.5_KDR ) * dPhS
            Rho_S_sq  =  ( RS  *  sin ( ThS ) ) ** 2
              Z_S_sq  =  ( RS  *  cos ( ThS ) ) ** 2
            select case ( nD )
            case ( 2 )
              dVS  =  2 * Pi * RS ** 2  *  sin ( ThS ) &
                        * dRS * dThS
            case ( 3 )
              dVS  =  RS ** 2  *  sin ( ThS ) &
                      *  dRS  *  dThS  *  dPhS
            end select 
            VS  =  VS  +  dVS
            if ( Rho_S_sq / A_1 ** 2 + Z_S_sq / A_3 ** 2  <=  1.0_KDR ) &
              BVF ( iV )  =  BVF ( iV )  +  dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iV )  =  BVF ( iV )  /  VS

    end do !-- iV

    do iV  =  1,  nV
      if ( BVF ( iV )  >  0.0_KDR ) then
        N   ( iV )  =  D  *  BVF ( iV )
        V_1 ( iV )  =  0.0_KDR
      else
        N   ( iV )  =  AP  *  D_OS  *  ( R_C ( iV ) / R_OS ) ** ( -1.5_KDR )
        V_1 ( iV )  =  - sqrt ( 2.0_KDR * M / R_C ( iV ) )
      end if
      V_2 ( iV )  =  0.0_KDR
      V_3 ( iV )  =  0.0_KDR
    end do !-- iV


  end subroutine SetFluidKernel

  
end module LinMestelShu_Form
