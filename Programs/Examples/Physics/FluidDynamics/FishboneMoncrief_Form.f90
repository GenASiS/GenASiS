module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidCentralExcisionForm ) :: FishboneMoncriefForm
    real ( KDR ), private :: &
      AngularMomentumParameter, &
      CentralMass, &
      RadiusInner, &
      AdiabaticIndex, &
      DensityMax, &
      AtmosphereParameter
  contains
    procedure, private, pass :: &
      Initialize_FM
    generic, public :: &
      Initialize => Initialize_FM
    final :: &
      Finalize
  end type FishboneMoncriefForm

    private :: &
      InitializeFluidCentralExcision, &
      SetProblem

      private :: &
        SetFluid, &
        SetBaryonDensityMin

        private :: &
          SetFluidKernel

contains


  subroutine Initialize_FM ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    character ( * ), intent ( in ) :: &
      Name

    if ( FM % Type == '' ) &
      FM % Type = 'a FishboneMoncrief'

    call InitializeFluidCentralExcision ( FM, Name )
    call SetProblem ( FM )

  end subroutine Initialize_FM


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

  end subroutine Finalize


  subroutine InitializeFluidCentralExcision ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      R_Min, &
      R_Out, &
      FinishTime, &
      SpecificAngularMomentum, &
      G, c

    associate &
      ( Kappa => FM % AngularMomentumParameter, &  !-- Between 1 and 2
        M     => FM % CentralMass, &
        R_In  => FM % RadiusInner, &
        L     => SpecificAngularMomentum )

    G  =  CONSTANT % GRAVITATIONAL
    c  =  CONSTANT % SPEED_OF_LIGHT

    Kappa  =  1.85_KDR
    M      =  3.0_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( Kappa, 'AngularMomentumParameter' )
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )

    R_In   =  6.0_KDR * G * M  /  c ** 2
    call PROGRAM_HEADER % GetParameter ( R_In,  'RadiusInner' )

    R_Min  =  ( 2.0_KDR / 3.0_KDR ) * R_In
    R_Out  =  R_In * Kappa / ( 2.0_KDR - Kappa )
    L      =  sqrt ( Kappa * G * M * R_In )

    FinishTime  =  1.0e-3_KDR * UNIT % SECOND

    call FM % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', Name = Name, &
             FinishTimeOption = FinishTime, &
             LimiterParameterOption = 1.8_KDR, &
             RadiusMaxOption = R_Out, &
             RadiusMinOption = R_Min, &
             CentralMassOption = FM % CentralMass, &
             nCellsPolarOption = 128 )

    call Show ( 'FishboneMoncrief parameters' )
    call Show ( Kappa, 'Kappa' ) 
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( L,     UNIT % KILOMETER * UNIT % SPEED_OF_LIGHT, 'L' )

    end associate !-- Kappa, etc.

  end subroutine InitializeFluidCentralExcision


  subroutine SetProblem ( FM )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    real ( KDR ) :: &
      amu
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => FM % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    associate &
      ( Gamma => FM % AdiabaticIndex, &
        N_Max => FM % DensityMax, &
        AP    => FM % AtmosphereParameter )

    amu  =  CONSTANT % ATOMIC_MASS_UNIT

    Gamma  =  1.4_KDR
    N_Max  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS / amu
    AP     =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( N_Max, 'DensityMax' )
    call PROGRAM_HEADER % GetParameter ( AP,    'AtmosphereParameter' )

    call Show ( Gamma, 'Gamma' )
    call Show ( N_Max, UNIT % FEMTOMETER ** (-3), 'N_Max' )
    call Show ( amu * N_Max, UNIT % MASS_DENSITY_CGS, 'Rho_Max' )
    call Show ( AP, 'AtmosphereParameter' )

    end associate !-- Gamma, etc.

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( FM, F, G )

    call SetBaryonDensityMin ( FM, F )

    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine SetFluid ( FM, F, G )

    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    real ( KDR ) :: &
      EnthalpyMax, &
      PolytropicParameter, &
      GC, c, amu
    real ( KDR ), dimension ( : ), allocatable :: &
      Enthalpy

    associate &
      ( Kappa => FM % AngularMomentumParameter, &
        M     => FM % CentralMass, &
        R_In  => FM % RadiusInner, &
        Gamma => FM % AdiabaticIndex, &
        N_Max => FM % DensityMax, &
        AP    => FM % AtmosphereParameter, &
        W_Max => EnthalpyMax, &
        K     => PolytropicParameter )

    GC   =  CONSTANT % GRAVITATIONAL
    c    =  CONSTANT % SPEED_OF_LIGHT
    amu  =  CONSTANT % ATOMIC_MASS_UNIT

    call F % SetAdiabaticIndex ( Gamma )

    W_Max = GC * M / ( c ** 2  *  R_In )  &
            *  ( 0.5_KDR * ( Kappa  +  1.0_KDR / Kappa )  -  1.0_KDR )

    K  =  ( Gamma - 1.0_KDR ) / Gamma  &
          *  W_Max  /  ( amu * N_Max ) ** ( Gamma - 1.0_KDR )

    allocate ( Enthalpy ( F % nValues ) )

    call SetFluidKernel &
           ( G % Value ( :, G % CENTER_U ( 1 ) ), &
             G % Value ( :, G % CENTER_U ( 2 ) ), &
             Kappa, M, R_In, Gamma, N_Max, K, AP, GC, c, amu, &
             Enthalpy, &
             F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             F % Value ( :, F % INTERNAL_ENERGY ), &
             F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

    end associate !-- Kappa, etc.

  end subroutine SetFluid


  subroutine SetBaryonDensityMin ( FM, F )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( I => FM % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( PSC => PS % Chart )
    class is ( Chart_SLD_Form )

    call CO % Initialize &
           ( PS % Communicator, nOutgoing = [ 1 ], nIncoming = [ 1 ] )

    associate &
      ( My_N_Min => CO % Outgoing % Value ( 1 ), &
           N_Min => CO % Incoming % Value ( 1 ) )
 
    My_N_Min  =  minval ( F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                          mask = PSC % IsProperCell )

    call CO % Reduce ( REDUCTION % MIN )

    call F % SetBaryonDensityMin ( N_Min )

    end associate !-- My_N_Min, etc.
    end select !-- PSC
    end select !-- PS
    end select !-- I

  end subroutine SetBaryonDensityMin


  subroutine SetFluidKernel &
               ( R, Theta, Kappa, M, R_In, Gamma, N_Max, K, AP, GC, c, amu, &
                 W, N, E, V_1, V_2, V_3 )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      Theta
    real ( KDR ), intent ( in ) :: &
      Kappa, &
      M, &
      R_In, &
      Gamma, &
      N_Max, &
      K, &
      AP, &
      GC, c, amu
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      W, &
      N, &
      E, &
      V_1, V_2, V_3

    W  =  max ( 0.0_KDR, &
                GC * M / ( c ** 2  *  R_In )  &
                * ( R_In / R  -  1.0_KDR  +  0.5_KDR * Kappa &
                    -  0.5_KDR * Kappa  *  R_In ** 2  &
                       /  ( R * sin ( Theta ) ) ** 2 ) )
    N  =  ( ( ( Gamma - 1.0_KDR ) / ( Gamma * K ) * W ) &
          ** ( 1.0_KDR / ( Gamma - 1.0_KDR ) ) ) &
          / amu

    E  =  K  *  ( amu * N ) ** Gamma  /  ( Gamma - 1.0_KDR )

    V_1 = 0.0_KDR
    V_2 = 0.0_KDR

    where ( N > 0.0_KDR )
      V_3  =  sqrt ( Kappa * GC * M * R_In )  /  ( R * sin ( Theta ) ) ** 2
    elsewhere
      N    =    AP * N_Max  *  ( R / R_In ) ** ( - 1.5_KDR )
      V_1  =  - sqrt ( 2.0_KDR * GC * M / R )
      V_3  =    0.0_KDR
    end where

  end subroutine SetFluidKernel


end module FishboneMoncrief_Form
