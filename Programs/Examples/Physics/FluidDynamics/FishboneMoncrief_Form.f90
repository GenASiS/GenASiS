module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: FishboneMoncriefForm
    real ( KDR ), private :: &
      AngularMomentumParameter, &
      CentralMass, &
      RadiusInner, &
      AdiabaticIndex, &
      DensityMax, &
      AtmosphereParameter
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type FishboneMoncriefForm

    private :: &
      SetFluid

contains


  subroutine Initialize ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      RadiusMin, &
      RadiusOuter, &
      SpecificAngularMomentum, &
      FinishTime


    if ( FM % Type == '' ) &
      FM % Type = 'a FishboneMoncrief'

    call FM % InitializeTemplate ( Name )

    associate &
      ( Kappa  => FM % AngularMomentumParameter, &  !-- Between 1 and 2
        M      => FM % CentralMass, &
        R_In   => FM % RadiusInner, &
        Gamma  => FM % AdiabaticIndex, &
        N_Max  => FM % DensityMax, &
        AP     => FM % AtmosphereParameter, &
        R_Out  => RadiusOuter, &
        L      => SpecificAngularMomentum, &
        G      => CONSTANT % GRAVITATIONAL, &
        c      => CONSTANT % SPEED_OF_LIGHT, &
        amu    => CONSTANT % ATOMIC_MASS_UNIT )

    Kappa  =  1.85_KDR
    M      =  3.0_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( Kappa, 'AngularMomentumParameter' )
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )

    R_In   =  6.0_KDR * G * M  /  c ** 2
    Gamma  =  1.4_KDR
    N_Max  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS / amu
    AP     =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter ( R_In,  'RadiusInner' )
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( N_Max, 'DensityMax' )
    call PROGRAM_HEADER % GetParameter ( AP,    'AtmosphereParameter' )

    R_Out  =  R_In * Kappa / ( 2.0_KDR - Kappa )
    L      =  sqrt ( Kappa * G * M * R_In )

    call Show ( 'FishboneMoncrief parameters' )
    call Show ( Kappa, 'Kappa' ) 
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( L,     UNIT % KILOMETER * UNIT % SPEED_OF_LIGHT, 'L' )
    call Show ( Gamma, 'Gamma' )
    call Show ( N_Max, UNIT % FEMTOMETER ** (-3), 'N_Max' )
    call Show ( amu * N_Max, UNIT % MASS_DENSITY_CGS, 'Rho_Max' )
    call Show ( AP, 'AtmosphereParameter' )

    RadiusMin  =  ( 2.0_KDR / 3.0_KDR ) * R_In
    call PROGRAM_HEADER % GetParameter ( RadiusMin, 'RadiusMin' )
    end associate !-- Kappa, etc.

    FinishTime  =  1.0e-3_KDR * UNIT % SECOND


    !-- Integrator

    allocate ( FluidCentralExcisionForm :: FM % Integrator )
    select type ( FCE => FM % Integrator )
    type is ( FluidCentralExcisionForm )
    call FCE % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN', &
             FinishTimeOption = FinishTime, &
             LimiterParameterOption = 1.8_KDR, &
             RadiusMaxOption = RadiusOuter, &
             RadiusMinOption = RadiusMin, &
             CentralMassOption = FM % CentralMass, &
             nCellsPolarOption = 128 )
    end select !-- FCE


    !-- Initial conditions

    call SetFluid ( FM )


  end subroutine Initialize


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    call FM % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( FM )

    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

    integer ( KDI ) :: &
      iB  !-- iBoundary
    real ( KDR ) :: &
      EnthalpyMax, &
      PolytropicParameter
    real ( KDR ), dimension ( : ), allocatable :: &
      Enthalpy
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    associate &
      ( Kappa => FM % AngularMomentumParameter, &
        M     => FM % CentralMass, &
        R_In  => FM % RadiusInner, &
        Gamma => FM % AdiabaticIndex, &
        N_Max => FM % DensityMax, &
        AP    => FM % AtmosphereParameter, &
        W_Max => EnthalpyMax, &
        K     => PolytropicParameter, &
        GC    => CONSTANT % GRAVITATIONAL, &
        c     => CONSTANT % SPEED_OF_LIGHT, &
        amu   => CONSTANT % ATOMIC_MASS_UNIT )

    select type ( FCE => FM % Integrator )
    class is ( FluidCentralExcisionForm )

    select type ( FA => FCE % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    !-- Fluid

    call F % SetAdiabaticIndex ( Gamma )

    W_Max = GC * M / ( c ** 2  *  R_In )  &
            *  ( 0.5_KDR * ( Kappa  +  1.0_KDR / Kappa )  -  1.0_KDR )

    K  =  ( Gamma - 1.0_KDR ) / Gamma  &
          *  W_Max  /  ( amu * N_Max ) ** ( Gamma - 1.0_KDR )

    allocate ( Enthalpy ( F % nValues ) )

    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     W => Enthalpy, &
            N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
            R => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ) )

    W = max ( 0.0_KDR, &
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

    call F % ComputeFromPrimitive ( G )

    end associate !-- W, etc.
    end select !-- PS
    end select !-- FA
    end select !-- FCE
    end associate !-- Kappa, etc.
    nullify ( F, G )

  end subroutine SetFluid


end module FishboneMoncrief_Form
