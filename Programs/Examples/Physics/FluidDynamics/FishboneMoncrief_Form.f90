module FishboneMoncrief_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_CE_Form ) :: FishboneMoncriefForm
    real ( KDR ), private :: &
      CentralMass, &
      AngularMomentumParameter, &
      RadiusInner, &
      RadiusOuter, &
      SpecificAngularMomentum, &
      AdiabaticIndex, &
      DensityMax, &
      AtmosphereParameter
  contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type FishboneMoncriefForm

    private :: &
      InitializeUniverse, &
      SetInitial

      private :: &
        SetFluid

contains


  subroutine Initialize_H ( U, Name )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name

    if ( U % Type == '' ) &
      U % Type = 'a FishboneMoncrief'

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


  impure elemental subroutine Finalize ( FM )
    
    type ( FishboneMoncriefForm ), intent ( inout ) :: &
      FM

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( FishboneMoncriefForm ), intent ( in ) :: &
      U

    call U % Universe_F_CE_Form % ShowParameters ( )

    call Show ( U % CentralMass, &
                UNIT % SOLAR_MASS, &
                'CentralMass' )
    call Show ( U % AngularMomentumParameter, &
                'AngularMomentumParameter' )
    call Show ( U % RadiusInner, &
                UNIT % KILOMETER,  &
                'RadiusInner' )
    call Show ( U % RadiusOuter, &
                UNIT % KILOMETER,  &
                'RadiusOuter' )
    call Show ( U % SpecificAngularMomentum, &
                UNIT % KILOMETER * UNIT % SPEED_OF_LIGHT, &
                'SpecificAngularMomentum' )
    call Show ( U % AdiabaticIndex, &
                'AdiabaticIndex' )
    call Show ( U % DensityMax, &
                UNIT % MASS_DENSITY_CGS, &
                'DensityMax' )
    call Show ( U % AtmosphereParameter, &
                'AtmosphereParameter' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( FM, Name )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      FM
    character ( * ), intent ( in )  :: &
      Name

    real ( KDR ) :: &
      G, &
      c, &
      R_Min, &
      T_Finish

     G  =  CONSTANT % GRAVITATIONAL
     c  =  CONSTANT % SPEED_OF_LIGHT

    associate &
      ( M      =>  FM % CentralMass, &
        Kappa  =>  FM % AngularMomentumParameter, &  !-- Between 1 and 2
        R_In   =>  FM % RadiusInner, &
        R_Out  =>  FM % RadiusOuter, &
        L      =>  FM % SpecificAngularMomentum )

    M      =  3.0_KDR * UNIT % SOLAR_MASS
    Kappa  =  1.85_KDR
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )
    call PROGRAM_HEADER % GetParameter ( Kappa, 'AngularMomentumParameter' )

    R_In   =  6.0_KDR * G * M  /  c ** 2
    call PROGRAM_HEADER % GetParameter ( R_In,  'RadiusInner' )

    R_Out  =  R_In * Kappa / ( 2.0_KDR - Kappa )
    L      =  sqrt ( Kappa * G * M * R_In )

    R_Min  =  ( 2.0_KDR / 3.0_KDR ) * R_In
    call PROGRAM_HEADER % GetParameter ( R_Min, 'RadiusMin' )

    T_Finish  =  1.0e-3_KDR  *  UNIT % SECOND

    call FM % Initialize &
           ( FluidType = 'IDEAL', &
             GravitationType = 'NEWTON_CM', &
             Name = Name, &
             FinishTimeOption = T_Finish, &
             RadiusMaxOption = R_Out, &
             RadiusExcisionOption = R_Min, &
             CentralMassOption = M )

    end associate !-- M, etc.

    FM % Integrator % SetInitial    =>  SetInitial
    FM % Integrator % System        =>  FM

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( FM  =>  I % System )
      class is ( FishboneMoncriefForm )
    select type ( I => FM % Integrator )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )

    associate &
      ( Gamma   =>  FM % AdiabaticIndex, &
        RhoMax  =>  FM % DensityMax, &
        AP      =>  FM % AtmosphereParameter )

    Gamma   =  1.4_KDR
    RhoMax  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS
    AP      =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter ( Gamma,  'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( RhoMax, 'DensityMax' )
    call PROGRAM_HEADER % GetParameter ( AP,     'AtmosphereParameter' )

    call F % SetAdiabaticIndex &
           ( Gamma )

    call SetFluid ( FM, F )

    if ( allocated ( FM % SA_Fluid ) ) then
      select type ( F_SA  =>  FM % SA_Fluid % FieldSet_SA )
      class is ( Fluid_P_I_Form )
        call F_SA % SetAdiabaticIndex &
               ( Gamma )
        call F_SA % SetFiducialParameters &
               ( FiducialBaryonDensity = F % FiducialBaryonDensity, &
                 FiducialPressure = F % FiducialPressure )
      end select !-- F_SA
    end if

    if ( allocated ( FM % AA_Fluid ) ) then
      select type ( F_AA  =>  FM % AA_Fluid % FieldSet_AA )
      class is ( Fluid_P_I_Form )
        call F_AA % SetAdiabaticIndex &
               ( Gamma )
        call F_AA % SetFiducialParameters &
               ( FiducialBaryonDensity = F % FiducialBaryonDensity, &
                 FiducialPressure = F % FiducialPressure )
      end select !-- F_AA
    end if

    end associate !-- Gamma, etc.
    end select !-- F
    end select !-- I
    end select !-- FM

  end subroutine SetInitial


  subroutine SetFluid ( FM, F )

    class ( FishboneMoncriefForm ), intent ( inout ), target :: &
      FM
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    real ( KDR ) :: &
      EnthalpyMax, &
      PolytropicParameter, &
      EnergyDensityMin, &
      GC, &
      c, &
      amu
    real ( KDR ), dimension ( : ), allocatable :: &
      Enthalpy

    associate &
      ( Kappa   =>  FM % AngularMomentumParameter, &
        M       =>  FM % CentralMass, &
        R_In    =>  FM % RadiusInner, &
        Gamma   =>  FM % AdiabaticIndex, &
        RhoMax  =>  FM % DensityMax, &
        AP      =>  FM % AtmosphereParameter, &
        W_Max   =>  EnthalpyMax, &
        K       =>  PolytropicParameter, &
        E_Min   =>  EnergyDensityMin )

     GC  =  CONSTANT % GRAVITATIONAL
      c  =  CONSTANT % SPEED_OF_LIGHT
    amu  =  CONSTANT % ATOMIC_MASS_UNIT

    W_Max = GC * M / ( c ** 2  *  R_In )  &
            *  ( 0.5_KDR * ( Kappa  +  1.0_KDR / Kappa )  -  1.0_KDR )

    K  =  ( Gamma - 1.0_KDR ) / Gamma  &
          *  W_Max  /  RhoMax ** ( Gamma - 1.0_KDR )

    associate &
      ( G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    allocate ( Enthalpy ( size ( FV, dim = 1 ) ) )

    associate &
      (     W => Enthalpy, &
            N => FV ( :, F % BARYON_DENSITY_C ), &
          V_1 => FV ( :, F % VELOCITY_U_1 ), &
          V_2 => FV ( :, F % VELOCITY_U_2 ), &
          V_3 => FV ( :, F % VELOCITY_U_3 ), &
            E => FV ( :, F % ENERGY_DENSITY_C ), &
            R => GV ( :, G % CENTER_U_1 ), &
        Theta => GV ( :, G % CENTER_U_2 ) )

    !-- Set atmosphere to determine fluid min parameters

    E_Min  =  1.0e-20_KDR  *  UNIT % ENERGY_DENSITY_NUCLEAR
    call PROGRAM_HEADER % GetParameter ( E_Min, 'EnergyDensityMin' )

    N  =  AP  *  RhoMax / amu  *  ( R / R_In ) ** ( - 1.5_KDR )
    E  =  E_Min

    call F % SetBaryonDensityMin ( )
    call F % SetEnergyDensityMin ( )
    call F % SetFiducialParameters &
           ( FiducialBaryonDensity = AP * RhoMax / amu, &
             FiducialPressure = E_Min * ( Gamma - 1.0_KDR ) )

    N  =  0.0_KDR
    E  =  0.0_KDR

    !-- Set disk

    W  =  max ( 0.0_KDR, &
                GC * M / ( c ** 2  *  R_In )  &
                * ( R_In / R  -  1.0_KDR  +  0.5_KDR * Kappa &
                    -  0.5_KDR * Kappa  *  R_In ** 2  &
                       /  ( R * sin ( Theta ) ) ** 2 ) )

    N  =  ( ( Gamma - 1.0_KDR ) / ( Gamma * K ) * W ) &
          ** ( 1.0_KDR / ( Gamma - 1.0_KDR ) )

    V_1  =  0.0_KDR
    V_2  =  0.0_KDR

    !-- Reset atmosphere

    where ( N  >  amu  *  F % BaryonDensityMin )
      V_3  =  sqrt ( Kappa * GC * M * R_In )  /  ( R * sin ( Theta ) ) ** 2
      E    =  K  *  N ** Gamma  /  ( Gamma - 1.0_KDR )
    elsewhere
      N    =    AP  *  RhoMax  *  ( R / R_In ) ** ( - 1.5_KDR )
      V_1  =  - sqrt ( 2.0_KDR * GC * M / R )
      V_3  =    0.0_KDR
      E    =    E_Min
    end where

    N  =  N / amu

    end associate !-- W, etc.
    end associate !-- FV, etc.
    end associate !-- G
    end associate !-- Kappa, etc.

  end subroutine SetFluid


end module FishboneMoncrief_Form
