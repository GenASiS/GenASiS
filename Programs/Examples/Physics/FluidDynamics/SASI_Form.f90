module SASI_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SASIForm
    real ( KDR ), private :: &
      CentralMass, &
      ShockRadius, &
      AdiabaticIndex, &
      AccretionRate
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SASIForm

    private :: &
      SetFluid

contains


  subroutine Initialize ( SF, Name )

    class ( SASIForm ), intent ( inout ) :: &
      SF
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      RadiusIn, &
      RadiusOuter, &
      FinishTime


    if ( SF % Type == '' ) &
      SF % Type = 'a SASI'

    call SF % InitializeTemplate ( Name )

    associate &
      ( M      => SF % CentralMass, &
        R_S    => SF % ShockRadius, &
        Gamma  => SF % AdiabaticIndex, &
        AR     => SF % AccretionRate, &
        R_In   => RadiusIN, &
        R_Out  => RadiusOuter, &
        Pi     => CONSTANT % PI )

    M      =  3.0_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )
    R_S    =  200.0_KDR * UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( R_S, 'ShockRadius' )
    Gamma  =  1.3_KDR 
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    AR     =  0.3_KDR * UNIT % SOLAR_MASS / UNIT % SECOND
    call PROGRAM_HEADER % GetParameter ( AP, 'AccretionRate' )

    R_Out  = 2 * R_S
    R_In   = 0.1_KDR * R_S

    call Show ( 'SASI parameters' )
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( R_S,   UNIT % KILOMETER, 'ShockRadius' )
    call Show ( AR, UNIT % SOLAR_MASS / UNIT % SECOND, 'AccretionRate' )
    call Show ( Gamma, 'Gamma' )

    end associate !-- M, etc.

    FinishTime  =  1.0e-3_KDR * UNIT % SECOND


    !-- Integrator

    allocate ( FluidCentralExcisionForm :: SF % Integrator )
    select type ( FCE => SF % Integrator )
    type is ( FluidCentralExcisionForm )
    call FCE % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN', &
             FinishTimeOption = FinishTime, &
             LimiterParameterOption = 1.8_KDR, &
             RadiusMaxOption = RadiusOuter, &
             RadiusMinOption = RadiusIn, &
             CentralMassOption = SF % CentralMass, &
             nCellsPolarOption = 128 )
    end select !-- FCE


    !-- Initial conditions

    call SetFluid ( SF )


  end subroutine Initialize


  impure elemental subroutine Finalize ( SF )
    
    type ( SASIForm ), intent ( inout ) :: &
      SF

    call SF % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( SF )

    type ( SASIForm ), intent ( inout ) :: &
      SF

    integer ( KDI ) :: &
      iB  !-- iBoundary
    real ( KDR ) :: &
      Kappa, &
      u_sr, &
      rho_sr, &
      p_sr
    class ( RootFindingForm ) :: &
      Root
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    associate &
      ( M     => SF % CentralMass, &
        R_S   => SF % RadiusShock, &
        Gamma => SF % AdiabaticIndex, &
        AR    => SF % AccretionRate, &
        GC    => CONSTANT % GRAVITATIONAL, &
        FourPi = > 4 * CONSTANT % PI )

    select type ( FCE => SF % Integrator )
    class is ( FluidCentralExcisionForm )

    select type ( FA => FCE % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    !-- Fluid

    call F % SetAdiabaticIndex ( Gamma )

    u_sr = - ( Gamma - 1.0_KDR ) / ( Gamma + 1.0_KDR ) * AR / FourPi &
              * sqrt ( 2 * GC * M ) * R_S ** ( - 1.0_KDR / 2 )
    rho_sr = ( Gamma + 1.0_KDR ) / ( Gamma - 1.0_KDR ) * AR / FourPi &
               * 1.0_KDR / sqrt ( 2 * GC * M ) * R_S ** ( - 3.0_KDR / 2 )
    p_sr = - 2.0_KDR / ( Gamma + 1.0_KDR ) * AR / FourPi &
               * sqrt ( 2 * GC * M ) * R_S ** ( - 5.0_KDR / 2 )

    Kappa = p_sr / rho_sr ** ( Gamma )
    

    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
           AR => SF % AccretionRate, &
            R => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
        FourPi)

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    where ( R > R_S )
      V_1  =  sqrt ( 2 * Kappa * GC * M * R_In )  / R  
      N    =  AR / ( 
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


end module SASI_Form
