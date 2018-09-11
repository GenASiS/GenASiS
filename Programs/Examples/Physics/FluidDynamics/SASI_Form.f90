module SASI_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SASIForm
    real ( KDR ), private :: &
      CentralMass, &
      ShockRadius, &
      AdiabaticIndex, &
      PolytropicConstant, &
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
        Kappa  => SF % PolytropicConstant, &
        FourPi => 4 * CONSTANT % PI )

    M      =  1.2_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )
    R_S    =  200.0_KDR * UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( R_S, 'ShockRadius' )
    Gamma  =  4.0_KDR / 3.0_KDR 
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    AR     =  0.36_KDR * UNIT % SOLAR_MASS / UNIT % SECOND
    call PROGRAM_HEADER % GetParameter ( AP, 'AccretionRate' )

    R_Out  = 2 * R_S
    R_In   = 0.1_KDR * R_S

    SF % PolytropicConstant &
      = 2 / ( Gamma + 1.0_KDR ) &
          * ( ( Gamma - 1.0_KDR ) / ( Gamma + 1.0_KDR ) ) ** Gamma &
          ( AR / FourPi ) ** ( 1.0_KDR - Gamma ) &
          ( 2 * CONSTANT % GRAVITATIONAL * M ) ** ( 1.0_KDR + Gamma ) / 2 &
          R_S ** ( ( 3.0_KDR * Gamma - 5.0_KDR ) / 2 )

    call Show ( 'SASI parameters' )
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( R_S,   UNIT % KILOMETER, 'ShockRadius' )
    call Show ( AR, UNIT % SOLAR_MASS / UNIT % SECOND, 'AccretionRate' )
    call Show ( Gamma, 'Gamma' )
    call Show ( Kappa, 'PolytropicConstant' )

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
    real ( KDR ) :: &
      u
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

    call SetVelocityProfile ( SF, u )

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
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ) )

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    where ( R < R_S )
      N    = AR / ( FourPi * u * R ** 2 ) / CONSTANT % ATOMIC_MASS_UNIT
      V_1  = - u
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

  
  subroutine SetVelocityProfile ( SF, u )
    type ( SASIForm ), intent ( inout ) :: &
      SF
    real ( KDR ), dimension ( : ), intent ( out ) ::
      u

    logical ( KDL ) :: &
      Converged
    integer ( KDI ) :: &
      nIterations
    real ( KDR ) :: &
      a
    real ( KDR ), allocatable, dimension ( : ) :: & 
      u1, &
      u2, &
      F1, &
      F2, &
      root
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( M      => SF % CentralMass, &
        R_S    => SF % RadiusShock, &
        Gamma  => SF % AdiabaticIndex, &
        Kappa  => SF % PolytropicConstant, &
        AR     => SF % AccretionRate, &
        G_C     => CONSTANT % GRAVITATIONAL, &
        FourPi => 4 * CONSTANT % PI )

    a = 2 * Gamma * Kappa / ( Gamma - 1.0_KDR ) &
          * ( AR / FourPi ) ** ( Gamma - 1.0_KDR )

   select type ( FCE => SF % Integrator )
   class is ( FluidCentralExcisionForm )
   
   select type ( PS => FCE % PositionSpace )
   class is ( Atlas_SC_Form )
   G => PS % Geometry ( )
   
   associate ( R => G % Value ( :, G % CENTER_U ( 1 ) ) )

   allocate &
     ( u1 ( size ( R ) ), &
       u2 ( size ( R ) ), &
       F1 ( size ( R ) ), & 
       F2 ( size ( R ) ), &
       root ( size ( R ) ) )
        

   u1 = ( (Gamma - 1.0_KDR ) * G_C * M / ( Gamma * Kappa ) ) &
          ** ( 1.0_KDR / ( 1.0_KDR - Gamma ) ) * AR / FourPi

   u2 = ( Gamma - 1.0_KDR ) / ( Gamma + 1.0_KDR ) &
          * sqrt ( 2 * G_C * M / R_S )

   nIterations = 0
   
   Converged = .false.

   do while ( .not. converged )
     nIterations = nIterations + 1

     !--- Use secant method (see mathworld.com) to find 
     !      root of Bernoulli equation: 
     
     where ( R < R_S )
       F1 = R * u1 ** 2 &
            + a * r ** ( 3.0_KDR - 2 * Gamma ) * u1 ** ( 1.0_KDR - Gamma )
            - 2 * G_C * M
       F2 = R * u2 ** 2 &
            + a * r ** ( 3.0_KDR - 2 * Gamma ) * u2 ** ( 1.0_KDR - Gamma )
            - 2 * G_C * M
       root = u1 - F1 / ( F1 - F2 ) * (u1 - u2 )
     elsewhere
       root = 1.0_KDR
     end where

     if ( all ( abs ( ( root - u1 ) / root ) <= 1.0e-8_KDR ) ) then
       Converged = .true.
     else
        u2 = u1 
        u1 = root
     end if   
   end do

   u = root

   end associate !-- R
   end select !-- G
   end select !-- FCE
   end associate !-- M, etc.

  end subroutine SetVelocityProfile

end module SASI_Form
