module SASI_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SASIForm
    real ( KDR ) :: &
      CentralMass, &
      ShockRadius, &
      AdiabaticIndex, &
      PolytropicConstant, &
      AccretionRate, &
      MachNumber
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SASIForm

    private :: &
      SetFluid, &
      SetVelocityProfile, &
      Apply_BC_CSL_PowerLaw
        
        private :: &
          IntroduceDensityRingsPerturbation

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
        Kappa  => SF % PolytropicConstant, &
        AR     => SF % AccretionRate, &
        M_N    => SF % MachNumber, &
        R_In   => RadiusIn, &
        R_Out  => RadiusOuter, &
        G_C    => CONSTANT % GRAVITATIONAL, &
        FourPi => 4 * CONSTANT % PI )

    M      =  1.2_KDR * UNIT % SOLAR_MASS
    call PROGRAM_HEADER % GetParameter ( M, 'CentralMass' )
    R_S    =  200.0_KDR * UNIT % KILOMETER
    call PROGRAM_HEADER % GetParameter ( R_S, 'ShockRadius' )
    Gamma  =  4.0_KDR / 3.0_KDR 
    call PROGRAM_HEADER % GetParameter ( Gamma, 'AdiabaticIndex' )
    AR     =  0.36_KDR * UNIT % SOLAR_MASS / UNIT % SECOND
    call PROGRAM_HEADER % GetParameter ( AR, 'AccretionRate' )
    M_N    = 3.0e2_KDR
    call PROGRAM_HEADER % GetParameter ( M_N, 'MachNumber' )

    R_Out  = 2 * R_S
    R_In   = 0.2_KDR * R_S

    Kappa = 2 / ( Gamma + 1.0_KDR ) &
              * ( ( Gamma - 1.0_KDR ) / ( Gamma + 1.0_KDR ) ) ** Gamma &
              * ( AR / FourPi ) ** ( 1.0_KDR - Gamma ) &
              * ( 2 * G_C * M ) ** ( ( 1.0_KDR + Gamma ) / 2 ) &
              * R_S ** ( ( 3.0_KDR * Gamma - 5.0_KDR ) / 2 )

    call Show ( 'SASI parameters' )
    call Show ( M,     UNIT % SOLAR_MASS, 'M' )
    call Show ( R_In,  UNIT % KILOMETER, 'R_In' )
    call Show ( R_Out, UNIT % KILOMETER, 'R_Out' )
    call Show ( R_S,   UNIT % KILOMETER, 'Shock Radius' )
    call Show ( AR, UNIT % SOLAR_MASS / UNIT % SECOND, 'Accretion Rate' )
    call Show ( M_N, 'Mach Number' )
    call Show ( Gamma, 'Gamma' )
    call Show ( Kappa, 'Polytropic Constant' )

    end associate !-- M, etc.

    FinishTime  =  150.0e-3_KDR * UNIT % SECOND


    !-- Integrator

    allocate ( FluidCentralExcisionForm :: SF % Integrator )
    select type ( FCE => SF % Integrator )
    type is ( FluidCentralExcisionForm )
    call FCE % Initialize &
           ( Name, FluidType = 'IDEAL', &
             GeometryType = 'NEWTONIAN', &
             UseCustomBoundaryInnerOption = .true., &
             FinishTimeOption = FinishTime, &
             LimiterParameterOption = 1.8_KDR, &
             RadiusMaxOption = RadiusOuter, &
             RadiusMinOption = RadiusIn, &
             CentralMassOption = SF % CentralMass )

    !-- Custom boundary condition
    
    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_CE_Form )
    PS % Apply_BC_CSL_Custom => Apply_BC_CSL_PowerLaw

    !-- Initial conditions

    call SetFluid ( SF )

    !call IntroduceDensityRingsPerturbation ( SF )

    !-- Cleanup

    end select !-- PS
    end select !-- FCE

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
      iV  !-- iValue
    real ( KDR ) :: &
      B, &
      C
    real ( KDR ), allocatable, dimension ( : ) :: &
      u
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    associate &
      ( M      => SF % CentralMass, &
        R_S    => SF % ShockRadius, &
        Gamma  => SF % AdiabaticIndex, &
        Kappa  => SF % PolytropicConstant, &
        AR     => SF % AccretionRate, &
        M_N    => SF % MachNumber, &
        G_C    => CONSTANT % GRAVITATIONAL, &
        FourPi => 4 * CONSTANT % PI )

    select type ( FCE => SF % Integrator )
    class is ( FluidCentralExcisionForm )

    select type ( FA => FCE % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_P_I ( )

    !-- Fluid

    call F % SetAdiabaticIndex ( Gamma )
    
    select type ( PS => FCE % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (     N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
            P => F % Value ( :, F % PRESSURE ), & 
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
            R => G % Value ( :, G % CENTER_U ( 1 ) ) )

    allocate ( u ( size ( R ) ) )

    call SetVelocityProfile ( SF, u )

    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    do iV = 1, size ( R )
      if ( R ( iV ) < R_S ) then 
        V_1 ( iV )  = ( - u ( iV ) )
        N ( iV )    = AR / ( FourPi * abs ( V_1 ( iV ) ) * R ( iV ) ** 2 )
        P ( iV )    = Kappa * N ( iV ) ** Gamma
      else
        V_1 ( iV )  = - sqrt ( 2.0_KDR * G_C * M / R ( iV ) )
        N ( iV )    = AR / ( FourPi * abs ( V_1 ( iV ) ) * R ( iV ) ** 2 )
        P ( iV )    = N ( iV ) / Gamma * ( V_1 ( iV ) / M_N ) ** 2
      end if

      if ( R ( iV ) <= 0.21 * R_S .and. R ( iV ) > 0.19 * R_S ) then
         B = P ( iV ) * R ( iV ) ** 4
         C = N ( iV ) * R ( iV ) ** 3
      end if
    end do
    
    do iV = 1, size ( R )
      if ( R ( iV ) < 0.2*R_S ) then
         P ( iV ) = B * R ( iV ) ** ( - 4 )
         N ( iV ) = C * R ( iV ) ** ( -3 )
      end if
   end do

   E = P / ( Gamma - 1.0_KDR )

   N = N / CONSTANT % ATOMIC_MASS_UNIT 

    call F % ComputeFromPrimitive ( G )

    end associate !-- N, etc.
    end select !-- PS
    end select !-- FA
    end select !-- FCE
    end associate !-- M, etc.
    nullify ( F, G )

  end subroutine SetFluid

  
  subroutine SetVelocityProfile ( SF, u )
    type ( SASIForm ), intent ( inout ) :: &
      SF
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      u

    logical ( KDL ) :: &
      Converged
    integer ( KDI ) :: &
      nIterations, &
      iV
    real ( KDR ) :: &
      a, &
      F_MP, &
      F_U
    real ( KDR ), allocatable, dimension ( : ) :: & 
      u_U, &
      u_L, &
      MP, &
      root
    class ( GeometryFlatForm ), pointer :: &
      G

    associate &
      ( M      => SF % CentralMass, &
        R_S    => SF % ShockRadius, &
        Gamma  => SF % AdiabaticIndex, &
        Kappa  => SF % PolytropicConstant, &
        AR     => SF % AccretionRate, &
        G_C    => CONSTANT % GRAVITATIONAL, &
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
     ( u_U ( size ( R ) ), &
       u_L ( size ( R ) ), &
       root ( size ( R ) ), &
       MP   ( size ( R ) ) )
  
   u_L = 1.0_KDR
   u_U = 1.0_KDR
   MP  = 1.0_KDR

   where ( R < R_S )
     u_L = 100.0_KDR * UNIT % KILOMETER / UNIT % SECOND
     u_U = ( Gamma - 1.0_KDR ) / ( Gamma + 1.0_KDR ) &
            * sqrt ( 2 * G_C * M / R_S )
   end where

   nIterations = 0
   
   Converged = .false.
   root = 1.0_KDR

   !do while ( .not. converged )
     
   !  call Show ( nIterations, '>>> n1' )
     !--- Use Bisection method to find 
     !      root of Bernoulli equation: 
     do iV = 1, size ( R )
       Converged = .false.
       nIterations = 0
       do while ( .not. Converged )
         nIterations = nIterations + 1
         if ( R ( IV ) < R_S ) then
           MP ( iV )  = ( u_U ( iV ) + u_L ( iV ) ) / 2
           F_MP = R ( iV ) * MP ( iV ) ** 2 &
                  + a * R ( iV ) ** ( 3.0_KDR - 2 * Gamma ) &
                  * MP ( iV ) ** ( 1.0_KDR - Gamma ) &
                  - 2 * G_C * M
           F_U = R ( iV ) * u_U ( iV ) ** 2 &
                  + a * R ( iV ) ** ( 3.0_KDR - 2 * Gamma ) &
                  * u_U ( iV ) ** ( 1.0_KDR - Gamma ) &
                  - 2 * G_C * M
           if ( sign ( 1.0_KDR , F_MP ) == sign ( 1.0_KDR, F_U ) ) then
             u_U ( iV ) = MP ( iV )
           else
             u_L ( iV ) = MP ( iV )
           end if
           if ( u_U ( iV ) - u_L ( iV ) <= 1.0e-8_KDR ) then
             root ( iV ) = MP ( iV )
             Converged = .true.
           end if
         else
          Converged = .true. 
         end if
       end do
     end do
   !  if ( all ( ( u_U - u_L ) <= 1.0e-8_KDR ) ) then
   !    root = MP
   !    Converged = .true.
   !  end if
   !end do


   u = root

   end associate !-- R
   end select !-- G
   end select !-- FCE
   end associate !-- M, etc.
call Show ( 'end of Velocity' )
  end subroutine SetVelocityProfile 


  subroutine Apply_BC_CSL_PowerLaw ( CSL, F, iDimension, iConnection )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      CSL
    class ( StorageForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      !-- FIXME: hard coding OUTFLOW for the moment
      iS, &  !-- iSelected
      iF     !-- iField
    class ( GeometryFlatForm ), pointer :: &
      G

    !-- FIXME: hard coding OUTFLOW for the moment

    do iS = 1, F % nVariables
      iF = F % iaSelected ( iS )
      call CSL % CopyBoundary ( F % Value ( :, iF ), iDimension, iConnection )
    end do !-- iS
    
    !-- Begin PowerLaw

    select type ( PS => CSL % Atlas )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( F )
    class is ( Fluid_P_I_Form )

    select type ( CSL )
    class is ( Chart_SLD_Form )

      associate &
        ( Connectivity => CSL % Atlas % Connectivity, &
          iD => iDimension, &
          iC => iConnection )
      if ( trim ( CSL % CoordinateSystem ) /= 'SPHERICAL' ) then
        call Show ( 'PowerLaw only implemented for spherical coordinates', &
                    CONSOLE % ERROR )
        call Show ( 'SASI_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Apply_BC_CSL_PowerLaw', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iD /= 1 ) then
        call Show ( 'PowerLaw only implemented for radial direction', &
                    CONSOLE % ERROR )
        call Show ( 'SASI_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Apply_BC_CSL_PowerLaw', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iC /= Connectivity % iaInner ( iD ) ) then
        call Show ( 'PowerLaw only implemented for inner boundary', &
                    CONSOLE % ERROR )
        call Show ( 'SASI_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Apply_BC_CSL_PowerLaw', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( CSL % iaBrick ( iD ) /= 1 ) &
        return
      end associate !-- Connectivity, etc.

call Show ( '>>> Implement power law' )
call Show ( '>>>>>> See e.g. Chart_SL_Template' ) 
call Show ( '>>>>>>>>> CopyBoundaryTemplate' )
call Show ( '>>>>>>>>> CopyBoundaryKernel' )
!    call F % ComputeFromPrimitive ( G )

    class default
      call Show ( 'Chart type not implemented', CONSOLE % ERROR )
      call Show ( 'SASI_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Apply_BC_CSL_PowerLaw', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- CSL

    class default
      call Show ( 'Power law boundary only for ideal gas', CONSOLE % ERROR )
      call Show ( 'SASI_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Apply_BC_CSL_PowerLaw', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- F

    end select !-- PS
    nullify ( G )

  end subroutine Apply_BC_CSL_PowerLaw


  subroutine IntroduceDensityRingsPerturbation ( SF )
    type ( SASIForm ), intent ( inout ) :: &
      SF

    integer ( KDI ) :: &
      iV
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    associate &
      ( M      => SF % CentralMass, &
        R_S    => SF % ShockRadius, &
        Gamma  => SF % AdiabaticIndex, &
        Kappa  => SF % PolytropicConstant, &
        AR     => SF % AccretionRate, &
        M_N    => SF % MachNumber, &
        G_C    => CONSTANT % GRAVITATIONAL, &
        FourPi => 4 * CONSTANT % PI )

   select type ( FCE => SF % Integrator )
   class is ( FluidCentralExcisionForm )
   
   select type ( PS => FCE % PositionSpace )
   class is ( Atlas_SC_Form )
   G => PS % Geometry ( )

   select type ( FA => FCE % Current_ASC )
   class is ( Fluid_ASC_Form )
   F => FA % Fluid_P_I ( )
   
   associate &
      (     N => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
            P => F % Value ( :, F % PRESSURE ), & 
          V_1 => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
          V_2 => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
          V_3 => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
            E => F % Value ( :, F % INTERNAL_ENERGY ), &
            R => G % Value ( :, G % CENTER_U ( 1 ) ) )
   
   do iV = 1, size ( R )
     if ( R ( iV ) > 1.1_KDR * R_S .and. R ( iV ) < 1.15_KDR * R_S ) then
        N ( iV )   = 1.2_KDR * N ( iV )
     !   V_1 ( iV ) = - AR / ( FourPi * N ( iV ) * R ( iV ) ** 2 )
        P ( iV )   = N ( iV ) * CONSTANT % ATOMIC_MASS_UNIT &
                      / Gamma * ( V_1 ( iV ) / M_N ) ** 2
        E ( iV )   = P ( iV ) / ( Gamma - 1.0_KDR )
     end if 

     if ( R ( iV ) > 1.5_KDR * R_S .and. R ( iV ) < 1.55_KDR * R_S ) then
        N ( iV )  = 1.3_KDR * N ( iV )
       ! V_1 ( iV ) = - AR / ( FourPi * N ( iV ) * R ( iV ) ** 2 )
        P ( iV )  = N ( iV ) * CONSTANT % ATOMIC_MASS_UNIT &
                      / Gamma * ( V_1 ( iV ) / M_N ) ** 2
        E ( iV )  = P ( iV ) / ( Gamma - 1.0_KDR )
     end if 
     
   end do

   end associate !-- N, etc
   end select !-- FA
   end select !-- PS
   end select !-- FCE
   end associate !-- M, etc

   nullify ( F, G )

  end subroutine IntroduceDensityRingsPerturbation

end module SASI_Form
