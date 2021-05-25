module LinMestelShu_Form

! From Lin, Mestel, and Shu 1965
! "The Gravitational Collapse of a Uniform Spheroid"

  use GenASiS
  use ODE_Solve_Command

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: LinMestelShuForm
    Real ( KDR ) :: &
      Eccentricity, &
      SemiMajor, &
      SemiMinor, &
      DensityInitial
    type ( Real_1D_Form ) :: &
      Parameters
    class ( ODEForm ), allocatable :: &
      ODEF
    type ( Fluid_ASC_Form ), allocatable :: &
      Reference, &
      Difference
    procedure ( DerivativeInterface ), public, pointer, nopass :: &
      DerivativeEvaluator => null () 
  contains
    procedure, public, pass :: &
      Initialize 
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type LinMestelShuForm

    private :: &
      SetFluid, &
      SetReference
      
      private :: &
        SetFluidKernel, &
        SetReferenceKernel

  class ( LinMestelShuForm ), private, pointer :: &
      LinMestelShu => null ( )

  interface 
    subroutine DerivativeInterface ( Parameters, X, Y, dYdX )
      use GenASiS
      class ( * ), intent ( in ) :: &
        Parameters
      real ( KDR ), intent ( in ) :: &
        X
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        Y
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        dYdX
    end subroutine DerivativeInterface
  end interface

contains


  subroutine Initialize ( LMS, Name )

    class ( LinMestelShuForm ), intent ( inout ), target :: &
      LMS
    character ( * ), intent ( in ) :: &
      Name

    class ( Fluid_D_Form ), pointer :: &
      F
    character ( LDN ):: &
      OS_Parameter
    real ( KDR ) :: &
      RadiusOppenheimerSnyder, & 
      DensityFactor, &
      RadiusFactor, &
      TimeScale, &
      Beta

    if ( LMS % Type == '' ) &
      LMS % Type = ' a LinMestelShu'

    LinMestelShu => LMS

    call LMS % InitializeTemplate ( Name )

    !-- Integrator

    allocate ( FluidCentralCoreForm :: LMS % Integrator )
    select type ( FCC => LMS % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'DUST' ,&
             GeometryType = 'NEWTONIAN', &
             DimensionlessOption = .true., &
             GravityFactorOption = 0.01_KDR, &
             LimiterParameterOption = 1.0_KDR )

    FCC % SetReference => SetReference

   select type ( PS => FCC % PositionSpace )
   class is ( Atlas_SC_Form )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )

   associate &
      ( R0   => LMS % SemiMajor, &
        Z0   => LMS % SemiMinor, &
        e0   => LMS % Eccentricity, &
        D0   => LMS % DensityInitial, &
        R_OS => RadiusOppenheimerSnyder, &
        DF   => DensityFactor, &
        RF   => RadiusFactor, &
        PI   => CONSTANT % PI )

   !-- Initial conditions
   e0 = 0.9_KDR
   R0 = PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR

   call PROGRAM_HEADER % GetParameter ( e0, 'Eccentricity' )
   call PROGRAM_HEADER % GetParameter ( R0, 'SemiMajor' )

   Z0 = R0 * sqrt ( 1.0_KDR - e0 ** 2 )

   call PROGRAM_HEADER % GetParameter ( OS_Parameter, 'OS_Parameter' )

   select case ( trim ( OS_Parameter ) )
   case ( 'RADIUS' )
     R_OS = R0
     !--   Find the density needed to have the same mass as the sphere in
     !-- OppenheimerSnyder (OS) collapse for a spheroid with an initial 
     !-- semi-major axis equal to the inital radius of the OS sphere.
     D0   = sqrt ( 1.0_KDR - e0 ** 2 )
   case default
     !-- default to density of unity or user input
     OS_Parameter = 'DENSITY'
     D0 = 1.0_KDR
     call PROGRAM_HEADER % GetParameter ( D0, 'DensityInitial' )
     !-- Find the radius of the OS sphere with M=3/(4piR0^3) and rho=1.0
     R_OS = R0 * ( 1.0_KDR - e0 ** 2 ) ** ( 3.0_KDR / 2.0_KDR )
   end select

   DF = 10.0_KDR                               

   RF = DF ** ( - 1.0_KDR / 3.0_KDR )
    
   TimeScale = sqrt ( 3.0_KDR / ( 8.0_KDR * PI * D0 ) )

   Beta = acos ( sqrt ( RF ) )

   FCC % FinishTime = ( Beta + sin ( 2 * Beta ) / 2 ) * TimeScale
    
   call Show ( 'LinMestelShu parameters' )
   call Show ( LMS % DensityInitial, 'DensityInitial' )
   call Show ( 3.0_KDR / ( 4 * PI * R0 ** 2 * Z0 ), 'Mass' )
   call Show ( LMS % SemiMajor, 'SemiMajorInitial' )
   call Show ( LMS % SemiMinor, 'SemiMinorInitial' )
   call Show ( LMS % Eccentricity, 'EccentricityInitial' )
   call Show ( OS_Parameter, 'OS_Parameter' )
   call Show ( R_OS, 'OS_Radius (M)' )
   call Show ( DensityFactor, 'DensityFactor' )
   call Show ( RadiusFactor, 'RadiusFactor' )
   call Show ( PI / 2  *  TimeScale, 'CollapseTime' )
   call Show ( FCC % FinishTime, 'Reset FinishTime' )

   end associate !-- R0, etc.
   
   F => FA % Fluid_D ( )
   call SetFluid ( LMS, F, Time = 0.0_KDR )

   !-- ODE Solver
   LMS % DerivativeEvaluator => DerivativeFunction
    
   call LMS % Parameters % Initialize ( 3 )
   LMS % Parameters % Value ( 1 ) = 2 * CONSTANT % PI
   LMS % Parameters % Value ( 2 ) = LMS % Eccentricity
   LMS % Parameters % Value ( 3 ) = LMS % DensityInitial

   allocate ( LMS % ODEF )
    
   call LMS % ODEF % Initialize &
                      ( LMS % Parameters, &
                        LMS % DerivativeEvaluator )
   
   !-- Reference Solution


   allocate ( LMS % Reference )
   allocate ( LMS % Difference )

   call LMS % Reference % Initialize &
           ( PS, NameShortOption = 'Reference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
   call LMS % Difference % Initialize &
           ( PS, NameShortOption = 'Difference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F )

  end subroutine Initialize


  subroutine ComputeError ( LMS )
    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS

    integer ( KDI ) :: &
      nM, &
      iV, &
      iE, &
      iM !-- iMessage
    real ( KDR ) :: &
      L1
    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference         
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( FCC => LMS % Integrator )
    type is ( FluidCentralCoreForm )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    F_R => LinMestelShu % Reference % Fluid_D ( )
    call SetReferenceKernel ( LinMestelShu, F_R, FCC % Time )

    F_D => LinMestelShu % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )


    select type ( PS => LMS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    
    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )
    
    associate ( D   => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
                R   => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ) )


    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )

    CO % Outgoing % Value ( 1 ) &
           = sum ( abs ( R ), mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 )&
           = sum ( abs ( D ), mask = C % IsProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    end associate !-- D, etc. 
    
    !-- L1 error

    associate ( IN   => CO % Incoming % Value )

    L1 = IN ( 2 ) / IN ( 1 )
    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- IN 

    end select !-- C
    end select !-- PS

    end select !-- FA
    end select !-- FCC
    nullify ( F, F_R, F_D )


  end subroutine


  impure elemental subroutine Finalize ( LMS )

    type ( LinMestelShuForm ), intent ( inout ) :: &
      LMS 
    
    call LMS % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( LMS, F, Time )

    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => LMS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )

    call SetFluidKernel &
           (  PS, &
              a_1 = LMS % SemiMajor, &
              a_3 = LMS % SemiMinor, &
              D0  = LMS % DensityInitial, &
              N   = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
              V_1 = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
              V_2 = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
              V_3 = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )
    call C % ExchangeGhostData ( F )

    end select    !-- C
    end select    !-- PS
    nullify ( G )

  end subroutine SetFluid 


  subroutine SetReference ( FCC )
    class ( IntegratorTemplate ), intent ( inout ) :: &
      FCC 

    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference

    select type ( FCC )
    class is ( FluidCentralCoreForm )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    F_R => LinMestelShu % Reference % Fluid_D ( )
    call SetReferenceKernel ( LinMestelShu, F_R, FCC % Time )

    F_D => LinMestelShu % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    F_D % Value = abs ( F_D % Value )

    end select !-- FA
    end select !-- FCC
    nullify ( F, F_R, F_D )

  end subroutine SetReference

  
  subroutine SetFluidKernel ( A, a_1, a_3, D0, N, V_1, V_2, V_3 )

    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      a_1, &
      a_3, &
      D0
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      V_1, &
      V_2, &
      V_3

    class ( GeometryFlatForm ), pointer :: &
      G
    integer ( KDI ) :: &
      i, &
      iC, &  !-- iCell
      iS, jS, kS     !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS      !--  X_Subcell
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF, &  !-- VolumeFraction
      rho_sq, & !-- X^2 + Y^2 
      Z_sq
    real ( KDR ) :: &
      dVS, &  !-- dVolumeSubcell
      VS, &   !-- VolumeSubcell
      rho_sq_in_in, &
      rho_sq_in_out, &
      rho_sq_out_in, &
      rho_sq_out_out, &
      Z_sq_in_in, &
      Z_sq_in_out, &
      Z_sq_out_in, &
      Z_sq_out_out, &
      rho_S_sq, &
      Z_s_sq

    N   = 0.0_KDR
    V_1 = 0.0_KDR
    V_2 = 0.0_KDR
    V_3 = 0.0_KDR

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )

    associate ( dV => G % Value ( :, G % VOLUME ) )

    allocate ( BVF ( size ( dV ) ) )
    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : C % nDimensions ) = 20

    do iC = 1, size ( dV )
      associate &
        (  X => &
             G % Value ( iC, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) ), &
          dX_L => &
            G % Value ( iC, G % WIDTH_LEFT_U ( 1 ) &
                              : G % WIDTH_LEFT_U ( 3 ) ), &
          dX_R => &
            G % Value ( iC, G % WIDTH_RIGHT_U ( 1 ) &
                              : G % WIDTH_RIGHT_U ( 3 ) ), &
          Pi => CONSTANT % PI )

      if ( .not. C % IsProperCell ( iC ) ) cycle

       X_I  =  X  - dX_L
       X_O  =  X  + dX_R

       rho_sq_in_in   = ( X_I ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
       rho_sq_in_out  = ( X_O ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
       rho_sq_out_in  = ( X_I ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2
       rho_sq_out_out = ( X_O ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2

       Z_sq_in_in   = ( X_I ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
       Z_sq_in_out  = ( X_O ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
       Z_sq_out_in  = ( X_I ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2
       Z_sq_out_out = ( X_O ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2

       if ( rho_sq_in_in / a_1 ** 2 + Z_sq_in_in / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_in_out / a_1 ** 2 &
                    + Z_sq_in_out / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_out_in / a_1 ** 2 &
                    + Z_sq_out_in / a_3 ** 2 <= 1.0_KDR &
            .and. rho_sq_out_out / a_1 ** 2 &
                    + Z_sq_out_out / a_3 ** 2 <= 1.0_KDR ) then 
         BVF ( iC ) = 1.0_KDR
         cycle
       end if

      if ( rho_sq_in_in / a_1 ** 2 + Z_sq_in_in / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_in_out / a_1 ** 2 &
                   + Z_sq_in_out / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_out_in / a_1 ** 2 &
                   + Z_sq_out_in / a_3 ** 2 > 1.0_KDR &
           .and. rho_sq_out_out / a_1 ** 2 &
                   + Z_sq_out_out / a_3 ** 2 > 1.0_KDR ) then 
        BVF ( iC ) = 0.0_KDR
        cycle
      end if
      
      dXS  =  ( X_O  -  X_I ) / nSubcells

      VS = 0.0_KDR
      do kS = 1, nSubcells ( 3 )
        do jS = 1, nSubcells ( 2 )
          do iS = 1, nSubcells ( 1 )
            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            rho_S_sq = ( XS ( 1 ) * sin ( XS ( 2 ) ) ) ** 2
            Z_S_sq = ( XS ( 1 ) * cos ( XS ( 2 ) ) ) ** 2
            select case ( C % nDimensions )
            case ( 2 )
              dVS = 2 * Pi * XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                      * dXS ( 1 ) * dXS ( 2 )
            case ( 3 )
              dVS = XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                     * dXS ( 1 ) * dXS ( 2 ) * dXS ( 3 )
            end select 
            VS = VS + dVS
            if ( rho_S_sq / a_1 ** 2 + Z_S_sq / a_3 ** 2 <= 1.0_KDR ) &
              BVF ( iC ) = BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

      end associate !-- X, etc.
    end do !-- iC

    end associate !-- dV
    end select !-- C

    nullify ( G )

    N = D0 * BVF

  end subroutine SetFluidKernel

  
  subroutine SetReferenceKernel ( LMS, F, Time )
    class ( LinMestelShuForm ), intent ( inout ) :: &
      LMS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

    class ( GeometryFlatForm ), pointer :: &
      G
    integer ( KDI ) :: &
      iC
    real ( KDR ) :: &
      X_1, &
      X_2, &
      New_a,& 
      New_c, &
      H
    real ( KDR ), dimension ( 4 ) :: &
      Y_Start
    real ( KDR ), dimension ( : ), allocatable :: &
      R_XY, &   !-- Cylindrical coordinates R ** 2 = X ** 2 + Y ** 2
      Z

    select type ( PS => LMS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )
    
    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )

    associate &
      ( Rho   => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        R     => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Rho_0 => LMS % DensityInitial, &
        a     => LMS % SemiMajor, &
        c     => LMS % SemiMinor )

    X_1 = 0.0_KDR
    X_2 = Time

    H = X_2 / 10

    Y_Start ( 1 ) = 1.0_KDR
    Y_Start ( 2 ) = 0.0_KDR
    Y_Start ( 3 ) = 1.0_KDR
    Y_Start ( 4 ) = 0.0_KDR

    if ( Time > 0.0_KDR ) then
      !-- Solve system of ODEs
      call LMS % ODEF % IntegrateODE &
                         ( Y_Start, X_1, X_2, &
                           LMS % ODEF % RequestedAccuracy, &
                           H, RK4Option = .false. ) 
    end if
    

    New_a = Y_Start ( 1 ) * a
    New_c = Y_Start ( 3 ) * c 
    
    Rho = 0.0_KDR

    R_XY = R * sin ( Theta )
    Z    = R * cos ( Theta )
 
    where ( ( R_XY / New_a) ** 2 + ( Z / New_c ) ** 2 <= 1.0_KDR ) 

      Rho = Rho_0 / ( Y_Start ( 1 ) ** 2 * Y_Start ( 3 ) )

    end where

    end associate !-- Rho, etc.
    end select !-- C
    end select !-- PS

  end subroutine SetReferenceKernel


  subroutine DerivativeFunction ( Parameters, X, Y, dYdX )
    class ( * ), intent ( in ) :: &
      Parameters
    real ( KDR ), intent ( in ) :: &
      X
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      dYdX

    real ( KDR ) :: &
      e, &
      A, &
      C, &
      sqrtTiny

    sqrtTiny = 1.0e-10_KDR

    select type ( P => Parameters ) 
    type is ( Real_1D_Form )
    associate &
      ( TwoPi => P % Value ( 1 ), &
        e0    => P % Value ( 2 ), &
        rho_0 => P % Value ( 3 ) )

    e &
      = max ( sqrt ( max ( 1.0_KDR &
                             - ( Y ( 3 ) / Y ( 1 ) ) ** 2 &
                               * ( 1.0_KDR - e0 ** 2 ), &
                           SqrtTiny ) ), SqrtTiny )
    A &
      = TwoPi * sqrt ( ( 1.0_KDR - max ( e ** 2, SqrtTiny ) ) ) &
          / max ( e ** 3, SqrtTiny) &
          * ( asin ( e ) - e * sqrt ( ( 1 - max ( e ** 2, SqrtTiny ) ) ) )

    C &
      = 2 * TwoPi / max ( e ** 2, SqrtTiny ) &
          * ( 1.0_KDR - sqrt ( 1.0_KDR - max ( e ** 2, SqrtTiny ) ) &
                          * asin ( e ) / e )

    dYdX ( 1 ) = Y ( 2 )
    dYdX ( 3 ) = Y ( 4 )

    if ( ( Y ( 1 ) * Y ( 3 ) ) > 0.0_KDR ) then
      dYdX ( 2 ) = - rho_0 * A &
                       / max ( ( Y ( 1 ) * Y ( 3 ) ), sqrtTiny )
    else
      dYdX ( 2 ) = - rho_0 * A &
                       / min ( ( Y ( 1 ) * Y ( 3 ) ), - sqrtTiny )
    end if

    dYdX ( 4 ) = - rho_0 * C &
                     / max ( ( Y ( 1 ) ** 2 ), sqrtTiny )

    end associate
    end select !-- P
  
  end subroutine DerivativeFunction

end module LinMestelShu_Form


