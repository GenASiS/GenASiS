module SpheroidCollapse_Form

! From Lin, Mestel, and Shu 1965
! "The Gravitational Collapse of a Uniform Spheroid"

  use GenASiS
  use ODE_Solve_Command

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: SpheroidCollapseForm
    Real ( KDR ) :: &
      SemiMajor, &
      SemiMinor, &
      Eccentricity, &
      InitialDensity
    real ( KDR ), dimension ( : ), allocatable :: &
      R0, &
      Z0
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
  end type SpheroidCollapseForm

    private :: &
      SetFluid, &
      SetReference
      
      private :: &
        SetFluidKernel, &
        SetReferenceKernel

  class ( SpheroidCollapseForm ), private, pointer :: &
      SpheroidCollapse => null ( )

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


  subroutine Initialize ( SC, Name )

    class ( SpheroidCollapseForm ), intent ( inout ), target :: &
      SC
    character ( * ), intent ( in ) :: &
      Name

    class ( Fluid_D_Form ), pointer :: &
      F
    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ) :: &
      TimeScale

    if ( SC % Type == '' ) &
      SC % Type = ' a SpheroidCollapse'

    SpheroidCollapse => SC

    call SC % InitializeTemplate ( Name )

    !-- Integrator

    allocate ( FluidCentralCoreForm :: SC % Integrator )
    select type ( FCC => SC % Integrator )
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
   G => PS % Geometry ( )

   select type ( FA => FCC % Current_ASC )
   class is ( Fluid_ASC_Form )

    !-- Initial conditions
    SC % Eccentricity   = 0.6_KDR
    SC % SemiMajor      = 8.0_KDR !PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    SC % InitialDensity = 1.0_KDR

    call PROGRAM_HEADER % GetParameter ( SC % Eccentricity, 'Eccentricity' )
    call PROGRAM_HEADER % GetParameter ( SC % InitialDensity, 'rho_0' )
    call PROGRAM_HEADER % GetParameter ( SC % SemiMajor, 'SemiMajor' )

    SC % SemiMinor = SC % SemiMajor * sqrt ( 1.0_KDR - SC % Eccentricity ** 2 )
    
    TimeScale &
      = sqrt ( 3.0_KDR / ( 8.0_KDR * CONSTANT % PI * SC % InitialDensity ) )

    FCC % FinishTime = 1.51_KDR * TimeScale
    call Show ( FCC % FinishTime, 'Reset FinishTime' )

    F => FA % Fluid_D ( )
    call SetFluid ( SC, F, Time = 0.0_KDR )

   !-- ODE Solver
   SC % DerivativeEvaluator => DerivativeFunction
    
   call SC % Parameters % Initialize ( 3 )
   SC % Parameters % Value ( 1 ) = 2 * CONSTANT % PI
   SC % Parameters % Value ( 2 ) = SC % Eccentricity
   SC % Parameters % Value ( 3 ) = SC % InitialDensity

   allocate ( SC % ODEF )
    
   call SC % ODEF % Initialize &
                      ( SC % Parameters, &
                        SC % DerivativeEvaluator )
   
   !-- Reference Solution

   associate &
     ( R     => G % Value ( :, G % Center_U ( 1 ) ), &
       Theta => G % Value ( :, G % Center_U ( 2 ) ) )

   allocate &
     ( SC % R0 ( size ( R ) ), &
       SC % Z0 ( size ( R ) ) ) 

   SC % R0 = R * sin ( Theta )
   SC % Z0 = R * cos ( Theta )

   allocate ( SC % Reference )
   allocate ( SC % Difference )

   call SC % Reference % Initialize &
           ( PS, NameShortOption = 'Reference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
   call SC % Difference % Initialize &
           ( PS, NameShortOption = 'Difference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    !-- Cleanup

    end associate !-- Rho
    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F, G )

  end subroutine Initialize


  subroutine ComputeError( SC )
    class ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC

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
    
    select type ( FCC => SC % Integrator )
    type is ( FluidCentralCoreForm )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    F_R => SpheroidCollapse % Reference % Fluid_D ( )
    call SetReferenceKernel ( SpheroidCollapse, F_R, FCC % Time )

    F_D => SpheroidCollapse % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )


    select type ( PS => SC % Integrator % PositionSpace )
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


  impure elemental subroutine Finalize ( SC )

    type ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC 
    
    call SC % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( SC, F, Time )

    class ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => SC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )

    call SetFluidKernel &
           (  PS, &
              a_1 = SC % SemiMajor, &
              a_3 = SC % SemiMinor, &
              D0  = SC % InitialDensity, &
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

    F_R => SpheroidCollapse % Reference % Fluid_D ( )
    call SetReferenceKernel ( SpheroidCollapse, F_R, FCC % Time )

    F_D => SpheroidCollapse % Difference % Fluid_D ( )
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

  
  subroutine SetReferenceKernel ( SC, F, Time )
    class ( SpheroidCollapseForm ), intent ( inout ) :: &
      SC
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
      R, &   !-- Cylindrical coordinates R ** 2 = X ** 2 + Y ** 2
      Z

    select type ( PS => SC % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )
    
    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )

    associate &
      ( Rho   => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        R_xy  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Rho_0 => SC % InitialDensity, &
        R0    => SC % R0, &
        Z0    => SC % Z0, &
        a     => SC % SemiMajor, &
        c     => SC % SemiMinor )

    X_1 = 0.0_KDR
    X_2 = Time

    H = X_2 / 10

    Y_Start ( 1 ) = 1.0_KDR
    Y_Start ( 2 ) = 0.0_KDR
    Y_Start ( 3 ) = 1.0_KDR
    Y_Start ( 4 ) = 0.0_KDR

    if ( Time > 0.0_KDR ) then
      !-- Solve system of ODEs
      call SC % ODEF % IntegrateODE &
                         ( Y_Start, X_1, X_2, &
                           SC % ODEF % RequestedAccuracy, &
                           H, RK4Option = .false. ) 
    end if
    
    allocate &
      ( R ( size ( G % Value ( :, G % CENTER_U ( 1 ) ) ) ), &
        Z ( size ( G % Value ( :, G % CENTER_U ( 1 ) ) ) ) )

    New_a = Y_Start ( 1 ) * a
    New_c = Y_Start ( 3 ) * c 
    
    Rho = 0.0_KDR

    R = R_xy * sin ( Theta )
    Z = R_xy * cos ( Theta )
 
    where ( ( R / New_a) ** 2 + ( Z / New_c ) ** 2 <= 1.0_KDR ) 

      Rho = Rho_0 / abs ( ( Y_Start ( 1 ) ** 2 * Y_Start ( 3 ) ) )

    end where

    end associate !-- Rho, etc.
    end select !-- C
    end select !-- PS
    deallocate ( R, Z )

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

end module SpheroidCollapse_Form


