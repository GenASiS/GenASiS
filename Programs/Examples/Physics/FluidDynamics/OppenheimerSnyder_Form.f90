module OppenheimerSnyder_Form

  use GenASiS
  use ODE_Solve_Command

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: OppenheimerSnyderForm
    real ( KDR ) :: &
      DensityInitial, &
      RadiusInitial, &
      Density, &
      Radius
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
  end type OppenheimerSnyderForm

    private :: &
      SetFluid, &
      SetReference

      private :: &
        SetFluidKernel, &
        SetReferenceKernel

  class ( OppenheimerSnyderForm ), private, pointer :: &
      OppenheimerSnyderCollapse => null ( )

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


  subroutine Initialize ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ), target :: &
      OS
    character ( * ), intent ( in ) :: &
      Name

    real ( KDR ) :: &
      TimeScale, &
      DensityFactor, &
      RadiusFactor, &
      Beta
    class ( Fluid_D_Form ), pointer :: &
      F
    class ( GeometryFlatForm ), pointer :: &
      G

    if ( OS % Type == '' ) &
      OS % Type = 'a OppenheimerSnyder'

    OppenheimerSnyderCollapse => OS

    call OS % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidCentralCoreForm :: OS % Integrator )
    select type ( FCC => OS % Integrator )
    type is ( FluidCentralCoreForm )
    call FCC % Initialize &
           ( Name, FluidType = 'DUST', &
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

    associate &
      ( R0 => OS % RadiusInitial, &
        D0 => OS % DensityInitial, &
        DF => DensityFactor, &
        RF => RadiusFactor, &
        PI => CONSTANT % PI )

    R0  =  PS % Chart % MaxCoordinate ( 1 ) / 1.1_KDR
    D0  =  1.0_KDR
    DF  =  10.0_KDR
    call PROGRAM_HEADER % GetParameter ( R0, 'RadiusInitial' )
    call PROGRAM_HEADER % GetParameter ( D0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter ( DF, 'DensityFactor' )

    RF = DF ** ( - 1.0_KDR / 3.0_KDR ) 

    TimeScale         =  sqrt ( 3.0 / ( 8.0 * PI * D0 ) )
    Beta              =  acos ( sqrt ( RF ) )
    FCC % FinishTime  =  ( Beta  +  0.5 * sin ( 2 * Beta ) )  *  TimeScale

    call Show ( 'OppenheimerSnyder parameters' )
    call Show ( OS % DensityInitial, 'DensityInitial' )
    call Show ( OS % RadiusInitial, 'RadiusInitial' )
    call Show ( DensityFactor, 'DensityFactor' )
    call Show ( RadiusFactor, 'RadiusFactor' )
    call Show ( PI / 2  *  TimeScale, 'CollapseTime' )
    call Show ( FCC % FinishTime, 'Reset FinishTime' )

    end associate !-- R0, etc.

    F => FA % Fluid_D ( )
    call SetFluid ( OS, F, Time = 0.0_KDR )


   !-- ODE Solver
   OS % DerivativeEvaluator => DerivativeFunction
    
   call OS % Parameters % Initialize ( 2 )
   OS % Parameters % Value ( 1 ) = 1.0_KDR
   OS % Parameters % Value ( 2 ) &
          = OS % DensityInitial * 4.0_KDR / 3.0_KDR * CONSTANT % PI

   allocate ( OS % ODEF )
    
   call OS % ODEF % Initialize &
                      ( OS % Parameters, &
                        OS % DerivativeEvaluator )
   
   !-- Reference Solution

   allocate ( OS % Reference )
   allocate ( OS % Difference )

   call OS % Reference % Initialize &
           ( PS, NameShortOption = 'Reference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
   call OS % Difference % Initialize &
           ( PS, NameShortOption = 'Difference', &
             FluidType = 'DUST', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FCC
    nullify ( F, G )

  end subroutine Initialize


  subroutine ComputeError ( OS )
    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

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
    
    select type ( FCC => OS % Integrator )
    type is ( FluidCentralCoreForm )
    select type ( FA => FCC % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    F_R => OppenheimerSnyderCollapse % Reference % Fluid_D ( )
    call SetReferenceKernel ( OppenheimerSnyderCollapse, F_R, FCC % Time )

    F_D => OppenheimerSnyderCollapse % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )


    select type ( PS => OS % Integrator % PositionSpace )
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


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    call OS % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( OS, F, Time )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time
    
    real ( KDR ) :: &
      Radius, &
      Density
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => OS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    if ( Time == 0.0_KDR ) then
      OS % Radius   =  OS % RadiusInitial
      OS % Density  =  OS % DensityInitial
    end if

    call SetFluidKernel &
           (    R = G % Value ( :, G % CENTER_U ( 1 ) ), &
             dR_L = G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             dR_R = G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             Density = OS % Density, &
             RadiusDensity = OS % Radius, &
              N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

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

    F_R => OppenheimerSnyderCollapse % Reference % Fluid_D ( )
    call SetReferenceKernel ( OppenheimerSnyderCollapse, F_R, FCC % Time )

    F_D => OppenheimerSnyderCollapse % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    F_D % Value = abs ( F_D % Value )

    end select !-- FA
    end select !-- FCC
    nullify ( F, F_R, F_D )

  end subroutine SetReference


  subroutine SetFluidKernel &
               ( R, dR_L, dR_R, Density, RadiusDensity, N, VX, VY, VZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      dR_L, &
      dR_R
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      VX, VY, VZ

    N  = 0.0_KDR
    VX = 0.0_KDR
    VY = 0.0_KDR
    VZ = 0.0_KDR

    associate &
      ( R_In  => R - dR_L, &
        R_Out => R + dR_R )

    where ( R_Out <= RadiusDensity )
      N = Density
    end where
    where ( R_In < RadiusDensity .and. R_Out > RadiusDensity )
      N = Density * ( RadiusDensity ** 3  -  R_In ** 3 ) &
                    / ( R_Out ** 3  -  R_In ** 3 )
    end where

!    N = Density / ( 1 + exp ( ( R - RadiusDensity ) &
!                              / ( 3 * ( dR_L + dR_R ) ) ) ) 
 
    end associate !-- R_In, etc.

  end subroutine SetFluidKernel


  subroutine SetReferenceKernel ( OS, F, Time )
    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      Time

    integer ( KDI ) :: &
      iC
    real ( KDR ) :: &
      X_1, &
      X_2, &
      H, &
      R_New
    real ( KDR ), dimension ( 2 ) :: &
      Y_Start
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => OS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )
    
    select type ( C => PS % Chart ) 
    class is ( Chart_SLD_Form )

    associate &
      ( Rho   => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        R     => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Rho_0 => OS % DensityInitial, &
        R0    => OS % RadiusInitial )

    X_1 = 0.0_KDR
    X_2 = Time

    H = X_2 / 2

    Y_Start ( 1 ) = 1.0_KDR
    Y_Start ( 2 ) = 0.0_KDR

    if ( Time > 0.0_KDR ) then 
      !-- Solve system of ODEs
      call OS % ODEF % IntegrateODE &
                         ( Y_Start, X_1, X_2, H, &
                           OS % ODEF % RequestedAccuracy, &
                           RK4Option = .false. ) 
    end if
 
    R_New = R0 * Y_Start ( 1 ) 
    
    Rho = 0.0_KDR
 
    where ( R <= R_New ) 

      Rho = Rho_0 / ( Y_Start ( 1 ) ** 3 )

    end where

    end associate !-- Rho, etc.
    end select !-- C
    end select !-- PS
    nullify ( G )

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
      SqrtTiny

    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    select type ( P => Parameters ) 
    type is ( Real_1D_Form )
    associate &
      ( G => P % Value ( 1 ), &
        M => P % Value ( 2 ) )

    dYdX ( 1 ) = Y ( 2 )
    dYdX ( 2 ) = - G * M / max ( Y ( 1 ) ** 2, SqrtTiny )

    end associate
    end select !-- P
  
  end subroutine DerivativeFunction


end module OppenheimerSnyder_Form
