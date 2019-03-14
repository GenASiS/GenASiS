module OppenheimerSnyder_Form

  !-- For example, Misner, Thorne, Wheeler p. 663, Eqs. (25.28)-(25.29)

  use GenASiS
  use ODE_Solve_Command

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: OppenheimerSnyderForm
    real ( KDR ) :: &
      RadiusInitial, &
      DensityInitial, &
      TimeScale
!      Density, &
!      Radius
  !   type ( Real_1D_Form ) :: &
  !     Parameters
  !   class ( ODEForm ), allocatable :: &
  !     ODEF
    type ( Fluid_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  !   procedure ( DerivativeInterface ), public, pointer, nopass :: &
  !     DerivativeEvaluator => null () 
  contains
    procedure, private, pass :: &
      Initialize_OS
    generic, public :: &
      Initialize => Initialize_OS
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type OppenheimerSnyderForm

    private :: &
!       SetReference
      InitializeFluidCentralCore, &
      InitializeDiagnostics, &
      SetProblem

      private :: &
        SetFluid

        private :: &
          SetFluidKernel

!       private :: &
!         SetFluidKernel, &
!         SetReferenceKernel

!   interface 
!     subroutine DerivativeInterface ( Parameters, X, Y, dYdX )
!       use GenASiS
!       class ( * ), intent ( in ) :: &
!         Parameters
!       real ( KDR ), intent ( in ) :: &
!         X
!       real ( KDR ), dimension ( : ), intent ( in ) :: &
!         Y
!       real ( KDR ), dimension ( : ), intent ( out ) :: &
!         dYdX
!     end subroutine DerivativeInterface
!   end interface

contains


  subroutine Initialize_OS ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    character ( * ), intent ( in ) :: &
      Name

    if ( OS % Type == '' ) &
      OS % Type = 'an OppenheimerSnyder'

    call InitializeFluidCentralCore ( OS, Name )
    call InitializeDiagnostics ( OS )
    call SetProblem ( OS )

!    !-- ODE Solver
!    OS % DerivativeEvaluator => DerivativeFunction
    
!    call OS % Parameters % Initialize ( 2 )
!    OS % Parameters % Value ( 1 ) = 1.0_KDR
!    OS % Parameters % Value ( 2 ) &
!           = OS % DensityInitial * 4.0_KDR / 3.0_KDR * CONSTANT % PI

!    allocate ( OS % ODEF )
    
!    call OS % ODEF % Initialize &
!                       ( OS % Parameters, &
!                         OS % DerivativeEvaluator )
   
  end subroutine Initialize_OS


  subroutine ComputeError ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

!     integer ( KDI ) :: &
!       nM, &
!       iV, &
!       iE, &
!       iM !-- iMessage
!     real ( KDR ) :: &
!       L1
!     class ( Fluid_D_Form ), pointer :: &
!       F, &
!       F_R, &  !-- F_Reference
!       F_D     !-- F_Difference         
!     type ( CollectiveOperation_R_Form ) :: &
!       CO
    
!     select type ( FCC => OS % Integrator )
!     type is ( FluidCentralCoreForm )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_D ( )

!     F_R => OppenheimerSnyderCollapse % Reference % Fluid_D ( )
!     call SetReferenceKernel ( OppenheimerSnyderCollapse, F_R, FCC % Time )

!     F_D => OppenheimerSnyderCollapse % Difference % Fluid_D ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )


!     select type ( PS => OS % Integrator % PositionSpace )
!     class is ( Atlas_SC_Form )
    
!     select type ( C => PS % Chart ) 
!     class is ( Chart_SLD_Form )
    
!     associate ( D   => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
!                 R   => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ) )


!     call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )

!     CO % Outgoing % Value ( 1 ) &
!            = sum ( abs ( R ), mask = C % IsProperCell )
!     CO % Outgoing % Value ( 2 )&
!            = sum ( abs ( D ), mask = C % IsProperCell )

!     call CO % Reduce ( REDUCTION % SUM )

!     end associate !-- D, etc. 
    
!     !-- L1 error

!     associate ( IN   => CO % Incoming % Value )

!     L1 = IN ( 2 ) / IN ( 1 )
!     call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
!                 nTrailingLinesOption = 2 )

!     end associate !-- IN 

!     end select !-- C
!     end select !-- PS

!     end select !-- FA
!     end select !-- FCC
!     nullify ( F, F_R, F_D )

  end subroutine


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    if ( allocated ( OS % Difference ) ) &
      deallocate ( OS % Difference )
    if ( allocated ( OS % Reference ) ) &
      deallocate ( OS % Reference )

  end subroutine Finalize


!   subroutine SetReference ( FCC )
!     class ( IntegratorTemplate ), intent ( inout ) :: &
!       FCC 

!     class ( Fluid_D_Form ), pointer :: &
!       F, &
!       F_R, &  !-- F_Reference
!       F_D     !-- F_Difference

!     select type ( FCC )
!     class is ( FluidCentralCoreForm )
!     select type ( FA => FCC % Current_ASC )
!     class is ( Fluid_ASC_Form )
!     F => FA % Fluid_D ( )

!     F_R => OppenheimerSnyderCollapse % Reference % Fluid_D ( )
!     call SetReferenceKernel ( OppenheimerSnyderCollapse, F_R, FCC % Time )

!     F_D => OppenheimerSnyderCollapse % Difference % Fluid_D ( )
!     call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

!     F_D % Value = abs ( F_D % Value )

!     end select !-- FA
!     end select !-- FCC
!     nullify ( F, F_R, F_D )

!   end subroutine SetReference


  subroutine InitializeFluidCentralCore ( OS, Name )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    character ( * ), intent ( in )  :: &
      Name

    call OS % Initialize &
           ( FluidType = 'DUST', &
             GeometryType = 'NEWTONIAN', &
             Name = Name, &
             DimensionlessOption = .true., &
             GravityFactorOption = 0.01_KDR, &
             LimiterParameterOption = 1.0_KDR )
!     OS % SetReference => SetReference

  end subroutine InitializeFluidCentralCore


  subroutine InitializeDiagnostics ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    select type ( PS => OS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )

    allocate ( OS % Reference )
    allocate ( OS % Difference )
    call OS % Reference % Initialize &
           ( PS, 'DUST', OS % Units, NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call OS % Difference % Initialize &
           ( PS, 'DUST', OS % Units, NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    end select !-- PS

  end subroutine InitializeDiagnostics


  subroutine SetProblem ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    real ( KDR ) :: &
      Mass, &
      DensityFactor, &
      RadiusFactor, &
      Eta
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_D_Form ), pointer :: &
      F

    select type ( I => OS % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    associate &
      ( R_Max => PS % Chart % MaxCoordinate ( 1 ), &
          R_0 => OS % RadiusInitial, &
          D_0 => OS % DensityInitial, &
          Tau => OS % TimeScale, &
            M => Mass, &
           DF => DensityFactor, &
           RF => RadiusFactor, &
           PI => CONSTANT % PI )

      M  =  1.0_KDR
    D_0  =  1.0e-3_KDR
     DF  =  1.0e2_KDR
    call PROGRAM_HEADER % GetParameter (   M, 'Mass' )
    call PROGRAM_HEADER % GetParameter ( D_0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter (  DF, 'DensityFactor' )

               Tau  =  sqrt ( 3.0 / ( 8.0 * PI * D_0 ) )
               R_0  =  ( 3.0 * M / ( 4.0 * PI * D_0 ) ) ** ( 1.0_KDR / 3.0_KDR )
                RF  =  DF ** ( - 1.0_KDR / 3.0_KDR ) 
               Eta  =  acos ( 2.0 * RF  -  1.0 )
    I % FinishTime  =  0.5 * Tau * ( Eta  +  sin ( Eta ) )

    call Show ( 'OppenheimerSnyder parameters' )
    call Show ( M, 'Mass' )
    call Show ( D_0, 'DensityInitial' )
    call Show ( R_0, 'RadiusInitial' )
    call Show ( DF, 'DensityFactor' )
    call Show ( RF, 'RadiusFactor' )
    call Show ( PI / 2  *  Tau, 'CollapseTime' )
    call Show ( I % FinishTime, 'Reset FinishTime' )

    if ( R_0 > R_Max  ) then
      call Show ( 'RadiusInitial too large', CONSOLE % ERROR )
      call Show ( R_0, 'RadiusInitial', CONSOLE % ERROR )
      call Show ( PS % Chart % MaxCoordinate ( 1 ), 'RadiusMax', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    G => PS % Geometry ( )
    F => FA % Fluid_D ( )
    call SetFluid ( OS, F, G, Time = 0.0_KDR )

    end associate !-- D0, etc.
    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine SetFluid ( OS, F, G, Time )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    real ( KDR ), intent ( in ) :: &
      Time
    
    real ( KDR ) :: &
      Radius, &
      Density

    if ( Time == 0.0_KDR ) then
      Radius   =  OS % RadiusInitial
      Density  =  OS % DensityInitial
    end if

    call SetFluidKernel &
           (    R = G % Value ( :, G % CENTER_U ( 1 ) ), &
             dR_L = G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
             dR_R = G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ), &
             Density = Density, &
             RadiusDensity = Radius, &
              N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
             VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
             VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    call F % ComputeFromPrimitive ( G )

  end subroutine SetFluid


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


!   subroutine SetReferenceKernel ( OS, F, Time )
!     class ( OppenheimerSnyderForm ), intent ( inout ) :: &
!       OS
!     class ( Fluid_D_Form ), intent ( inout ) :: &
!       F
!     real ( KDR ), intent ( in ) :: &
!       Time

!     integer ( KDI ) :: &
!       iC
!     real ( KDR ) :: &
!       X_1, &
!       X_2, &
!       H, &
!       R_New
!     real ( KDR ), dimension ( 2 ) :: &
!       Y_Start
!     class ( GeometryFlatForm ), pointer :: &
!       G

!     select type ( PS => OS % Integrator % PositionSpace )
!     class is ( Atlas_SC_Form )
!     G => PS % Geometry ( )
    
!     select type ( C => PS % Chart ) 
!     class is ( Chart_SLD_Form )

!     associate &
!       ( Rho   => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
!         R     => G % Value ( :, G % CENTER_U ( 1 ) ), &
!         Rho_0 => OS % DensityInitial, &
!         R0    => OS % RadiusInitial )

!     X_1 = 0.0_KDR
!     X_2 = Time

!     H = X_2 / 2

!     Y_Start ( 1 ) = 1.0_KDR
!     Y_Start ( 2 ) = 0.0_KDR

!     if ( Time > 0.0_KDR ) then 
!       !-- Solve system of ODEs
!       call OS % ODEF % IntegrateODE &
!                          ( Y_Start, X_1, X_2, H, &
!                            OS % ODEF % RequestedAccuracy, &
!                            RK4Option = .false. ) 
!     end if
 
!     R_New = R0 * Y_Start ( 1 ) 
    
!     Rho = 0.0_KDR
 
!     where ( R <= R_New ) 

!       Rho = Rho_0 / ( Y_Start ( 1 ) ** 3 )

!     end where

!     end associate !-- Rho, etc.
!     end select !-- C
!     end select !-- PS
!     nullify ( G )

!   end subroutine SetReferenceKernel


!   subroutine DerivativeFunction ( Parameters, X, Y, dYdX )
!     class ( * ), intent ( in ) :: &
!       Parameters
!     real ( KDR ), intent ( in ) :: &
!       X
!     real ( KDR ), dimension ( : ), intent ( in ) :: &
!       Y
!     real ( KDR ), dimension ( : ), intent ( out ) :: &
!       dYdX

!     real ( KDR ) :: &
!       SqrtTiny

!     SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

!     select type ( P => Parameters ) 
!     type is ( Real_1D_Form )
!     associate &
!       ( G => P % Value ( 1 ), &
!         M => P % Value ( 2 ) )

!     dYdX ( 1 ) = Y ( 2 )
!     dYdX ( 2 ) = - G * M / max ( Y ( 1 ) ** 2, SqrtTiny )

!     end associate
!     end select !-- P
  
!   end subroutine DerivativeFunction


end module OppenheimerSnyder_Form
