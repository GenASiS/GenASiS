module OppenheimerSnyder_Form

  !-- For example, Misner, Thorne, Wheeler p. 663, Eqs. (25.28)-(25.29)

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidCentralCoreForm ) :: OppenheimerSnyderForm
    real ( KDR ) :: &
      DensityInitial, &
      RadiusInitial, &
      TimeScale
    type ( RootFinderForm ), allocatable :: &
      RootFinder
    type ( Fluid_ASC_Form ), allocatable :: &
      Reference, &
      Difference
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
      SetReference, &
      InitializeFluidCentralCore, &
      InitializeDiagnostics, &
      SetProblem

      private :: &
        EvaluateZeroEta, &
        SetFluid

        private :: &
          SetFluidKernel

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

  end subroutine Initialize_OS


  subroutine ComputeError ( OS )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    real ( KDR ) :: &
      L1
    class ( Fluid_D_Form ), pointer :: &
      F_D, &
      F_R
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    select type ( PS => OS % Integrator % PositionSpace )
    class is ( Atlas_SC_Form )
    
    select type ( PSC => PS % Chart ) 
    class is ( Chart_SLD_Form )
    
    F_D => OS % Difference % Fluid_D ( )
    F_R => OS % Reference  % Fluid_D ( )

    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )

    associate &
      ( D => F_D % Value ( :, F_D % COMOVING_BARYON_DENSITY ), &
        R => F_R % Value ( :, F_R % COMOVING_BARYON_DENSITY ), &
        Norm_D => CO % Incoming % Value ( 1 ), &
        Norm_R => CO % Incoming % Value ( 2 ) )

    CO % Outgoing % Value ( 1 ) &
      = sum ( abs ( D ), mask = PSC % IsProperCell )
    CO % Outgoing % Value ( 2 ) &
      = sum ( abs ( R ), mask = PSC % IsProperCell )

    call CO % Reduce ( REDUCTION % SUM )

    L1 = Norm_D / Norm_R
    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end associate !-- D, etc.
    end select !-- PSC
    end select !-- PS
    nullify ( F_D, F_R )

  end subroutine ComputeError


  impure elemental subroutine Finalize ( OS )
    
    type ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS

    if ( allocated ( OS % Difference ) ) &
      deallocate ( OS % Difference )
    if ( allocated ( OS % Reference ) ) &
      deallocate ( OS % Reference )
    if ( allocated ( OS % RootFinder ) ) &
      deallocate ( OS % RootFinder )

  end subroutine Finalize


  subroutine SetReference ( I )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      I

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_D_Form ), pointer :: &
      F, &
      F_R, &  !-- F_Reference
      F_D     !-- F_Difference

    select type ( I )
    class is ( Integrator_C_PS_Form )

    select type ( OS => I % Universe )
    class is ( OppenheimerSnyderForm )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )
    F => FA % Fluid_D ( )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    F_R => OS % Reference % Fluid_D ( )
    call SetFluid ( OS, F_R, G )

    F_D => OS % Difference % Fluid_D ( )
    call MultiplyAdd ( F % Value, F_R % Value, -1.0_KDR, F_D % Value )

    end select !-- PS
    end select !-- FA
    end select !-- OS
    end select !-- I
    nullify ( G, F, F_R, F_D )

  end subroutine SetReference


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
     OS % Integrator % SetReference => SetReference

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
      Pi, &
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
           RF => RadiusFactor )

     Pi  =  CONSTANT % PI
      M  =  1.0_KDR
    D_0  =  1.0e-3_KDR
     DF  =  1.0e2_KDR
    call PROGRAM_HEADER % GetParameter (   M, 'Mass' )
    call PROGRAM_HEADER % GetParameter ( D_0, 'DensityInitial' )
    call PROGRAM_HEADER % GetParameter (  DF, 'DensityFactor' )

               Tau  =  sqrt ( 3.0 / ( 8.0 * Pi * D_0 ) )
               R_0  =  ( 3.0 * M / ( 4.0 * Pi * D_0 ) ) ** ( 1.0_KDR / 3.0_KDR )
                RF  =  DF ** ( - 1.0_KDR / 3.0_KDR ) 
               Eta  =  acos ( 2.0 * RF  -  1.0 )
    I % FinishTime  =  0.5 * Tau * ( Eta  +  sin ( Eta ) )

    call Show ( 'OppenheimerSnyder parameters' )
    call Show ( M, 'Mass' )
    call Show ( D_0, 'DensityInitial' )
    call Show ( R_0, 'RadiusInitial' )
    call Show ( DF, 'DensityFactor' )
    call Show ( RF, 'RadiusFactor' )
    call Show ( Pi / 2  *  Tau, 'CollapseTime' )
    call Show ( I % FinishTime, 'Reset FinishTime' )

    if ( R_0 > R_Max  ) then
      call Show ( 'RadiusInitial too large', CONSOLE % ERROR )
      call Show ( R_0, 'RadiusInitial', CONSOLE % ERROR )
      call Show ( PS % Chart % MaxCoordinate ( 1 ), 'RadiusMax', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    allocate ( OS % RootFinder )
    associate ( RF => OS % RootFinder )
    call RF % Initialize ( OS )
    RF % EvaluateZero  =>  EvaluateZeroEta

    G => PS % Geometry ( )
    F => FA % Fluid_D ( )
    call SetFluid ( OS, F, G )

    end associate !-- RF
    end associate !-- D0, etc.
    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine EvaluateZeroEta ( OS, Eta, Zero )

    class ( * ), intent ( in ) :: &
      OS
    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ), intent ( out ) :: &
      Zero

    select type ( OS )
    class is ( OppenheimerSnyderForm )

    associate &
      ( Tau  => OS % TimeScale, &
        Time => OS % Integrator % Time )

    Zero  =  2.0 * Time / Tau  -  ( Eta  +  sin ( Eta ) )

    end associate !-- Tau, etc.
    end select !-- OS
    
  end subroutine EvaluateZeroEta


  subroutine SetFluid ( OS, F, G )

    class ( OppenheimerSnyderForm ), intent ( inout ) :: &
      OS
    class ( Fluid_D_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G
    
    real ( KDR ) :: &
      Pi, &
      Eta, &
      Radius, &
      Density

    Pi  =  CONSTANT % PI

    associate &
      ( RF => OS % RootFinder, &
        D_0 => OS % DensityInitial, &
        R_0 => OS % RadiusInitial )

    call RF % Solve ( [ 0.0_KDR, Pi ], Eta )

    Radius   =  0.5 * R_0 * ( 1 + cos ( Eta ) )
    Density  =  D_0 * ( R_0 / Radius ) ** 3

    end associate !-- RF, etc.

  
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


end module OppenheimerSnyder_Form
