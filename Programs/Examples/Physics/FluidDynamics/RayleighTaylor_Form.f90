module RayleighTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidBoxForm ) :: RayleighTaylorForm
    real ( KDR ), private :: &
      Acceleration, &
      DensityAbove, DensityBelow, &
      PressureBase, &   
      AdiabaticIndex
  contains
    procedure, private, pass :: &
      Initialize_RT
    generic, public :: &
      Initialize => Initialize_RT
    final :: &
      Finalize
  end type RayleighTaylorForm

    private :: &
      InitializeFluidBox, &
      SetProblem

      private :: &
        SetFluid

        private :: &
          SetFluidKernel_2D, &
          SetFluidKernel_3D
    
contains


  subroutine Initialize_RT ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    if ( RT % Type == '' ) &
      RT % Type = 'a RayleighTaylor'

    call InitializeFluidBox ( RT, Name )
    call SetProblem ( RT )

  end subroutine Initialize_RT


  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

  end subroutine Finalize


  subroutine InitializeFluidBox ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate

    allocate ( RT % BoundaryConditionsFace ( 3 ) )
    associate ( BCF => RT % BoundaryConditionsFace )
    select case ( trim ( PROGRAM_HEADER % Dimensionality ) )
    case ( '2D' )
      MinCoordinate = [ -0.25_KDR, -0.75_KDR, 0.0_KDR ]
      MaxCoordinate = [ +0.25_KDR, +0.75_KDR, 0.0_KDR ]
      nCells = [ 64, 192, 1 ]
      call BCF ( 1 ) % Initialize ( [ 'PERIODIC', 'PERIODIC' ] )
      call BCF ( 2 ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )
    case ( '3D' )
      MinCoordinate = [ -0.25_KDR, -0.25_KDR, -0.75_KDR ]
      MaxCoordinate = [ +0.25_KDR, +0.25_KDR, +0.75_KDR ]
      nCells = [ 64, 64, 192 ]
      call BCF ( 1 ) % Initialize ( [ 'PERIODIC', 'PERIODIC' ] )
      call BCF ( 2 ) % Initialize ( [ 'PERIODIC', 'PERIODIC' ] )
      call BCF ( 3 ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )
    end select
    end associate !-- BCF

    RT % Acceleration  =  0.1_KDR
    call PROGRAM_HEADER % GetParameter &
           ( RT % Acceleration, 'Acceleration' )

    call RT % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', Name = Name, &
             GravitySolverTypeOption = 'UNIFORM', &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             FinishTimeOption = 8.5_KDR, &
             UniformAccelerationOption = RT % Acceleration, &
             nCellsOption = nCells )

  end subroutine InitializeFluidBox


  subroutine SetProblem ( RT )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => RT % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    RT % DensityAbove   = 2.0_KDR
    RT % DensityBelow   = 1.0_KDR
    RT % PressureBase   = 2.5_KDR
    RT % AdiabaticIndex = 1.4_KDR

    call PROGRAM_HEADER % GetParameter &
           ( RT % DensityAbove, 'DensityAbove' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % DensityBelow, 'DensityBelow' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % Acceleration, 'Acceleration' )    

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( RT, F, G )

    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine SetFluid ( RT, F, G )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    select type ( I => RT % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

    call F % SetAdiabaticIndex ( RT % AdiabaticIndex )

    select case ( PS % nDimensions )
    case ( 2 )

      call SetFluidKernel_2D &
             (       X = G % Value ( :, G % CENTER_U ( 1 ) ), &
                     Y = G % Value ( :, G % CENTER_U ( 2 ) ), &
               N_Above = RT % DensityAbove, &
               N_Below = RT % DensityBelow, &
                     A = RT % Acceleration, &
                 Gamma =  F % AdiabaticIndex, &
                   P_0 = RT % PressureBase, &
                     N =  F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                     E =  F % Value ( :, F % INTERNAL_ENERGY ), &
                    VY = F % Value ( :, F % VELOCITY_U ( 2 ) ) )

    case ( 3 )

      call SetFluidKernel_3D &
             (       X = G % Value ( :, G % CENTER_U ( 1 ) ), &
                     Y = G % Value ( :, G % CENTER_U ( 2 ) ), &
                     Z = G % Value ( :, G % CENTER_U ( 3 ) ), &
               N_Above = RT % DensityAbove, &
               N_Below = RT % DensityBelow, &
                     A = RT % Acceleration, &
                 Gamma =  F % AdiabaticIndex, &
                   P_0 = RT % PressureBase, &
                     N =  F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                     E =  F % Value ( :, F % INTERNAL_ENERGY ), &
                    VZ =  F % Value ( :, F % VELOCITY_U ( 3 ) ) )

    end select !-- nDimensions

    call F % ComputeFromPrimitive ( G )

    end select !-- PS
    end select !-- I

  end subroutine SetFluid


  subroutine SetFluidKernel_2D &
               ( X, Y, N_Above, N_Below, A, Gamma, P_0, N, E, VY )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y
    real ( KDR ) :: &
      N_Above, N_Below, &
      A, &
      Gamma, &
      P_0
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      E, &
      VY

    real ( KDR ) :: &
      Pi

    Pi  =  CONSTANT % PI

    where ( Y > 0.0_KDR )     
      N  = N_Above
    elsewhere
      N  = N_Below
    end where

    E  =  ( P_0  -  N * A * Y )  /  ( Gamma  -  1.0_KDR )

    VY  =  ( 0.01_KDR / 4.0_KDR ) &
           *  ( 1.0_KDR  +  cos ( 4.0_KDR * Pi * X ) ) &
           *  ( 1.0_KDR  +  cos ( 3.0_KDR * Pi * Y ) )
    
  end subroutine SetFluidKernel_2D


  subroutine SetFluidKernel_3D &
               ( X, Y, Z, N_Above, N_Below, A, Gamma, P_0, N, E, VZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ) :: &
      N_Above, N_Below, &
      A, &
      Gamma, &
      P_0
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      E, &
      VZ

    real ( KDR ) :: &
      Pi

    Pi  =  CONSTANT % PI

    where ( Z > 0.0_KDR )     
      N  = N_Above
    elsewhere
      N  = N_Below
    end where
    
    E  =  ( P_0  -  N * A * Z )  /  ( Gamma  -  1.0_KDR )

    VZ  =  ( 0.01_KDR / 4.0_KDR ) &
           *  ( 1.0_KDR  +  cos ( 4.0_KDR * Pi &
                                  * sqrt ( X ** 2  +  Y ** 2 ) ) ) &
           *  ( 1.0_KDR  +  cos ( 3.0_KDR * Pi * Z ) )

  end subroutine SetFluidKernel_3D


end module RayleighTaylor_Form
