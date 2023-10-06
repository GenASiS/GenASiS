module RayleighTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: RayleighTaylorForm
    real ( KDR ), private :: &
      Acceleration, &
      DensityAbove, DensityBelow, &
      PressureBase, &   
      AdiabaticIndex
  contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type RayleighTaylorForm

    private :: &
      InitializeUniverse, &
      SetInitial

      private :: &
        SetFluid

        private :: &
          SetFluidKernel_2D, &
          SetFluidKernel_3D
    

contains


  subroutine Initialize_H ( U, Name, CommunicatorOption )

    class ( RayleighTaylorForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ) :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption

    if ( U % Type  ==  '' ) &
      U % Type  =  'a RayleighTaylor'

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( RayleighTaylorForm ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % Acceleration,   'Acceleration' )
    call Show ( U % DensityAbove,   'DensityAbove' )
    call Show ( U % DensityBelow,   'DensityBelow' )
    call Show ( U % PressureBase,   'PressureBase' )
    call Show ( U % AdiabaticIndex, 'AdiabaticIndex' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ), target :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate

    select case ( trim ( PROGRAM_HEADER % Dimensionality ) )
    case ( '2D' )
      MinCoordinate = [ -0.25_KDR, -0.75_KDR, 0.0_KDR ]
      MaxCoordinate = [ +0.25_KDR, +0.75_KDR, 0.0_KDR ]
      nCells = [ 64, 192, 1 ]
    case ( '3D' )
      MinCoordinate = [ -0.25_KDR, -0.25_KDR, -0.75_KDR ]
      MaxCoordinate = [ +0.25_KDR, +0.25_KDR, +0.75_KDR ]
      nCells = [ 64, 64, 192 ]
    end select

    RT % Acceleration  =  0.1_KDR
    call PROGRAM_HEADER % GetParameter &
           ( RT % Acceleration, 'Acceleration' )

    call RT % Initialize &
           ( FluidType = 'IDEAL', &
             GravitationType = 'NEWTON_UA', &
             Name = Name, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             FinishTimeOption = 8.5_KDR, &
             UniformAccelerationOption = RT % Acceleration, &
             nCellsOption = nCells )

    select type ( I  =>  RT % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    select case ( trim ( PROGRAM_HEADER % Dimensionality ) )
    case ( '2D' )
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 1 )
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    case ( '3D' )
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 1 )
      call F % SetBoundaryConditionsFace &
             ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 2 )
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 3 )
    end select !-- Dimensionality
    end associate !-- F
    end select !-- I
             
    RT % Integrator % SetInitial  =>  SetInitial
    RT % Integrator % System      =>  RT

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    select type ( RT  =>  I % System )
      class is ( RayleighTaylorForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )

    RT % DensityAbove    =  2.0_KDR
    RT % DensityBelow    =  1.0_KDR
    RT % PressureBase    =  2.5_KDR
    RT % AdiabaticIndex  =  1.4_KDR

    call PROGRAM_HEADER % GetParameter ( RT % DensityAbove, 'DensityAbove' )
    call PROGRAM_HEADER % GetParameter ( RT % DensityBelow, 'DensityBelow' )
    call PROGRAM_HEADER % GetParameter ( RT % AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter ( RT % Acceleration, 'Acceleration' )    

    call SetFluid ( RT, F )

    end select !-- F
    end select !-- I
    end select !-- RT

  end subroutine SetInitial


  subroutine SetFluid ( RT, F )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    call F % SetAdiabaticIndex ( RT % AdiabaticIndex )

    associate &
      ( G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    select case ( trim ( PROGRAM_HEADER % Dimensionality ) )
    case ( '2D' )

      call SetFluidKernel_2D &
             (       X = GV ( :, G % CENTER_U_1 ), &
                     Y = GV ( :, G % CENTER_U_2 ), &
               N_Above = RT % DensityAbove, &
               N_Below = RT % DensityBelow, &
                     A = RT % Acceleration, &
                 Gamma =  F % AdiabaticIndex, &
                   P_0 = RT % PressureBase, &
                     N = FV ( :, F % BARYON_DENSITY_C ), &
                     E = FV ( :, F % ENERGY_DENSITY_C ), &
                    VY = FV ( :, F % VELOCITY_U_2 ) )

    case ( '3D' )

      call SetFluidKernel_3D &
             (       X = GV ( :, G % CENTER_U_1 ), &
                     Y = GV ( :, G % CENTER_U_2 ), &
                     Z = GV ( :, G % CENTER_U_3 ), &
               N_Above = RT % DensityAbove, &
               N_Below = RT % DensityBelow, &
                     A = RT % Acceleration, &
                 Gamma =  F % AdiabaticIndex, &
                   P_0 = RT % PressureBase, &
                     N =  FV ( :, F % BARYON_DENSITY_C ), &
                     E =  FV ( :, F % ENERGY_DENSITY_C ), &
                    VZ =  FV ( :, F % VELOCITY_U_3 ) )

    end select !-- nDimensions

    end associate !-- FV, etc.
    end associate !-- G

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
