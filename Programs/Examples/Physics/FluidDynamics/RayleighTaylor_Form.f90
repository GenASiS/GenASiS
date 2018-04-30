module RayleighTaylor_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RayleighTaylorForm
    real ( KDR ), private :: &
      DensityAbove, DensityBelow, &
      PressureBase, &   
      AdiabaticIndex, &
      Acceleration
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RayleighTaylorForm

    private :: &
      SetFluid

contains


  subroutine Initialize ( RT, Name )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace

    if ( RT % Type == '' ) &
      RT % Type = 'a RayleighTaylor'

    call RT % InitializeTemplate ( Name )

    RT % DensityAbove   = 2.0_KDR
    RT % DensityBelow   = 1.0_KDR
    RT % PressureBase   = 2.5_KDR
    RT % AdiabaticIndex = 1.4_KDR
    RT % Acceleration   = 0.1_KDR

    call PROGRAM_HEADER % GetParameter &
           ( RT % DensityAbove, 'DensityAbove' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % DensityBelow, 'DensityBelow' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % AdiabaticIndex, 'AdiabaticIndex' )
    call PROGRAM_HEADER % GetParameter &
           ( RT % Acceleration, 'Acceleration' )    
    call PROGRAM_HEADER % GetParameter &
           ( RT % Acceleration, 'Acceleration' )

    !-- Integrator

    allocate ( FluidBoxForm :: RT % Integrator )
    select type ( FB => RT % Integrator )
    type is ( FluidBoxForm )

    associate ( BCF => BoundaryConditionsFace )
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

    call FB % Initialize &
           ( Name, FluidType = 'IDEAL', GeometryType = 'NEWTONIAN', &
             BoundaryConditionsFaceOption = BCF, &
             GravitySolverTypeOption = 'UNIFORM', &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             FinishTimeOption = 8.5_KDR, &
             UniformAccelerationOption = RT % Acceleration, &
             nCellsOption = nCells )

    end associate !-- BCF
    end select !-- FB


    !-- Initial conditions

    call SetFluid ( RT )


  end subroutine Initialize


  subroutine Finalize ( RT )

    type ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    call RT % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( RT )

    class ( RayleighTaylorForm ), intent ( inout ) :: &
      RT

    integer ( KDI ) :: &
      iB  !-- iBoundary
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( FB => RT % Integrator )
    class is ( FluidBoxForm )

    select type ( FA => FB % Current_ASC )
    class is ( Fluid_ASC_Form )

    !-- Initial conditions
       
    F => FA % Fluid_P_I ( )

    call F % SetAdiabaticIndex ( RT % AdiabaticIndex )

    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      ( X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        N  => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        E  => F % Value ( :, F % INTERNAL_ENERGY ), &
        VY => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        Pi => CONSTANT % PI )

    where ( Y > 0.0_KDR )     
      N  = RT % DensityAbove
    elsewhere
      N  = RT % DensityBelow
    end where

    E  =  ( RT % PressureBase  -  N  *  RT % Acceleration  *  Y ) &
          / ( RT % AdiabaticIndex  -  1.0_KDR )

    VY  =  ( 0.01_KDR / 4.0_KDR ) &
           *  ( 1.0_KDR  +  cos ( 4.0_KDR * Pi * X ) ) &
           *  ( 1.0_KDR  +  cos ( 3.0_KDR * Pi * Y ) )
    
    call F % ComputeFromPrimitive ( G )

    end associate !-- X, etc.
    end select !-- PS
    end select !-- FA
    end select !-- FB
    nullify ( F, G )
    
  end subroutine SetFluid

  
end module RayleighTaylor_Form
