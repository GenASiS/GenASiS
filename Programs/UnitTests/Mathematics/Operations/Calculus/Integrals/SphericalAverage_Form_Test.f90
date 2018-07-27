module SphericalAverage_Form_Test__Form

  use Basics
  use Manifolds
  use SphericalAverage_Form

  implicit none
  private

  character ( LDF ), public, parameter :: &
    ProgramName = 'SphericalAverage_Form_Test'

  type, public :: SphericalAverage_Form_Test_Form
    type ( GridImageStreamForm ), allocatable :: &
      GridImageStream, &
      GridImageStream_SA
    type ( Atlas_SC_CC_Form ), allocatable :: &
      Atlas
    type ( Atlas_SC_Form ), allocatable :: &
      Atlas_SA
    type ( Storage_ASC_Form ), allocatable :: &
      Integrand_ASC, &
      Average_ASC
    type ( SphericalAverageForm ), allocatable :: &
      SphericalAverage
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type SphericalAverage_Form_Test_Form

contains


  subroutine Initialize ( SAFT, Name )

    class ( SphericalAverage_Form_Test_Form ), intent ( inout ), target :: &
      SAFT
    character ( * ), intent ( in ) :: &
      Name

!     real ( KDR ), dimension ( 1 ) :: &
!       Integral
!     type ( Real_1D_Form ), dimension ( 1 ) :: &
!       Integrand
    real ( KDR ) :: &
      A1, A2, A3
    real ( KDR ), dimension ( : ), allocatable :: &
      X, Y, Z
    type ( Real_1D_Form ), dimension ( 1 ) :: &
      Edge
    class ( StorageForm ), pointer :: &
      I
    class ( GeometryFlatForm ), pointer :: &
      G

    !-- Stream

    allocate &
      ( SAFT % GridImageStream, &
        SAFT % GridImageStream_SA )
    associate &
      ( GIS    => SAFT % GridImageStream, &
        GIS_SA => SAFT % GridImageStream_SA )
    call GIS % Initialize &
           ( Name, CommunicatorOption = PROGRAM_HEADER % Communicator )
    call GIS_SA % Initialize &
           ( trim ( Name ) // '_SA' )
    
    !-- Atlases, etc.

    allocate &
      ( SAFT % Atlas, &
        SAFT % Atlas_SA )
    associate &
      ( A    => SAFT % Atlas, &
        A_SA => SAFT % Atlas_SA )

    call A % Initialize ( Name, PROGRAM_HEADER % Communicator )
    call A % CreateChart_CC ( )  
    call A % SetGeometry ( )
    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )

    call Edge ( 1 ) % Initialize &
           ( C % Edge ( 1 ) % Value ( 1 : C % nCells ( 1 ) + 1 ) )

    call A_SA % Initialize ( 'SphericalAverage', nDimensionsOption = 1 )
    call A_SA % CreateChart &
           ( CoordinateSystemOption = 'SPHERICAL', &
             nCellsOption           = C % nCells ( 1 : 1 ), &
             nGhostLayersOption     = [ 0 ] )
    call A_SA % SetGeometry ( EdgeOption = Edge )

    !-- Integrand
    
    allocate &
      ( SAFT % Integrand_ASC, &
        SAFT % Average_ASC )
    associate &
      ( IA => SAFT % Integrand_ASC, &
        AA => SAFT % Average_ASC )
    call IA % Initialize &
           ( A, 'Integrand', nFields = 3, &
             VariableOption = [ 'Sphere   ', 'Spheroid ', 'Ellipsoid' ], &
             WriteOption = .true. )
    call AA % Initialize &
           ( A_SA, 'Average', nFields = 3, &
             VariableOption = [ 'Sphere   ', 'Spheroid ', 'Ellipsoid' ], &
             WriteOption = .true. )

    A1  =  0.6_KDR * C % MaxCoordinate ( 1 )
    A2  =  0.4_KDR * C % MaxCoordinate ( 1 )
    A3  =  0.2_KDR * C % MaxCoordinate ( 1 )

    I => IA % Storage ( )

    associate &
      ( R     => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Phi   => G % Value ( :, G % CENTER_U ( 3 ) ), &
        Sphere    => I % Value ( :, 1 ), &
        Spheroid  => I % Value ( :, 2 ), &
        Ellipsoid => I % Value ( :, 3 ) )

    allocate ( X ( G % nValues ), Y ( G % nValues ), Z ( G % nValues ) )

    select case ( A % nDimensions )
    case ( 1 )
      X  =  R
      Y  =  0.0_KDR
      Z  =  0.0_KDR
    case ( 2 )
      X  =  R * Sin ( Theta )
      Y  =  0.0_KDR
      Z  =  R * Cos ( Theta )
    case ( 3 )
      X  =  R * Sin ( Theta ) * Cos ( Phi )
      Y  =  R * Sin ( Theta ) * Sin ( Phi )
      Z  =  R * Cos ( Theta )
    end select !-- nDimensions

    where ( ( X ** 2  +  Y ** 2  +  Z ** 2 ) / A1 ** 2  <  1.0_KDR )
      Sphere  =  1.0_KDR
    elsewhere
      Sphere  =  0.0_KDR
    end where

    where ( ( X ** 2  +  Y ** 2 ) / A1 ** 2  &
            +  Z ** 2 / A2 ** 2  <  1.0_KDR )
      Spheroid  =  1.0_KDR
    elsewhere
      Spheroid  =  0.0_KDR
    end where

    where ( X ** 2 / A1 ** 2  +  Y ** 2 / A2 ** 2 &
            +  Z ** 2 / A3 ** 2  <  1.0_KDR )
      Ellipsoid  =  1.0_KDR
    elsewhere
      Ellipsoid  =  0.0_KDR
    end where

!     !-- SphericalAverage

!     allocate ( SAFT % SphericalAverage )
!     associate ( VI => SAFT % SphericalAverage )

!     call VI % Compute ( C, Integrand, Integral )
!     call Show ( Integral ( 1 ), '*** SphericalAverage', nLeadingLinesOption = 1 )
!     call Show ( 4.0_KDR / 3.0_KDR  *  CONSTANT % PI  &
!                 *  C % MaxCoordinate ( 1 ) ** 3, &
!                 '*** Expected', nTrailingLinesOption = 1 )

    !-- Write

    call A % OpenStream ( GIS, '1', iStream = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call A % Write ( iStream = 1 )
    call GIS % Close ( )

    !-- Cleanup

!     end associate !-- VI
    end associate !-- R, etc.
    end associate !-- IA, etc.
    end select !-- C
    end associate !-- A, etc.
    end associate !-- GIS, etc.
    nullify ( I, G )
    
  end subroutine Initialize


  subroutine Finalize ( SAFT )

    type ( SphericalAverage_Form_Test_Form ), intent ( inout ) :: &
      SAFT

!     deallocate ( SAFT % SphericalAverage )
     deallocate &
       ( SAFT % Average_ASC, &
         SAFT % Integrand_ASC )
     deallocate &
       ( SAFT % Atlas_SA, &
         SAFT % Atlas )
     deallocate &
       ( SAFT % GridImageStream_SA, &
         SAFT % GridImageStream )

  end subroutine Finalize


end module SphericalAverage_Form_Test__Form



program SphericalAverage_Form_Test

  use Basics
  use SphericalAverage_Form_Test__Form

  implicit none
  
  type ( SphericalAverage_Form_Test_Form ), allocatable :: &
    SAFT
    
  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( ProgramName )
    
  allocate ( SAFT )
  call SAFT % Initialize ( PROGRAM_HEADER % Name )
  deallocate ( SAFT )

  deallocate ( PROGRAM_HEADER )

end program SphericalAverage_Form_Test
