program SurfaceIntegral_Form_Test

  use Basics
  use Manifolds
  use Geometries
  use Integrals

  implicit none

  integer ( KDI ) :: &
    iD  !-- iDimension
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( SurfaceIntegralForm ), allocatable :: &
    SI

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SurfaceIntegral_Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( G )
  call G % Initialize ( A )

  allocate ( SI )
  call SI % Initialize &
         ( G, nIntegrals = 1, IgnorabilityOption = CONSOLE % INFO_1 )

  call A % Show ( )
  call G % Show ( )

  associate ( C  =>  A % Chart_GS )

  do iD  =  1,  C % nDimensions
    associate &
      ( I_I => SI % Integrand ( 1, C % Connectivity % iaInner ( iD ) ), &
        I_O => SI % Integrand ( 1, C % Connectivity % iaOuter ( iD ) ) )
      select case ( iD )
      case ( 1, 2 ) !-- outward normal
        I_I % Value  =  -1.0_KDR
        I_O % Value  =  +1.0_KDR
      case ( 3 )  !-- must be continuous at periodic boundary
        I_I % Value  =  +1.0_KDR
        I_O % Value  =  +1.0_KDR
      end select !-- iD
    end associate !-- I_I, I_O
  end do !-- iD

  call SI % Compute ( )

  call Show ( 4.0_KDR  *  CONSTANT % PI  *  C % MaxCoordinate ( 1 ) ** 2, &
              'Expected' )   

  end associate !-- C

  deallocate ( SI )
  deallocate ( G )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program SurfaceIntegral_Form_Test
