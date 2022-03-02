program VolumeIntegral_Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use Integrals

  implicit none

  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( VolumeIntegralForm ), allocatable :: &
    VI

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'VolumeIntegral_Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( G )
  call G % Initialize ( A )

  allocate ( VI )
  call VI % Initialize &
         ( G, nIntegrals = 1, IgnorabilityOption = CONSOLE % INFO_1 )

  call A % Show ( )
  call G % Show ( )

  associate ( I   =>  VI % Integrand )
  associate ( IV  =>  I % Storage ( 1 ) % Value ( :, 1 ) )
  IV  =  1.0_KDR
  end associate !-- IV
  end associate !-- I

  call VI % Compute ( )

  associate ( C  =>  A % Chart_GS )
  call Show ( 4.0_KDR / 3.0_KDR  *  CONSTANT % PI  &
                *  C % MaxCoordinate ( 1 ) ** 3, &
              'Expected' )   
  end associate !-- C

  deallocate ( VI )
  deallocate ( G )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program VolumeIntegral_Form_Test
