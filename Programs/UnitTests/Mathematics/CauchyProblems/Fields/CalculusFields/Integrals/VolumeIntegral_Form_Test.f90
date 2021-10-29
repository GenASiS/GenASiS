program VolumeIntegral_Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use Integrals

  implicit none

  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( FieldSetForm ), allocatable :: &
    I  !-- Integrand
  type ( VolumeIntegralForm ) :: &
    VI

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'VolumeIntegral_Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( G )
  call G % Initialize ( A )

  allocate ( I )
  call I % Initialize ( A, NameOption = 'Integrand' )

  call A % Show ( )
  call G % Show ( )
  call I % Show ( )

  associate ( IV  =>  I % Storage ( 1 ) % Value ( :, 1 ) )
  IV  =  1.0_KDR
  end associate !-- IV

  call VI % Compute ( I, G, IgnorabilityOption = CONSOLE % INFO_1 )

  deallocate ( I )
  deallocate ( G )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program VolumeIntegral_Form_Test
