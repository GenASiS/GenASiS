program VolumeIntegral_Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use Integrals

  implicit none

  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( FieldSet_A_Form ), allocatable :: &
    IA  !-- Integrand
  type ( VolumeIntegralForm ) :: &
    VI

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'VolumeIntegral__Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( GA )
  call GA % Initialize ( A )

  allocate ( IA )
  call IA % Initialize ( A, NameOption = 'Integrand' )

  call  A % Show ( )
  call GA % Show ( )
  call IA % Show ( )

  associate ( IC  =>  IA % FieldSet_C ( 1 ) % Element )
  associate ( IV  =>  IC % Storage_FSC % Storage % Value ( :, 1 ) )
  IV  =  1.0_KDR
  end associate !-- IV
  end associate !-- IC

  call VI % Compute ( IA, GA, IgnorabilityOption = CONSOLE % INFO_1 )

  deallocate ( IA )
  deallocate ( GA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program VolumeIntegral_Form_Test
