program VolumeIntegral_Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Integrals

  implicit none

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'VolumeIntegral__Form_Test', DimensionalityOption = '2D' )

  deallocate ( PROGRAM_HEADER )

end program VolumeIntegral_Form_Test
