program DynamicDiffusion

  use Basics
  use DynamicDiffusion_Form

  implicit none

  type ( DynamicDiffusionForm ), allocatable :: &
    DD

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'DynamicDiffusion' )!, DimensionalityOption = '2D' )

  allocate ( DD )
  call DD % Initialize ( PROGRAM_HEADER % Name )
  call DD % Evolve ( )
  call DD % ComputeError ( )
  deallocate ( DD )

  deallocate ( PROGRAM_HEADER )

end program DynamicDiffusion
