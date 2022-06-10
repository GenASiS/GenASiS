program SineWaveDiffusion

  use Basics
  use SineWaveDiffusion_Form

  implicit none

  type ( SineWaveDiffusionForm ), allocatable :: &
    SWD

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'SineWaveDiffusion', DimensionalityOption = '1D' )

  allocate ( SWD )
  call SWD % Initialize ( PROGRAM_HEADER % Name )
  call SWD % Evolve ( )
  call SWD % ComputeError ( )
  deallocate ( SWD )

  deallocate ( PROGRAM_HEADER )

end program SineWaveDiffusion
