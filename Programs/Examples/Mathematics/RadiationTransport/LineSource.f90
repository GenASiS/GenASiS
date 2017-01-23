program LineSource

  use Basics
  use LineSource_Form

  implicit none

  type ( LineSourceForm ), allocatable :: &
    LS

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'LineSource', DimensionalityOption = '1D' )

  allocate ( LS )
  call LS % Initialize ( PROGRAM_HEADER % Name )
  call LS % Evolve ( )
  call LS % ComputeError ( )
  deallocate ( LS )

  deallocate ( PROGRAM_HEADER )

end program LineSource
