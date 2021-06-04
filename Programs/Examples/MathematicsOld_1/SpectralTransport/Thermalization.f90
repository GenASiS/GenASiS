program Thermalization

  use Basics
  use Thermalization_Form

  implicit none

  type ( ThermalizationForm ), allocatable :: &
    T

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Thermalization', DimensionalityOption = '2D_1D' )

  allocate ( T )
  call T % Initialize ( PROGRAM_HEADER % Name )
  call T % Evolve ( )
  deallocate ( T )

  deallocate ( PROGRAM_HEADER )

end program Thermalization
