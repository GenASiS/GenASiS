program ThermalizationNeutrinos

  use Basics
  use ThermalizationNeutrinos_Form

  implicit none

  type ( ThermalizationNeutrinosForm ), allocatable :: &
    T

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'ThermalizationNeutrinos', DimensionalityOption = '2D_1D' )

  allocate ( T )
  call T % Initialize ( PROGRAM_HEADER % Name )
!  call T % Evolve ( )
  deallocate ( T )

  deallocate ( PROGRAM_HEADER )

end program ThermalizationNeutrinos
