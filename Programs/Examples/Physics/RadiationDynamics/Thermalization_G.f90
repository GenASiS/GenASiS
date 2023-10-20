program Thermalization_G

  !-- Thermalization_Grey

  use GenASiS
  use Thermalization_Form

  implicit none

  type ( ThermalizationForm ), allocatable :: &
    T

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Thermalization_G', DimensionalityOption = '2D' )

  allocate ( T )
  call T % Initialize ( 'GREY', PROGRAM_HEADER % Name )
!call T % Show ( )
!  call T % Evolve ( )
  deallocate ( T )

  deallocate ( PROGRAM_HEADER )

end program Thermalization_G
