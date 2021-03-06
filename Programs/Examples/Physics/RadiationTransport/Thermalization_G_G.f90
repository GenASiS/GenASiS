program Thermalization_G_G

  !-- Thermalization_Generic_Grey

  use GenASiS
  use Thermalization_G__Form

  implicit none

  type ( Thermalization_G_Form ), allocatable :: &
    T

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Thermalization_G_G', DimensionalityOption = '2D' )

  allocate ( T )
  call T % Initialize ( 'GREY', PROGRAM_HEADER % Name )
  call T % Evolve ( )
  deallocate ( T )

  deallocate ( PROGRAM_HEADER )

end program Thermalization_G_G
