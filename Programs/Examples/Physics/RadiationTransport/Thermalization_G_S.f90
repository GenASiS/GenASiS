program Thermalization_G_S

  !-- Thermalization_Generic_Spectral

  use GenASiS
  use Thermalization_G__Form

  implicit none

  type ( Thermalization_G_Form ), allocatable :: &
    T

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Thermalization_G_S', DimensionalityOption = '2D_1D' )

  allocate ( T )
  call T % Initialize ( 'SPECTRAL', PROGRAM_HEADER % Name )
  call T % Evolve ( )
  deallocate ( T )

  deallocate ( PROGRAM_HEADER )

end program Thermalization_G_S
