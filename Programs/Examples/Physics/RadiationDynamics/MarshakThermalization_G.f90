program MarshakThermalization_G

  use Basics
  use MarshakThermalization_Form

  implicit none

  type ( MarshakThermalizationForm ), allocatable :: &
    MT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'MarshakThermalization_G', DimensionalityOption = '1D' )

  allocate ( MT )
  call MT % Initialize_MT ( 'GREY', PROGRAM_HEADER % Name )
  call MT % Evolve ( )
  deallocate ( MT )

  deallocate ( PROGRAM_HEADER )

end program MarshakThermalization_G
