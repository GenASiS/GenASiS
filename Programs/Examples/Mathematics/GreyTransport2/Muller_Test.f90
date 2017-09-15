program Muller_Test

  use Basics
  use Muller_Test_Form

  implicit none

  type ( MullerTestForm ), allocatable :: &
    MT

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'MullerTest', DimensionalityOption = '1D' )

  allocate ( MT )
  call MT % Initialize ( PROGRAM_HEADER % Name )
  call MT % Evolve ( )
  deallocate ( MT )

  deallocate ( PROGRAM_HEADER )

end program Muller_Test
