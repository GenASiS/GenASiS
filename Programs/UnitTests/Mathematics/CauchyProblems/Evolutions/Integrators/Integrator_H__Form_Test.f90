program Integrator_H__Form_Test

  !-- Integrator_Header__Form_Test

  use Basics
  use Integrators

  implicit none

  type ( Integrator_H_Form ), allocatable :: &
    I

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Integrator_H__Form_Test', DimensionalityOption = '2D' )

  allocate ( I )
  call I % Initialize ( )
  call I % Evolve ( )

  deallocate ( I )
  deallocate ( PROGRAM_HEADER )

end program Integrator_H__Form_Test
