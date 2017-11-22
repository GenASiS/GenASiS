program FluidCentralCore_Form_Test

  use Basics
  use FluidCentralCore_Form

  implicit none

  type ( FluidCentralCoreForm ), allocatable :: &
    FB

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'FluidCentralCore_Form_Test' )

  allocate ( FB )
  call FB % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'DUST', &
           GeometryType = 'NEWTONIAN', &
           GravitySolverTypeOption = 'MULTIPOLE' )
  deallocate ( FB )

  deallocate ( PROGRAM_HEADER )

end program FluidCentralCore_Form_Test
