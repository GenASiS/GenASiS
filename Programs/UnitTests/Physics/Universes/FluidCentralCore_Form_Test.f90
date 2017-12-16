program FluidCentralCore_Form_Test

  use Basics
  use FluidCentralCore_Form

  implicit none

  type ( FluidCentralCoreForm ), allocatable :: &
    FCC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'FluidCentralCore_Form_Test' )

  allocate ( FCC )
  call FCC % Initialize &
         ( PROGRAM_HEADER % Name, &
           FluidType = 'DUST', &
           GeometryType = 'NEWTONIAN' )
  deallocate ( FCC )

  deallocate ( PROGRAM_HEADER )

end program FluidCentralCore_Form_Test
