program FluidSymmetricCurvilinear_Form_Test

  use Basics
  use FluidSymmetricCurvilinear_Form

  implicit none

  type ( FluidSymmetricCurvilinearForm ), allocatable :: &
    FSC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'FluidSymmetricCurvilinear_Form_Test' )

  allocate ( FSC )
  call FSC % Initialize ( PROGRAM_HEADER % Name, FluidType = 'IDEAL' )
  deallocate ( FSC )

  deallocate ( PROGRAM_HEADER )

end program FluidSymmetricCurvilinear_Form_Test
