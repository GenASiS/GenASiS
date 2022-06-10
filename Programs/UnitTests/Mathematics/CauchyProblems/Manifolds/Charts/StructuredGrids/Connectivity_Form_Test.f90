program Connectivity_Form_Test

  use Basics
  use Connectivity_Form

  implicit none

  character ( LDF ), parameter :: &
    Name = 'Connectivity_Form_Test'

  integer ( KDI ) :: &
    iD  !-- iDimension
  type ( ConnectivityForm ), allocatable :: &
    Connectivity

  allocate ( PROGRAM_HEADER )
  
  call PROGRAM_HEADER % Initialize ( Name )

  do iD = 1, 3

    allocate ( Connectivity )
    call Connectivity % Initialize ( nDimensions = iD )
    call Connectivity % Show ( nLeadingLinesOption = 2 )
    deallocate ( Connectivity )

    allocate ( Connectivity )
    call Connectivity % Initialize ( nDimensions = iD, &
                                     IncludeEdgesOption = .true. )
    call Connectivity % Show ( nLeadingLinesOption = 2 )
    deallocate ( Connectivity )

  end do !-- iD

  deallocate ( PROGRAM_HEADER )

end program Connectivity_Form_Test
