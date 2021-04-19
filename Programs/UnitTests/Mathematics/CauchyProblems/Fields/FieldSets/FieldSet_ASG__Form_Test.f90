program FieldSet_ASG__Form_Test

  !-- FieldSet_AtlasSingleGrid__Form_Test

  use Basics
  use Manifolds
  use FieldSets

  implicit none

  type ( Atlas_SG_Form ), allocatable :: &
    A
  type ( FieldSet_ASG_Form ), allocatable :: &
    FSA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_ASG__Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSA )
  call FSA % Initialize ( A )

  call   A % Show ( )
  call FSA % Show ( )

  deallocate ( FSA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_ASG__Form_Test
