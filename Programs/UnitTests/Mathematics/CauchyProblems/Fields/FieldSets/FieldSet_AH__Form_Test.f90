program FieldSet_AH__Form_Test

  !-- FieldSet_AtlasHeader__Form_Test

  use Basics
  use Manifolds
  use FieldSets

  implicit none

  type ( Atlas_SG_Form ), allocatable :: &
    A
  type ( FieldSet_AH_Form ), allocatable :: &
    FSA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_AH__Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSA )
  call FSA % Initialize_H ( A )
  select type ( G  =>  A % Chart ( 1 ) % Element )
  class is ( Grid_S_Form )

  allocate ( FieldSet_GS_Form :: FSA % FieldSet_C ( 1 ) % Element )
  select type ( FSG  =>  FSA % FieldSet_C ( 1 ) % Element )
  class is ( FieldSet_GS_Form )
  call FSG % Initialize ( G )

  end select !-- FSG
  end select !-- G

  call   A % Show ( )
  call FSA % Show ( )

  deallocate ( FSA )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_AH__Form_Test
