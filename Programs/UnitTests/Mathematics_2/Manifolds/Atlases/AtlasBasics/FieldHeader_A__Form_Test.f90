program FieldHeader_A__Form_Test

  use Basics
  use AtlasBasics

  type ( AtlasHeaderForm ), allocatable :: &
    A
  type ( FieldHeader_A_Form ), allocatable :: &
    F

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'AtlasHeader_Form_Test', DimensionalityOption = '2D' )

  allocate ( A )
  call A % Initialize &
         ( 'Atlas', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call A % Show ( )

  allocate ( F )
  call F % Initialize ( A, 'Field' ) 

  deallocate ( F )
  deallocate ( A )
  deallocate ( PROGRAM_HEADER )

end program FieldHeader_A__Form_Test
