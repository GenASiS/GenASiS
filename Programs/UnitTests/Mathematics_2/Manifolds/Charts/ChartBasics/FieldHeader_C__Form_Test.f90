program FieldHeader_C__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
  type ( ChartHeaderForm ), allocatable :: &
    C
  type ( FieldHeader_C_Form ), allocatable :: &
    FC
  type ( ManifoldHeaderForm ), allocatable :: &
    M
  type ( FieldHeader_M_Form ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldHeader_C__Form_Test', DimensionalityOption = '2D' )

  IsPeriodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize &
         ( M, IsPeriodic, iChart = 1 )
  call M % Show ( )
  call C % Show ( )

  allocate ( FM )
  allocate ( FC )
  call FM % Initialize ( M, 'Field' ) 
  call FC % Initialize ( FM, C, 'Field' ) 

  deallocate ( FM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldHeader_C__Form_Test
