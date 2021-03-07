program FieldsHeader_C__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
  type ( ChartHeaderForm ), allocatable :: &
    C
  type ( FieldsHeader_C_Form ), allocatable :: &
    FC
  type ( ManifoldHeaderForm ), allocatable :: &
    M
  type ( FieldsHeader_M_Form ), allocatable :: &
    FM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldsHeader_C__Form_Test', DimensionalityOption = '2D' )

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
  call FM % Initialize ( M, 'Fields' ) 
  call FC % Initialize ( FM, C, 'Fields' ) 

  deallocate ( FM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldsHeader_C__Form_Test
