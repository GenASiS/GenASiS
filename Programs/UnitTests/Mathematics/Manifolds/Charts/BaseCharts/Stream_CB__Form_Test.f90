program Stream_CB__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  logical ( KDL ), dimension ( 3 ) :: &
    IsPeriodic
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    GM
  type ( FieldSet_CH_Form ), allocatable :: &
    GC
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( Chart_BH_Form ), allocatable :: &
    C
  type ( Stream_CB_Form ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_CB__Form_Test', DimensionalityOption = '2D' )

  IsPeriodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize_BH &
         ( M, IsPeriodic, iChart = 1 )
  call M % Show ( )
  call C % Show ( )

  allocate ( GM )
  allocate ( GC )
  call GM % Initialize ( M, 'Geometry' ) 
  call GC % Initialize ( GM, C, 'Geometry' ) 

  call CONSOLE % SetVerbosity ( 'INFO_4' )

  allocate ( G )
  call G % Initialize &
         ( GC, nValues = C % nValues, NameOption = 'Geometry_F' )

  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = C % Communicator )

  allocate ( SC )
  call SC % Initialize ( C, GIS, 'Stream' )

  deallocate ( SC )
  deallocate ( GIS )

  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( G )
  deallocate ( GC )
  deallocate ( GM )
  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Stream_CB__Form_Test
