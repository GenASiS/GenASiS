program FieldSet_CB__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  implicit none

  integer ( KDI ) :: &
    nFields = 5
  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_BH_Form ), allocatable :: &
    C
  type ( FieldSet_CB_Form ), allocatable :: &
    FSC
  type ( Stream_CH_Form ), allocatable :: &
    SC
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FSM
  type ( Stream_MH_Form ), allocatable :: &
    SM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_CB__Form_Test', DimensionalityOption = '2D' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize_BH &
         ( M, Periodic, iChart = 1 )
  call M % Show ( )
  call C % Show ( )

  call CONSOLE % SetVerbosity ( 'INFO_2' )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize ( M, 'Fields' ) 
  call FSC % Initialize ( FSM, C, 'Fields', nFields ) 

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  allocate ( SM )
  allocate ( SC )
  call SM % Initialize ( M, GIS, 'Stream' )
  call SC % Initialize ( SM, C, GIS, 'Stream' )

  call SM % AddFieldSet ( FSM )
  call SC % AddFieldSet ( FSC )

  call FSC % AddStream ( SC )
  call FSC % AddStream ( SC )  !-- Test the prevention of duplication

  deallocate ( SC )
  deallocate ( SM )
  deallocate ( GIS )
  deallocate ( FSC )
  deallocate ( FSM )

  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program FieldSet_CB__Form_Test

