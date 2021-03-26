program FieldSet_CB__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  implicit none

  integer ( KDI ) :: &
    nFields = 5
  logical ( KDL ) :: &
    Pinned, &
    UseDeviceGhost
  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_BH_Form ), allocatable :: &
    C
  type ( FieldSet_CB_Form ), allocatable :: &
    GC, &
    FSC
  type ( Stream_CH_Form ), allocatable :: &
    SC
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    GM, &
    FSM
  type ( Stream_MH_Form ), allocatable :: &
    SM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'FieldSet_CB__Form_Test', DimensionalityOption = '2D' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( GM )
  allocate ( C )
  allocate ( GC )
  allocate ( G )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call GM % Initialize &
         ( M, 'Geometry' ) 
  call C % Initialize &
         ( M, 'Global', Periodic )
  call GC % Initialize &
         ( C, GM )
  call G % Initialize &
         ( GC, nValues = C % nValues )
  call C % ComputeGeometry ( G )
  call M % Show ( )
  call C % Show ( )

  call CONSOLE % SetVerbosity ( 'INFO_2' )

  Pinned  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( Pinned, 'Pinned' )

  UseDeviceGhost  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter ( UseDeviceGhost, 'UseDeviceGhost' )

call Show ( Pinned, '>>> Pinned' )
call Show ( UseDeviceGhost, '>>> UseDeviceGhost' )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize ( M, 'Fields' ) 
  call FSC % Initialize &
         ( C, FSM, 'Fields', nFields, PinnedOption = Pinned, &
           UseDeviceGhostOption = UseDeviceGhost ) 
  call FSC % ExchangeGhostData ( )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  allocate ( SM )
  allocate ( SC )
  call SM % Initialize ( M, GIS, 'Stream' )
  call SC % Initialize ( C, SM )

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

