program FieldSet_CB__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics
  use BaseCharts

  implicit none

  integer ( KDI ) :: &
    nFields = 5
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
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
  call C % Initialize &
         ( M, 'Global', Periodic )

  call CONSOLE % SetVerbosity ( 'INFO_3' )

  DeviceMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( PinnedMemory, 'PinnedMemory' )

  DevicesCommunicate  =  OffloadEnabled ( )  .and.  GetNumberOfDevices ( ) >= 1
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize &
         ( M, 'Fields', &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate )
  call FSC % Initialize &
         ( C, FSM, nFields ) 

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  allocate ( SM )
  allocate ( SC )
  call SM % Initialize ( M, GIS, 'Stream' )
  call SC % Initialize ( C, SM )

  call SM % AddFieldSet ( FSM )

  call FSC % AddStream ( SC )
  call FSC % AddStream ( SC )  !-- Test the prevention of duplication

  call M % Show ( )
  call Show ( M % nCharts,    'nCharts',    M % IGNORABILITY )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call FSM % Show ( )
  call Show ( FSM % nStreams,  'nStreams',  FSM % IGNORABILITY )

  call C % Show ( )
  call FSC % Show ( )

  call InitializeField ( FSC )

  deallocate ( SC )
  deallocate ( SM )
  deallocate ( GIS )
  deallocate ( FSC )
  deallocate ( FSM )

  call CONSOLE % SetVerbosity ( 'INFO_1' )

  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

contains


  subroutine InitializeField ( FSC )

    class ( FieldSet_CB_Form ), intent ( inout ) :: &
      FSC

    call FSC % ExchangeGhostData ( )

  end subroutine InitializeField


end program FieldSet_CB__Form_Test

