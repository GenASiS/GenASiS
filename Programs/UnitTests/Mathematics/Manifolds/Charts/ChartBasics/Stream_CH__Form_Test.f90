program Stream_CH__Form_Test

  use Basics
  use ManifoldBasics
  use ChartBasics

  implicit none

  integer ( KDI ) :: &
    iFS, &  !-- iFieldSet
    nFields = 5
  logical ( KDL ), dimension ( 3 ) :: &
    Periodic
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
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
         ( 'FieldSet_CH__Form_Test', DimensionalityOption = '2D' )

  Periodic  =  .true.

  allocate ( M )
  allocate ( C )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )
  call C % Initialize_H &
         ( M, 'Global', Periodic )

  allocate ( FSM )
  allocate ( FSC )
  call FSM % Initialize ( M, 'Fields', nFields ) 
  call FSC % Initialize ( C, FSM ) 

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  allocate ( SM )
  allocate ( SC )
  call SM % Initialize ( M, GIS, 'Stream' )
  call SC % Initialize ( C, SM )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SM % AddFieldSet ( FSM )
  call SC % AddFieldSet ( FSC )
  call SC % AddFieldSet ( FSC )  !-- Test the prevention of duplication
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call M % Show ( )
  call Show ( M % nCharts,    'nCharts',    M % IGNORABILITY )
  call Show ( M % nFieldSets, 'nFieldSets', M % IGNORABILITY )
  call Show ( M % nStreams,   'nStreams',   M % IGNORABILITY )

  call FSM % Show ( )
  call Show ( FSM % nStreams,  'nStreams',  FSM % IGNORABILITY )

  call SM % Show ( )
  call Show ( SM % nFieldSets, 'nFieldSets', SM % IGNORABILITY )
  do iFS  =  1, SM % nFieldSets
    associate ( FS  =>  SM % FieldSet ( iFS ) % Pointer )
    call Show ( FS % Name, 'FieldSet', SM % IGNORABILITY )
    end associate !-- FS
  end do !-- iFS

  call C % Show ( )
  call Show ( C % nFieldSets, 'nFieldSets', C % IGNORABILITY )
  call Show ( C % nStreams,   'nStreams',   C % IGNORABILITY )

  call FSC % Show ( )
  call Show ( FSC % nStreams,  'nStreams',  FSC % IGNORABILITY )

  call SC % Show ( )
  call Show ( SC % nFieldSets, 'nFieldSets', SC % IGNORABILITY )
  do iFS  =  1, SC % nFieldSets
    associate ( FS  =>  SC % FieldSet ( iFS ) % Pointer )
    call Show ( FS % Name, 'FieldSet', SC % IGNORABILITY )
    end associate !-- FS
  end do !-- iFS

  deallocate ( SC )
  deallocate ( SM )
  deallocate ( GIS )

  deallocate ( FSC )
  deallocate ( FSM )
  deallocate ( C )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Stream_CH__Form_Test
