program Stream_CH__Form_Test

  !-- Stream_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( FieldSet_CH_Form ), allocatable :: &
    FSC
  type ( Stream_CH_Form ), allocatable :: &
    SC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize_H ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( FSC )
  call FSC % Initialize_H ( C )

  allocate ( SC )
  call SC % Initialize_H ( C, GIS )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SC % AddFieldSet ( FSC )
  call SC % AddFieldSet ( FSC )  !-- Test the prevention of duplication
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call  C % Show ( )
  call SC % Show ( )

  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Stream_CH__Form_Test
