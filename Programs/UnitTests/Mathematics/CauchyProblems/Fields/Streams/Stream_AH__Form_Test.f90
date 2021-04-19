program Stream_AH__Form_Test

  !-- Stream_AtlasHeader__Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Streams

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_H_Form ), allocatable :: &
    A
  type ( FieldSet_AH_Form ), allocatable :: &
    FSA
  type ( Stream_AH_Form ), allocatable :: &
    SA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_AH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize_H ( )
  allocate ( A % Chart ( 1 ) % Element )
  associate ( C  =>  A % Chart ( 1 ) % Element )
  call C % Initialize_H &
         ( PeriodicOption = [ .true., .true., .true. ], &
           iDimensionalityOption = 1 )

  allocate ( FSA )
  call FSA % Initialize_H ( A )
  allocate ( FSA % FieldSet_C ( 1 ) % Element )
  associate ( FSC  =>  FSA % FieldSet_C ( 1 ) % Element )
  call FSC % Initialize_H ( C )

  allocate ( SA )
  call SA % Initialize_H ( A )
  allocate ( SA % Stream_C ( 1 ) % Element )
  associate ( SC  =>  SA % Stream_C ( 1 ) % Element )
  call SC % Initialize_H ( C, GIS )

  end associate !-- SC
  end associate !-- FSC
  end associate !-- C

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SA % AddFieldSet ( FSA )
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   A % Show ( )
  call FSA % Show ( )
  call  SA % Show ( )

  deallocate ( SA )
  deallocate ( FSA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Stream_AH__Form_Test
