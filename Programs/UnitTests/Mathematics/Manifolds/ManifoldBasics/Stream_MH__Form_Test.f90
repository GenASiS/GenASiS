program Stream_MH__Form_Test

  use Basics
  use ManifoldBasics

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Manifold_H_Form ), allocatable :: &
    M
  type ( FieldSet_MH_Form ), allocatable :: &
    FSM
  type ( Stream_MH_Form ), allocatable :: &
    SM

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Stream_MH__Form_Test', DimensionalityOption = '2D' )

  allocate ( M )
  call M % Initialize &
         ( 'Manifold', CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( FSM )
  call FSM % Initialize ( M, 'Fields' ) 

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, CommunicatorOption = M % Communicator )

  allocate ( SM )
  call SM % Initialize ( M, GIS, 'Stream' )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call SM % AddFieldSet ( FSM )
  call SM % AddFieldSet ( FSM )  !-- Test the prevention of duplication
  call CONSOLE % SetVerbosity ( 'INFO_1' )

  call   M % Show ( )
  call FSM % Show ( )
  call  SM % Show ( )

  deallocate ( SM )
  deallocate ( GIS )

  deallocate ( FSM )
  deallocate ( M )
  deallocate ( PROGRAM_HEADER )

end program Stream_MH__Form_Test
