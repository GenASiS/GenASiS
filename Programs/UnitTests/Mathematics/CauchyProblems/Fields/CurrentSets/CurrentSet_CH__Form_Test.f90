program CurrentSet_CH__Form_Test

  !-- CurrentSet_ChartHeader__Form_Test

  use Basics
  use Manifolds
  use Streams
  use CurrentSets

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Chart_H_Form ), allocatable :: &
    C
  type ( Stream_CH_Form ), allocatable :: &
    SC
  type ( CurrentSet_CH_Form ), allocatable :: &
    CSC

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'CurrentSet_CH__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( C )
  call C % Initialize_H &
         ( PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SC )
  call SC % Initialize_H ( C, GIS )

  allocate ( CSC )
  call CSC % Initialize_H ( C )
!  call GC % SetStream ( SC )

  call  C % Show ( )
  ! call GC % Show ( )
  call SC % Show ( )

  deallocate ( CSC )
  deallocate ( SC )
  deallocate ( C )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program CurrentSet_CH__Form_Test
