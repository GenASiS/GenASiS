program Geometry_F_AH__Form_Test

  !-- Geometry_F_AtlasHeader__Form_Test

  use Basics
  use Manifolds
  use Streams
  use Geometries

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_H_Form ), allocatable :: &
    A
  type ( Stream_AH_Form ), allocatable :: &
    SA
  type ( Geometry_F_AH_Form ), allocatable :: &
    GA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Geometry_AH__Form_Test', DimensionalityOption = '2D' )

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

  allocate ( SA )
  call SA % Initialize_H ( A )
  allocate ( SA % Stream_C ( 1 ) % Element )
  associate ( SC  =>  SA % Stream_C ( 1 ) % Element )
  call SC % Initialize_H ( C, GIS )

  allocate ( GA )
  call GA % Initialize_H ( A )
  allocate ( GA % Geometry_C ( 1 ) % Element )
  associate ( GC  =>  GA % Geometry_C ( 1 ) % Element )
  call GC % Initialize_H ( C )

  end associate !-- GC
  end associate !-- SC
  end associate !--  C

  call GA % SetStream ( SA )

  call  A % Show ( )
  call GA % Show ( )
  call SA % Show ( )

  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program Geometry_F_AH__Form_Test
