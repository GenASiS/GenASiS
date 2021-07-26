program SphericalAverage_Form_Test

  use Basics
  use Manifolds
  use FieldSets
  use Geometries
  use Integrals

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A, &
    A_SA
  type ( Geometry_F_Form ), allocatable :: &
    G, &
    G_SA
  type ( FieldSetForm ), allocatable :: &
    FS
  type ( SphericalAverageForm ), allocatable :: &
    FS_SA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SphericalAverage_Form_Test' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  allocate ( A_SA )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )
  call A_SA % Initialize ( A )

  allocate ( G )
  allocate ( G_SA )
  call G % Initialize ( A )
  call G_SA % Initialize ( A_SA )

  allocate ( FS )
  allocate ( FS_SA )
  call FS % Initialize &
         ( A, &
           FieldOption = [ 'Sphere   ', 'Spheroid ', 'Ellipsoid' ], &
           NameOption = 'Integrand', &
           nFieldsOption = 3 )
  call FS_SA % Initialize ( G, FS, A_SA )

  call A     % Show ( )
  call A_SA  % Show ( )
  call G     % Show ( )
  call G_SA  % Show ( )
  call FS    % Show ( )
  call FS_SA % Show ( )

  deallocate ( FS_SA )
  deallocate ( FS )
  deallocate ( G_SA )
  deallocate ( G )
  deallocate ( A_SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )

end program SphericalAverage_Form_Test
