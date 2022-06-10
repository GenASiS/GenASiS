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
  type ( StreamForm ), allocatable :: &
    S, &
    S_SA
  type ( SphericalAverageForm ), allocatable :: &
    SA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'SphericalAverage_Form_Test' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator, &
           NameOption = 'PositionSpace' )

  allocate ( S )
  call S % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( S )

  allocate ( FS )
  call FS % Initialize &
         ( A, &
           FieldOption = [ 'Sphere   ', 'Spheroid ', 'Ellipsoid' ], &
           nFieldsOption = 3 )
  call S % AddFieldSet ( FS )

  call  A    % Show ( )
  call  G    % Show ( )
  call FS    % Show ( )
  call  S    % Show ( )

  if ( A % Chart_GS_CC % nDimensions  >  1 ) then

    allocate ( A_SA )
    call A_SA % Initialize ( A, nDimensions = 1 )

    allocate ( S_SA )
    call S_SA % Initialize &
           ( A_SA, GIS, NameOption = trim ( S % Name ) // '_SA' )

    allocate ( G_SA )
    call G_SA % Initialize ( A_SA, NameOption = trim ( G % Name ) // '_SA' )
    call G_SA % SetStream ( S_SA )

    allocate ( SA )
    call SA % Initialize ( G, FS, A_SA )

    associate ( FS_SA  =>  SA % FieldSet_SA )

    call S_SA % AddFieldSet ( FS_SA )

    call  A_SA % Show ( )
    call  G_SA % Show ( )
    call FS_SA % Show ( )
    call  S_SA % Show ( )

    end associate !-- FS_SA

  end if !-- nDimensions > 1

  call SetFields ( )

  if ( allocated ( A_SA ) ) &
    call SA % Compute ( IgnorabilityOption = CONSOLE % INFO_1 )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S    % Write ( )
  if ( allocated ( A_SA ) ) &
    call S_SA % Write ( )
  call GIS % Close ( )

  if ( allocated ( A_SA ) ) then
    deallocate ( SA )
    deallocate ( G_SA )
    deallocate ( S_SA )
    deallocate ( A_SA )
  end if
  deallocate ( FS )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetFields ( )

    real ( KDR ) :: &
      A1, A2, A3
    real ( KDR ), dimension ( : ), allocatable :: &
      X, Y, Z

    associate &
      ( C      =>  A % Chart_GS_CC, &
        nV     =>  G % Storage_GS % nValues, &
        R      =>  G % Storage_GS % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta  =>  G % Storage_GS % Value ( :, G % CENTER_U ( 2 ) ), &
        Phi    =>  G % Storage_GS % Value ( :, G % CENTER_U ( 3 ) ), &
        Sphere     =>  FS % Storage_GS % Value ( :, 1 ), &
        Spheroid   =>  FS % Storage_GS % Value ( :, 2 ), &
        Ellipsoid  =>  FS % Storage_GS % Value ( :, 3 ) )

    A1  =  0.6_KDR  *  C % MaxCoordinate ( 1 )
    A2  =  0.4_KDR  *  C % MaxCoordinate ( 1 )
    A3  =  0.2_KDR  *  C % MaxCoordinate ( 1 )

    allocate ( X ( nV ), Y ( nV ), Z ( nV ) )

    select case ( C % nDimensions )
    case ( 1 )
      X  =  R
      Y  =  0.0_KDR
      Z  =  0.0_KDR
    case ( 2 )
      X  =  R * Sin ( Theta )
      Y  =  0.0_KDR
      Z  =  R * Cos ( Theta )
    case ( 3 )
      X  =  R * Sin ( Theta ) * Cos ( Phi )
      Y  =  R * Sin ( Theta ) * Sin ( Phi )
      Z  =  R * Cos ( Theta )
    end select !-- nDimensions

    where ( ( X ** 2  +  Y ** 2  +  Z ** 2 ) / A1 ** 2  <  1.0_KDR )
      Sphere  =  1.0_KDR
    elsewhere
      Sphere  =  0.0_KDR
    end where

    where ( ( X ** 2  +  Y ** 2 ) / A1 ** 2  &
            +  Z ** 2 / A2 ** 2  <  1.0_KDR )
      Spheroid  =  1.0_KDR
    elsewhere
      Spheroid  =  0.0_KDR
    end where

    where ( X ** 2 / A1 ** 2  +  Y ** 2 / A2 ** 2 &
            +  Z ** 2 / A3 ** 2  <  1.0_KDR )
      Ellipsoid  =  1.0_KDR
    elsewhere
      Ellipsoid  =  0.0_KDR
    end where

    end associate !-- C, etc.

  end subroutine SetFields


end program SphericalAverage_Form_Test
