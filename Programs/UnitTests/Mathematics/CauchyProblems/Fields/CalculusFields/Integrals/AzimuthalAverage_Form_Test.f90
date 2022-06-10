program AzimuthalAverage_Form_Test

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
    A_AA
  type ( Geometry_F_Form ), allocatable :: &
    G, &
    G_AA
  type ( FieldSetForm ), allocatable :: &
    FS
  type ( StreamForm ), allocatable :: &
    S, &
    S_AA
  type ( AzimuthalAverageForm ), allocatable :: &
    AA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize ( 'AzimuthalAverage_Form_Test' )

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

  if ( A % Chart_GS_CC % nDimensions  >  2 ) then

    allocate ( A_AA )
    call A_AA % Initialize ( A, nDimensions = 2 )

    allocate ( S_AA )
    call S_AA % Initialize &
           ( A_AA, GIS, NameOption = trim ( S % Name ) // '_AA' )

    allocate ( G_AA )
    call G_AA % Initialize ( A_AA, NameOption = trim ( G % Name ) // '_AA' )
    call G_AA % SetStream ( S_AA )

    allocate ( AA )
    call AA % Initialize ( G, FS, A_AA )

    associate ( FS_AA  =>  AA % FieldSet_AA )

    call S_AA % AddFieldSet ( FS_AA )

    call  A_AA % Show ( )
    call  G_AA % Show ( )
    call FS_AA % Show ( )
    call  S_AA % Show ( )

    end associate !-- FS_AA

  end if !-- nDimensions > 1

  call SetFields ( )

  if ( allocated ( A_AA ) ) &
    call AA % Compute ( IgnorabilityOption = CONSOLE % INFO_1 )

  call GIS % Open ( GIS % ACCESS_CREATE )
  call S    % Write ( )
  if ( allocated ( A_AA ) ) &
    call S_AA % Write ( )
  call GIS % Close ( )

  if ( allocated ( A_AA ) ) then
    deallocate ( AA )
    deallocate ( G_AA )
    deallocate ( S_AA )
    deallocate ( A_AA )
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


end program AzimuthalAverage_Form_Test
