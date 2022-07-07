program Laplacian_M_ASCG__Form_Test

  !-- Laplacian_Multipole_AtlasSingleChartGrid__Form_Test

  use Basics
  use Manifolds
  use Fields
  use PoissonEquations

  implicit none

  integer ( KDI ) :: &
    nEquations, &
    MaxDegree
  logical ( KDL ) :: &
    DeviceMemory, &
    PinnedMemory, &
    DevicesCommunicate
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    S
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( Laplacian_M_ASCG_Form ), allocatable :: &
    L

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Laplacian_M_ASCG__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( S )
  call S % Initialize ( A, GIS )
  
  DeviceMemory  =  OffloadEnabled ( )  .and.  NumberOfDevices ( ) >= 1 
  call PROGRAM_HEADER % GetParameter ( DeviceMemory, 'DeviceMemory' )

  PinnedMemory        =  DeviceMemory
  DevicesCommunicate  =  DeviceMemory
  call PROGRAM_HEADER % GetParameter &
         ( PinnedMemory, 'PinnedMemory' )
  call PROGRAM_HEADER % GetParameter &
         ( DevicesCommunicate, 'DevicesCommunicate' )

  allocate ( G )
  call G % Initialize &
         ( A, &
           DeviceMemoryOption = DeviceMemory, &
           PinnedMemoryOption = PinnedMemory, &
           DevicesCommunicateOption = DevicesCommunicate )
  call G % SetStream ( S )

  nEquations = 1

  MaxDegree = 3
  call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

  allocate ( L )
  call L % Initialize ( G, MaxDegree, nEquations )

  call A % Show ( )
  call G % Show ( )
  call L % Show ( )

  call TestAssociatedLegendre ( )
  call TestAngularFunctions ( )
  call TestHomogeneousSphere ( )

  deallocate ( L )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestAssociatedLegendre ( )

    real ( KDR ) :: &
      X_Random, &
      Cos_X, Sin_X, &
      Pi

    call Show ( 'Testing normalized Associated Legendre polynomials' )

    call InitializeRandomSeed ( PROGRAM_HEADER % Communicator )
    call random_number ( X_Random )
    X_Random  =  -1.0_KDR  +  2.0 * X_Random
       Cos_X  =  X_Random
       Sin_X  =  sqrt ( 1.0_KDR  -  X_Random ** 2 )
          Pi  =  CONSTANT % PI

    call Show ( X_Random, 'X_Random' )

    call Show ( L % AssociatedLegendre ( X_Random, 0, 0 ), 'P_0_0 computed' )
    call Show ( sqrt ( 1.0_KDR / ( 4.0_KDR * Pi ) ), 'P_0_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 1, 0 ), 'P_1_0 computed' )
    call Show ( sqrt ( 3.0_KDR / ( 4.0_KDR * Pi ) ) * Cos_X, &
                'P_1_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 1, 1 ), 'P_1_1 computed' )
    call Show ( - sqrt ( 3.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X, &
                'P_1_1 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 0 ), 'P_2_0 computed' )
    call Show ( sqrt ( 5.0_KDR / ( 4.0_KDR * Pi ) ) &
                *  (    ( 3.0_KDR / 2.0_KDR )  *  Cos_X ** 2  &
                     -  ( 1.0_KDR / 2.0_KDR ) ), &
                'P_2_0 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 1 ), 'P_2_1 computed' )
    call Show ( - sqrt ( 15.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X * Cos_X, &
                'P_2_1 expected' )

    call Show ( L % AssociatedLegendre ( X_Random, 2, 2 ), 'P_2_2 computed' )
    call Show ( ( 1.0_KDR / 4.0_KDR ) * sqrt ( 15.0_KDR / ( 2.0_KDR * Pi ) ) &
                *  Sin_X ** 2, &
                'P_2_2 expected' )

  end subroutine TestAssociatedLegendre


  subroutine TestAngularFunctions ( )

    integer ( KDI ) :: &
     iR  !-- iRank
    integer ( KDI ), dimension ( : ), allocatable :: &
      Rank_AA
    real ( KDR ) :: &
      Pi
    character ( LDF ) :: &
      Name_AA
    type ( CommunicatorForm ), allocatable :: &
      Communicator_AA
    type ( GridImageStreamForm ), allocatable :: &
      GIS_AA
    type ( Atlas_SCG_Form ), allocatable :: &
      AA
    type ( Geometry_F_Form ), allocatable :: &
      GA
    type ( StreamForm ), allocatable :: &
      SA

    call Show ( 'Testing angular functions' )

    associate ( C  =>  A % Chart_GS )

    Name_AA  =  'AtlasAngular'

    allocate ( Rank_AA ( C % Communicator % Size  /  C % nBricks ( 1 ) ) )
    Rank_AA  =  [ ( iR  *  C % nBricks ( 1 ), &
                        iR = 0, size ( Rank_AA ) - 1 ) ]
 
    allocate ( Communicator_AA )
    call Communicator_AA % Initialize &
           ( C % Communicator, Rank_AA, Name_AA ) 

    if ( any ( C % Communicator % Rank == Rank_AA ) ) then

      Pi  =  CONSTANT % PI

      allocate ( GIS_AA )
      call GIS_AA % Initialize &
             ( Name_AA, CommunicatorOption = Communicator_AA )    

      allocate ( AA )
      call AA % Initialize &
             ( CommunicatorOption = Communicator_AA, &
               CoordinateLabelOption = [ 'r    ', 'Theta', 'Phi  ' ], &
               NameOption = Name_AA, &
               MinCoordinateOption = [ 0.0_KDR, 0.0_KDR, 0.0_KDR ], &
               MaxCoordinateOption = [ 1.0_KDR, Pi, 2.0_KDR * Pi ], &
               nCellsOption = [ 1, C % nCells ( 2 : 3 ) ], &
               nBricksOption = [ 1, C % nBricks ( 2 : 3 ) ], &
               nGhostLayersOption = [ 0, 0, 0 ], &  
               nDimensionsOption = 3 )

      allocate ( GA )
      call GA % Initialize ( AA, NameOption = 'GeometryAngular' )
    
      call AA % Show ( )
      call GA % Show ( )

      allocate ( SA )
      call SA % Initialize ( AA, GIS_AA )

      select case ( AA % Chart_GS % nDimensions )
      case ( 1 )
        call SA % CurveImage ( 1 ) % AddStorage ( L % AngularFunctions )
      case default
        call SA % GridImage ( 1 ) % AddStorage ( L % AngularFunctions )
      end select

      call GIS_AA % Open ( GIS_AA % ACCESS_CREATE )
      call SA % Write ( )
      call GIS_AA % Close ( )

      deallocate ( SA )
      deallocate ( GA )
      deallocate ( AA )
      deallocate ( GIS_AA )

    end if !-- Rank in Rank_AA

    deallocate ( Communicator_AA )

    end associate !-- C

  end subroutine TestAngularFunctions


  subroutine TestHomogeneousSphere ( )

    real ( KDR ) :: &
      RadiusDensity, &
      Density
    character ( LDL ), dimension ( 1 ) :: &
      Field
    type ( FieldSetForm ), allocatable :: &
      Source, &
      Reference
    type ( TimerForm ), pointer :: &
      T_CM, &
      T_W

    call Show ( 'Testing homogeneous sphere' )

    associate ( C  =>  A % Chart_GS )

    Field  =  [ 'HomogeneousSphere' ]

    allocate ( Source )
    call Source % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Source', &
             DeviceMemoryOption = L % DeviceMemory, &
             nFieldsOption = nEquations )
    
    allocate ( Reference )
    call Reference % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Reference', &
             nFieldsOption = nEquations )
    
    call S % AddFieldSet ( Source )
    call S % AddFieldSet ( Reference )
    call S % Show ( )

    RadiusDensity = C % MaxCoordinate ( 1 ) / 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    call SetHomogeneousSphere &
           ( Source, Reference, G, Density, RadiusDensity, iField = 1 )

    call Source % UpdateDevice ( )

    T_CM  =>  L % Timer ( Level = 1 )
    call T_CM % Start ( )
    call L % ComputeMoments ( Source, T_Option = T_CM )
    call T_CM % Stop ( )
    call L % ShowMoments ( )

    T_W  =>  S % TimerWrite ( Level = 1 )
    call T_W % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    call T_W % Stop ( )

    end associate !-- C

  end subroutine TestHomogeneousSphere


  subroutine SetHomogeneousSphere &
               ( Source, Reference, Geometry, &
                 Density, RadiusDensity, iField )

    class ( FieldSetForm ), intent ( inout ) :: &
      Source, &
      Reference
    class ( Geometry_F_Form ), intent ( in ) :: &
      Geometry
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity
    integer ( KDI ), intent ( in ) :: &
      iField

    !-- Geometry

    associate &
      ( GV  =>  Geometry % Storage_GS % Value )
    associate &
      ( R_E  =>  GV ( :, Geometry % EDGE_I_U ( 1 ) ), &
        R_W  =>  GV ( :, Geometry % WIDTH_U  ( 1 ) ), &
        R_C  =>  GV ( :, Geometry % CENTER_U ( 1 ) ) )

    !-- Source

    associate &
      ( SV  =>  Source % Storage_GS % Value )
    associate &
      ( D  =>  SV ( :, iField ) )

    call SetDensityKernel ( R_E, R_W, RadiusDensity, Density, D )

    end associate !-- D
    end associate !-- SV

    !-- Reference

    associate &
      ( RV  =>  Reference % Storage_GS % Value )
    associate &
      ( Phi  =>  RV ( :, iField ), &
        Pi   =>  CONSTANT % PI )

    where ( R_C  <  RadiusDensity )
      Phi  =  1.0_KDR / 6.0_KDR  *  Density  *  R_C ** 2  &
              -  1.0_KDR / 2.0_KDR  *  Density  *  RadiusDensity ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  Density  *  RadiusDensity ** 3  /  R_C
    end where

    end associate !-- Phi, etc.
    end associate !-- RV

    !-- Cleanup

    end associate !-- R_E, etc.
    end associate !-- GV

  end subroutine SetHomogeneousSphere


  subroutine SetDensityKernel ( R_E, R_W, RD, Density, D )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E, &
      R_W
    real ( KDR ), intent ( in ) :: &
      RD, &
      Density
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      D

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      R_I, R_O

    do iV  =  1, size ( R_E )
      R_I  =  R_E ( iV )
      R_O  =  R_E ( iV )  +  R_W ( iV )
      if ( R_O  <=  RD ) then
        D ( iV )  =  Density
      else if ( R_I  <  RD .and. R_O  >  RD ) then
        D ( iV )  =  Density * ( RD ** 3  -  R_I ** 3 ) &
                     / ( R_O ** 3  -  R_I ** 3 )
      else
        D ( iV )  =  0.0_KDR
      end if
    end do !-- iV

  end subroutine SetDensityKernel


end program Laplacian_M_ASCG__Form_Test
