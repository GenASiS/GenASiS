program Laplacian_M_A__Form_Test

  !-- Laplacian_Multipole_Atlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use PoissonEquations

  implicit none

  integer ( KDI ) :: &
    nEquations, &
    MaxDegree
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( Laplacian_M_A_Form ), allocatable :: &
    LA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Laplacian_M_A__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )
  call GA % SetStream ( SA )

  nEquations = 1

  MaxDegree = 3
  call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

  allocate ( LA )
  call LA % Initialize ( GA, MaxDegree, nEquations )

  call  A % Show ( )
  call GA % Show ( )
  call LA % Show ( )

  call TestAssociatedLegendre ( )
  call TestAngularFunctions ( )
  call TestHomogeneousSphere ( )

  deallocate ( LA )
  deallocate ( GA )
  deallocate ( SA )
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

    call Show ( LA % AssociatedLegendre ( X_Random, 0, 0 ), 'P_0_0 computed' )
    call Show ( sqrt ( 1.0_KDR / ( 4.0_KDR * Pi ) ), 'P_0_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 1, 0 ), 'P_1_0 computed' )
    call Show ( sqrt ( 3.0_KDR / ( 4.0_KDR * Pi ) ) * Cos_X, &
                'P_1_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 1, 1 ), 'P_1_1 computed' )
    call Show ( - sqrt ( 3.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X, &
                'P_1_1 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 0 ), 'P_2_0 computed' )
    call Show ( sqrt ( 5.0_KDR / ( 4.0_KDR * Pi ) ) &
                *  (    ( 3.0_KDR / 2.0_KDR )  *  Cos_X ** 2  &
                     -  ( 1.0_KDR / 2.0_KDR ) ), &
                'P_2_0 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 1 ), 'P_2_1 computed' )
    call Show ( - sqrt ( 15.0_KDR / ( 8.0_KDR * Pi ) ) * Sin_X * Cos_X, &
                'P_2_1 expected' )

    call Show ( LA % AssociatedLegendre ( X_Random, 2, 2 ), 'P_2_2 computed' )
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
    type ( Geometry_F_A_Form ), allocatable :: &
      GAA
    type ( Stream_A_Form ), allocatable :: &
      SAA

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

      allocate ( GAA )
      call GAA % Initialize ( AA, NameOption = 'GeometryAngular' )
    
      call  AA % Show ( )
      call GAA % Show ( )

      allocate ( SAA )
      call SAA % Initialize ( AA, GIS_AA )

      associate ( SC  =>  SAA % Stream_C ( 1 ) % Element )
      if ( allocated ( SC % CurveImage ) ) then
        call SC % CurveImage % AddStorage ( LA % AngularFunctions )
      else if ( allocated ( SC % GridImage ) ) then
        call SC % GridImage % AddStorage ( LA % AngularFunctions )
      end if
      end associate !-- SC

      call GIS_AA % Open ( GIS_AA % ACCESS_CREATE )
      call SAA % Write ( )
      call GIS_AA % Close ( )

      deallocate ( SAA )
      deallocate ( GAA )
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
    type ( FieldSet_A_Form ), allocatable :: &
      Source_A, &
      Reference_A

    call Show ( 'Testing homogeneous sphere' )

    associate ( C  =>  A % Chart_GS )

    Field  =  [ 'HomogeneousSphere' ]

    allocate ( Source_A )
    call Source_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Source', &
             DeviceMemoryOption = LA % DeviceMemory, &
             nFieldsOption = nEquations )
    
    allocate ( Reference_A )
    call Reference_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Reference', &
             nFieldsOption = nEquations )
    
    call SA % AddFieldSet ( Source_A )
    call SA % AddFieldSet ( Reference_A )
    call SA % Show ( )

    RadiusDensity = C % MaxCoordinate ( 1 ) / 10.0_KDR
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density = 1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    call SetHomogeneousSphere &
           ( Source_A, Reference_A, GA, Density, RadiusDensity, iField = 1 )

    call Source_A % UpdateDevice ( )

    call LA % ComputeMoments ( Source_A )
    call LA % ShowMoments ( )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

    end associate !-- C

  end subroutine TestHomogeneousSphere


  subroutine SetHomogeneousSphere &
               ( Source_A, Reference_A, Geometry_A, &
                 Density, RadiusDensity, iField )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Source_A, &
      Reference_A
    class ( Geometry_F_A_Form ), intent ( in ) :: &
      Geometry_A
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity
    integer ( KDI ), intent ( in ) :: &
      iField

    !-- Geometry

    select type ( GC  =>  Geometry_A % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )
    associate &
      ( GV  =>  GC % Storage_FSC % Storage % Value )
    associate &
      ( R_E  =>  GV ( :, GC % EDGE_I_U ( 1 ) ), &
        R_W  =>  GV ( :, GC % WIDTH_U  ( 1 ) ), &
        R_C  =>  GV ( :, GC % CENTER_U ( 1 ) ) )

    !-- Source

    associate &
      ( SC  =>  Source_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( SV  =>  SC % Storage_FSC % Storage % Value )
    associate &
      ( D  =>  SV ( :, iField ) )

    call SetDensityKernel ( R_E, R_W, RadiusDensity, Density, D )

    end associate !-- D
    end associate !-- SV
    end associate !-- SC

    !-- Reference

    associate &
      ( RC  =>  Reference_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( RV  =>  RC % Storage_FSC % Storage % Value )
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
    end associate !-- RC

    !-- Cleanup

    end associate !-- R_E, etc.
    end associate !-- GV
    end select !-- GC

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


end program Laplacian_M_A__Form_Test
