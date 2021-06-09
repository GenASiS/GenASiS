program Poisson_A__Form_Test

  !-- Poisson_Atlas__Form_Test

  use Basics
  use Algebra
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
  type ( Poisson_A_Form ), allocatable :: &
    PA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Poisson_A__Form_Test', DimensionalityOption = '2D' )

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

  nEquations = 3

  MaxDegree = 3
  call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

  allocate ( PA )
  call PA % Initialize ( GA, 'MULTIPOLE', MaxDegree, nEquations )

  call  A % Show ( )
  call GA % Show ( )
  call PA % Show ( )

  call TestHomogeneousSpheres ( )

  deallocate ( PA )
  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestHomogeneousSpheres ( )

    integer ( KDI ) :: &
      iE  !-- iEquation
    real ( KDR ), dimension ( nEquations ) :: &
      RadiusDensity, &
      Density
    character ( LDL ), dimension ( nEquations ) :: &
      Field
    type ( FieldSet_A_Form ), allocatable :: &
      Source_A, &
      Solution_A, &
      Reference_A, &
      Difference_A
    type ( Gradient_A_Form ), allocatable :: &
      Gradient_A

    call Show ( 'Testing homogeneous spheres' )

    associate ( C  =>  A % Chart_GS )

    Field  =  [ 'HomogeneousSphere_1', &
                'HomogeneousSphere_2', &
                'HomogeneousSphere_3' ]

    allocate ( Source_A )
    call Source_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Source', &
             DeviceMemoryOption = PA % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    
    allocate ( Solution_A )
    call Solution_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Solution', &
             DeviceMemoryOption = PA % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    call Solution_A % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
    call Solution_A % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
    
    allocate ( Reference_A )
    call Reference_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Reference', &
             nFieldsOption = nEquations )
    
    allocate ( Difference_A )
    call Difference_A % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Difference', &
             nFieldsOption = nEquations )
    
    allocate ( Gradient_A )
    call Gradient_A % Initialize ( GA, Solution_A )

    call SA % AddFieldSet ( Source_A )
    call SA % AddFieldSet ( Solution_A )
    call SA % AddFieldSet ( Reference_A )
    call SA % AddFieldSet ( Difference_A )
    call SA % AddFieldSet ( Gradient_A )
    call SA % Show ( )

    RadiusDensity  =  C % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ]
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density  =  1.0_KDR / ( 4.0_KDR  *  CONSTANT % PI  *  RadiusDensity ** 3 )

    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE  =  1, nEquations
      call SetHomogeneousSphere &
             ( Source_A, Reference_A, GA, &
               Density ( iE ), RadiusDensity ( iE ), iField = iE )
    end do !-- iE
    
    call Source_A % UpdateDevice ( )

    call PA % Solve ( Solution_A, Source_A )
    call ComputeError ( Difference_A, Solution_A, Reference_A )

    call Gradient_A % Compute ( iD = 1 )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

    end associate !-- C

  end subroutine TestHomogeneousSpheres


  subroutine ComputeError ( Difference_A, Solution_A, Reference_A )

    class ( FieldSet_A_Form ), intent ( inout ) :: &
      Difference_A, &         
      Solution_A, &
      Reference_A

    real ( KDR ) :: &
      L1_1, &
      L1_2, &
      L1_3
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    associate &
      ( SC  =>  Solution_A   % FieldSet_C ( 1 ) % Element, &
        RC  =>  Reference_A  % FieldSet_C ( 1 ) % Element, &
        DC  =>  Difference_A % FieldSet_C ( 1 ) % Element )
    associate &
      ( SV  =>  SC % Storage_FSC % Storage % Value, &
        RV  =>  RC % Storage_FSC % Storage % Value, &
        DV  =>  DC % Storage_FSC % Storage % Value )

    call MultiplyAdd ( SV, RV, -1.0_KDR, DV )

    select type ( C  =>  SC % Chart )
    class is ( Chart_GS_Form ) 

    call CO % Initialize &
           ( C % Communicator, [ 2 * nEquations ], [ 2 * nEquations ] )
    CO % Outgoing % Value ( 1 )  &
      =  sum ( abs ( DV ( :, 1 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 2 )  &
      =  sum ( abs ( DV ( :, 2 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 3 )  &
      =  sum ( abs ( DV ( :, 3 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 4 )  &
      =  sum ( abs ( RV ( :, 1 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 5 )  &
      =  sum ( abs ( RV ( :, 2 ) ), mask = C % ProperCell )
    CO % Outgoing % Value ( 6 )  &
      =  sum ( abs ( RV ( :, 3 ) ), mask = C % ProperCell )
    call CO % Reduce ( REDUCTION % SUM )

    associate &
      ( Norm_D_1  =>  CO % Incoming % Value ( 1 ), &
        Norm_D_2  =>  CO % Incoming % Value ( 2 ), &
        Norm_D_3  =>  CO % Incoming % Value ( 3 ), &
        Norm_R_1  =>  CO % Incoming % Value ( 4 ), &
        Norm_R_2  =>  CO % Incoming % Value ( 5 ), &
        Norm_R_3  =>  CO % Incoming % Value ( 6 ) )

    L1_1  =  Norm_D_1 / Norm_R_1
    L1_2  =  Norm_D_2 / Norm_R_2
    L1_3  =  Norm_D_3 / Norm_R_3

    end associate !-- Norm_D_1, etc.

    call Show ( L1_1, '*** L1_1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_2, '*** L1_2 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_3, '*** L1_3 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    ! Difference % Value = abs ( Difference % Value / Reference % Value )

    end select !-- C
    end associate !-- SV, etc.
    end associate !-- SC, etc.

  end subroutine ComputeError


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


end program Poisson_A__Form_Test
