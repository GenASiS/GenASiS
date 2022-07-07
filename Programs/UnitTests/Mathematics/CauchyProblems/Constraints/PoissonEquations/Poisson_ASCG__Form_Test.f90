program Poisson_ASCG__Form_Test

  !-- Poisson_AtlasSingleChartGrid__Form_Test

  use Basics
  use Algebra
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
  type ( Poisson_ASCG_Form ), allocatable :: &
    PA

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Poisson_ASCG__Form_Test', DimensionalityOption = '2D' )

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

  nEquations = 3

  MaxDegree = 3
  call PROGRAM_HEADER % GetParameter ( MaxDegree, 'MaxDegree' )

  allocate ( PA )
  call PA % Initialize ( G, 'MULTIPOLE', MaxDegree, nEquations )

  call  A % Show ( )
  call  G % Show ( )
  call PA % Show ( )

  call TestHomogeneousSpheres ( )

  deallocate ( PA )
  deallocate ( G )
  deallocate ( S )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestHomogeneousSpheres ( )

    integer ( KDI ) :: &
      iE, &  !-- iEquation
      iS, &  !-- iSolve
      nSolve
    real ( KDR ), dimension ( nEquations ) :: &
      RadiusDensity, &
      Density
    character ( LDL ), dimension ( nEquations ) :: &
      Field
    type ( FieldSetForm ), allocatable :: &
      Source, &
      Solution, &
      Reference, &
      Difference
    type ( GradientForm ), allocatable :: &
      Gradient
    type ( TimerForm ), pointer :: &
      T_P, &
      T_W

    call Show ( 'Testing homogeneous spheres' )

    associate ( C  =>  A % Chart_GS )

    Field  =  [ 'HomogeneousSphere_1', &
                'HomogeneousSphere_2', &
                'HomogeneousSphere_3' ]

    allocate ( Source )
    call Source % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Source', &
             DeviceMemoryOption = PA % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    
    allocate ( Solution )
    call Solution % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Solution', &
             DeviceMemoryOption = PA % Laplacian_M % DeviceMemory, &
             nFieldsOption = nEquations )
    call Solution % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    call Solution % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    
    allocate ( Reference )
    call Reference % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Reference', &
             nFieldsOption = nEquations )
    
    allocate ( Difference )
    call Difference % Initialize &
           ( A, &
             FieldOption = Field, &
             NameOption = 'Difference', &
             nFieldsOption = nEquations )
    
    allocate ( Gradient )
    call Gradient % Initialize ( G, Solution )

    call S % AddFieldSet ( Source )
    call S % AddFieldSet ( Solution )
    call S % AddFieldSet ( Reference )
    call S % AddFieldSet ( Difference )
    call S % AddFieldSet ( Gradient )
    call S % Show ( )

    RadiusDensity  =  C % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ]
    call PROGRAM_HEADER % GetParameter ( RadiusDensity, 'RadiusDensity' )

    Density  =  1.0_KDR  /  RadiusDensity ** 3

    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iE  =  1, nEquations
      call SetHomogeneousSphere &
             ( Source, Reference, G, &
               Density ( iE ), RadiusDensity ( iE ), iField = iE )
    end do !-- iE
    
    call Source % UpdateDevice ( )

    nSolve  =  10
    call PROGRAM_HEADER % GetParameter ( nSolve, 'nSolve' )
    call Show ( 'Poisson solve' )
    call Show ( nSolve, 'nSolve' )

    T_P  =>  PA % Timer ( Level = 1 )
    call T_P % Start ( )
    do iS  =  1,  nSolve
      call PA % Solve ( Solution, Source, T_Option = T_P )
    end do !-- iS
    call T_P % Stop ( )

    call ComputeError ( Difference, Solution, Reference )

    call Gradient % Compute ( iD = 1 )

    T_W  =>  S % TimerWrite ( Level = 1 )
    call T_W % Start ( )
    call GIS % Open ( GIS % ACCESS_CREATE )
    call S % Write ( )
    call GIS % Close ( )
    call T_W % Stop ( )

    end associate !-- C

  end subroutine TestHomogeneousSpheres


  subroutine ComputeError ( Difference, Solution, Reference )

    class ( FieldSetForm ), intent ( inout ) :: &
      Difference, &         
      Solution, &
      Reference

    real ( KDR ) :: &
      L1_1, &
      L1_2, &
      L1_3
    type ( CollectiveOperation_R_Form ) :: &
      CO
    
    associate &
      ( SV  =>  Solution % Storage_GS % Value, &
        RV  =>  Reference % Storage_GS % Value, &
        DV  =>  Difference % Storage_GS % Value )

    call MultiplyAdd ( SV, RV, -1.0_KDR, DV )

    select type ( A  =>  Solution % Atlas )
     class is ( Atlas_SCG_Form ) 
    associate &
      ( C  =>  A % Chart_GS )

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

    end associate !-- C
    end select !-- A
    end associate !-- SV, etc.

  end subroutine ComputeError


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
      ( D  =>  SV ( :, iField ), &
        FourPi  =>  4.0_KDR  *  CONSTANT % PI )

    call SetDensityKernel ( R_E, R_W, RadiusDensity, Density, D )
    D  =  FourPi  *  D

    end associate !-- D
    end associate !-- SV

    !-- Reference

    associate &
      ( RV  =>  Reference % Storage_GS % Value )
    associate &
      ( Phi  =>  RV ( :, iField ), &
        FourPi   =>  4.0_KDR  *  CONSTANT % PI )

    where ( R_C  <  RadiusDensity )
      Phi  =  1.0_KDR / 6.0_KDR  *  FourPi  *  Density  *  R_C ** 2  &
              -  1.0_KDR / 2.0_KDR  *  FourPi  *  Density  &
                                    *  RadiusDensity ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  FourPi  *  Density  &
                                   *  RadiusDensity ** 3  /  R_C
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

    do iV  =  1, size ( D )
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


end program Poisson_ASCG__Form_Test
