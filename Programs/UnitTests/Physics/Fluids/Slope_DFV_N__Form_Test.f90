program Slope_DFV_N__Form_Test

  !-- Slope_DivergenceFiniteVolume_Newton__Form_Test

  use Basics
  use Mathematics
  use Gravitations
  use Fluids

  implicit none

  integer ( KDI ) :: &
    nCompute
  logical ( KDL ) :: &
    DivergenceParts
  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_CC_Form ), allocatable :: &
    A
  type ( StreamForm ), allocatable :: &
    Sm
  type ( Gravitation_N_SG_Form ), allocatable :: &
    G
  type ( Units_F_Form ), dimension ( : ), allocatable :: &
    U
  type ( Fluid_D_Form ), allocatable :: &
    F
  class ( DivergencePart_CS_Form ), allocatable :: &
    DT
  type ( DivergencePartElement ), dimension ( : ), allocatable :: &
    DP_1D
  type ( RiemannSolver_HLL_Form ), allocatable :: &
    RS
  type ( Slope_DFV_N_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Slope_DFV_N__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( RadiusMax = 10.0_KDR, &
           RadiusCore = 10.0_KDR / 8.0_KDR, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( Sm )
  call Sm % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A, GravitationalConstant = 1.0_KDR )
  call G % SetStream ( Sm )
  call Sm % AddFieldSet ( G % Source )

  allocate ( U ( 1 ) )
  call U ( 1 ) % Initialize ( )

  allocate ( F )
  call F % Initialize ( G, U )
  call F % SetStream ( Sm )

  DivergenceParts  =  .false.
  call PROGRAM_HEADER % GetParameter ( DivergenceParts, 'DivergenceParts' )

  if ( DivergenceParts ) then
    allocate ( DP_1D ( 1 ) )
    allocate ( DivergencePart_F_D_V_Form :: DP_1D ( 1 ) % Element )
    associate ( DP  =>  DP_1D ( 1 ) % Element )
    call DP % Initialize ( F )
    end associate !-- DP
  else
    allocate ( DivergencePart_F_D_T_Form :: DT )
    call DT % Initialize ( F )
  end if

  allocate ( RS )
  call RS % Initialize ( F )

  allocate ( S )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  if ( DivergenceParts ) then
    call S % Initialize &
           ( RS, DP_1D, &
             iVelocity_F = F % VELOCITY_U, &
             iMomentum_B = [ 2, 3, 4 ], &
             iBaryonMass_F = F % BARYON_MASS, &
             iBaryonDensity_F = F % BARYON_DENSITY_C, &
             iEnergy_B = 0 )
  else
    call S % Initialize &
           ( RS, DT, &
             iVelocity_F = F % VELOCITY_U, &
             iMomentum_B = [ 2, 3, 4 ], &
             iBaryonMass_F = F % BARYON_MASS, &
             iBaryonDensity_F = F % BARYON_DENSITY_C, &
             iEnergy_B = 0 )
  end if
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call S % SetStream ( Sm )

  call A % Show ( )
  call G % Show ( )
  call F % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
  call S % Show ( )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  call Sm % Show ( )

  nCompute  =  1
  call PROGRAM_HEADER % GetParameter ( nCompute, 'nCompute' )

  call TestSlope ( )

  call CONSOLE % SetVerbosity ( 'INFO_2' )
  deallocate ( S )
  call CONSOLE % SetVerbosity ( 'INFO_1' )
  deallocate ( RS )
  if ( allocated ( DP_1D ) ) &
    deallocate ( DP_1D )
  if ( allocated ( DT ) ) &
    deallocate ( DT )
  deallocate ( F )
  deallocate ( G )
  deallocate ( Sm )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine TestSlope ( )

    integer ( KDI ) :: &
      iHS
    real ( KDR ), dimension ( 3 ) :: &
      Radius, &
      Density

    call Show ( 'Testing homogeneous spheres' )

    associate ( C  =>  A % Chart_GS )
    
    Radius  =  C % MaxCoordinate ( 1 ) / [ 1.1_KDR, 2.0_KDR, 10.0_KDR ]
    call PROGRAM_HEADER % GetParameter ( Radius, 'Radius' )

    Density  =  1.0_KDR  /  Radius ** 3
    call PROGRAM_HEADER % GetParameter ( Density, 'Density' )

    do iHS  =  1, 3

      call SetHomogeneousSphere &
             ( F, G, Density ( iHS ), Radius ( iHS )  )
    
      call CONSOLE % SetVerbosity ( 'INFO_7' )
      call F % UpdateDevice ( )
      call F % ComputeFromPrimitive ( F )
      call G % Solve &
             ( F, &
               iBaryonMass = F % BARYON_MASS, &
               iBaryonDensity = F % BARYON_DENSITY_B )
      call S % Compute ( )
      call CONSOLE % SetVerbosity ( 'INFO_1' )

      call Show ( Radius ( iHS ), 'Radius', nLeadingLinesOption = 2 )
      call Show ( Density ( iHS ), 'Density' )

      call GIS % Open ( GIS % ACCESS_CREATE )
      call Sm % Write ( )
      call GIS % Close ( )

    end do !-- iHS

    end associate !-- C

  end subroutine TestSlope


  subroutine SetHomogeneousSphere &
               ( Fluid, Geometry, Density, Radius )

    class ( Fluid_D_Form ), intent ( inout ) :: &
      Fluid
    class ( Geometry_F_Form ), intent ( in ) :: &
      Geometry
    real ( KDR ), intent ( in ) :: &
      Density, &
      Radius

    !-- Geometry

    associate &
      ( GV  =>  Geometry % Storage_GS % Value )
    associate &
      ( R_E  =>  GV ( :, G % EDGE_I_U ( 1 ) ), &
        R_W  =>  GV ( :, G % WIDTH_U  ( 1 ) ), &
        R_C  =>  GV ( :, G % CENTER_U ( 1 ) ) )

    !-- Fluid

    associate &
      ( FV  =>  Fluid % Storage_GS % Value )
    associate &
      ( M  =>  FV ( :, F % BARYON_MASS ), &
        D  =>  FV ( :, F % BARYON_DENSITY_C ) )

    call SetDensityKernel ( R_E, R_W, Radius, Density, M, D )

    end associate !-- M, etc.
    end associate !-- FV

    !-- Cleanup

    end associate !-- R_E, etc.
    end associate !-- GV

  end subroutine SetHomogeneousSphere


  subroutine SetDensityKernel ( R_E, R_W, RD, Density, M, D )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_E, &
      R_W
    real ( KDR ), intent ( in ) :: &
      RD, &
      Density
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      M, &
      D

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      R_I, R_O

    do iV  =  1, size ( R_E )

      M ( iV )  =  1.0_KDR

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


end program Slope_DFV_N__Form_Test
