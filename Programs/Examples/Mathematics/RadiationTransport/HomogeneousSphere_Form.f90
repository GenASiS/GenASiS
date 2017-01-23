module HomogeneousSphere_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_C__Form
  use Interactions_CSL__Form
  use Interactions_ASC__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_Template ) :: HomogeneousSphereForm
    type ( Interactions_ASC_Form ), allocatable :: &
      Interactions_ASC
    type ( RadiationMoments_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type HomogeneousSphereForm

    private :: &
      SetRadiation

      private :: &
        SetRadiationKernel

      real ( KDR ), private :: &
      EquilibriumDensity, &
      EffectiveOpacity, &
      TransportOpacity

contains


  subroutine Initialize ( HS, Name )

    class ( HomogeneousSphereForm ), intent ( inout ) :: &
      HS
    character ( * ), intent ( in )  :: &
      Name

    class ( RadiationMomentsForm ), pointer :: &
      RM
    integer ( KDI ) :: &
      nCellsRadius
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    real ( KDR ) :: &
      MaxRadius
    character ( LDL ) :: &
      CoordinateSystem

    if ( HS % Type == '' ) &
      HS % Type = 'a HomogeneousSphere'

    EquilibriumDensity = 1.0_KDR
    EffectiveOpacity   = 5.0_KDR
    TransportOpacity   = EffectiveOpacity

    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: HS % PositionSpace )
    select type ( PS => HS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    MaxRadius = 7.0_KDR
    call PROGRAM_HEADER % GetParameter ( MaxRadius, 'MaxRadius' )

    CoordinateSystem = 'CARTESIAN'
    call PROGRAM_HEADER % GetParameter ( CoordinateSystem, 'CoordinateSystem' )

    nCellsRadius = 64
    call PROGRAM_HEADER % GetParameter ( nCellsRadius, 'nCellsRadius' )

    select case ( CoordinateSystem )
      case ( 'CARTESIAN' )
         MinCoordinate = [ - MaxRadius, - MaxRadius, - MaxRadius ]
         MaxCoordinate = [ + MaxRadius, + MaxRadius, + MaxRadius]
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 1 )
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 3 )
         nCells = [ 32, 32, 32 ]
        call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
      case ( 'SPHERICAL' )
         associate ( Pi => CONSTANT % PI )
         MinCoordinate = [   0.0_KDR, 0.0_KDR, 0.0_KDR ]
         MaxCoordinate = [ MaxRadius,      Pi, 2.0_KDR * Pi ]
         end associate !-- Pi

         call PS % SetBoundaryConditionsFace &
                ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
         nCells = [ nCellsRadius, 1, 1 ]
         call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

         if ( PS % nDimensions > 1 ) then
           call PS % SetBoundaryConditionsFace &
                  ( [ 'REFLECTING', 'REFLECTING' ], iDimension = 2 )
           nCells = [ nCellsRadius, nCellsRadius, 1 ]
           call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
         end if

         if ( PS % nDimensions > 2 ) then
           call PS % SetBoundaryConditionsFace &
                  ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )
           nCells = [ 32, 32, 32 ]
           call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
         end if
      case ( 'CYLINDRICAL' )
         associate ( Pi => CONSTANT % PI )
         MinCoordinate = [   0.0_KDR, - MaxRadius, 0.0_KDR ]
         MaxCoordinate = [ MaxRadius, + MaxRadius, 2.0_KDR * Pi ]
         end associate !-- Pi

         call PS % SetBoundaryConditionsFace &
                ( [ 'REFLECTING', 'OUTFLOW   ' ], iDimension = 1 )
         nCells = [ nCellsRadius, 1, 1 ]
         call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

         if ( PS % nDimensions > 1 ) then
           call PS % SetBoundaryConditionsFace &
                  ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
           nCells = [ nCellsRadius, 2 * nCellsRadius, 1 ]
           call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
         end if

         if ( PS % nDimensions > 2 ) then
           call PS % SetBoundaryConditionsFace &
                  ( [ 'PERIODIC', 'PERIODIC' ], iDimension = 3 )
           nCells = [ 32, 32, 32 ]
           call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
         end if    
      case DEFAULT 
        call show ( CoordinateSystem, 'CoordinateSystem' )
        call Show ( 'Coordinate System not implemented', CONSOLE % ERROR )
        call Show ( 'HomogeneousSphere_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    call PS % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells )

    !-- Geometry of PositionSpace

    allocate ( HS % Geometry_ASC )
    associate ( GA => HS % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Interactions
    allocate ( HS % Interactions_ASC )
    associate ( IA => HS % Interactions_ASC )
    
    call IA % Initialize &
               ( PS, EquilibriumDensityOption = EquilibriumDensity, &
                 EffectiveOpacityOption = EffectiveOpacity, &
                 TransportOpacityOption = TransportOpacity )
    
    select type ( IC => IA % Chart )
    class is ( Interactions_CSL_Form )
    select type ( I => IC % Field )
    class is ( Interactions_C_Form )
    associate &
      ( ED => I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
        EO => I % Value ( :, I % EFFECTIVE_OPACITY ), &
        TO => I % Value ( :, I % TRANSPORT_OPACITY ) )

    call SetInteractions &
           ( HS, CoordinateSystem, EquilibriumDensity, EffectiveOpacity, &
             TransportOpacity, RM, ED, EO, TO )
    end associate !-- ED, etc.
    end select !-- I
    end select !-- IC

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: HS % Current_ASC )
    select type ( RMA => HS % Current_ASC )  !-- FluidAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )
    call RMA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_Form :: HS % Step )
    select type ( S => HS % Step )
    class is ( Step_RK2_C_Form )
    call S % Initialize ( Name )
    S % ApplyRelaxation  =>  ApplyRelaxation_Interactions
    S % ApplySources  => ApplySourcesCurvilinear_RadiationMoments
    end select !-- S
    
    !-- Initial conditions
    
    RM => RMA % RadiationMoments ( )
    call SetRadiation ( HS, RM )
    
    !-- Initialize template
    
    call HS % InitializeTemplate_C ( Name, FinishTimeOption = 25.0_KDR )
    
    !-- Cleanup
           
    end select !-- RMA
    end associate !-- IA
    end select !-- PS
    nullify ( RM )
    
  end subroutine Initialize 


  impure elemental subroutine Finalize ( HS )

    type ( HomogeneousSphereForm ), intent ( inout ) :: &
      HS

    if ( allocated ( HS % Difference ) ) &
      deallocate ( HS % Difference )
    if ( allocated ( HS % Reference ) ) &
      deallocate ( HS % Reference )

    call HS % FinalizeTemplate_C ( )

  end subroutine Finalize


  subroutine SetRadiation ( HS, RM )

    class ( HomogeneousSphereForm ), intent ( in ) :: &
      HS
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    
    class ( GeometryFlatForm ), pointer :: &
      G

    select type ( PS => HS % PositionSpace )
      class is ( Atlas_SC_Form )
      G => PS % Geometry ( )

      call SetRadiationKernel &
             ( J  = RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
               HX = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
               HY = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
               HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ) )
 
      call RM % ComputeFromPrimitive ( G )
      end select !-- PS

  end subroutine SetRadiation


  subroutine SetInteractions &
               ( HS, CoordinateSystem, EquilibriumDensity, EffectiveOpacity, &
                 TransportOpacity, RM, ED, EO, TO )
    class ( HomogeneousSphereForm ), intent ( in ) :: &
      HS
    character ( LDL ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), intent ( in ) :: &
      EquilibriumDensity,&
      EffectiveOpacity, &
      TransportOpacity
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      ED, &
      EO, &
      TO

    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      R
    integer ( KDI ) :: &
     iV, & !-- iValue
     nValues
    
    select type ( PS => HS % PositionSpace )
      class is ( Atlas_SC_Form )
      G => PS % Geometry ( )
      
      nValues = size ( G % Value ( :, G % CENTER ( 1 ) ) )
      allocate ( R ( nValues ) )

      select case ( CoordinateSystem ) 
        case ( 'CARTESIAN' )
          associate &
            ( X =>  G % Value ( :, G % CENTER ( 1 ) ), &
              Y =>  G % Value ( :, G % CENTER ( 2 ) ), &
              Z =>  G % Value ( :, G % CENTER ( 3 ) ) )
  
          !$OMP parallel do private ( iV )
          do iV = 1, nValues
            R ( iV ) = sqrt ( X ( iV ) ** 2 + Y ( iV ) ** 2 + Z ( iV ) ** 2 )
          end do !-- iV
          !$OMP end parallel do
          end associate !-- X, etc.
        case ( 'CYLINDRICAL' )
          associate &
            ( Rho =>  G % Value ( :, G % CENTER ( 1 ) ), &
              Z   =>  G % Value ( :, G % CENTER ( 2 ) ) )
          !$OMP parallel do private ( iV )
          do iV = 1, nValues
            R ( iV ) = sqrt ( ( Rho ( iV ) ** 2 + Z ( iV ) ** 2 ) )
          end do !-- iV
          !$OMP end parallel do
          end associate !-- Rho, etc. 
        case ( 'SPHERICAL' )
          !$OMP parallel do private ( iV )
          do iV = 1, nValues
            R ( iV ) = G % Value ( iV, G % CENTER ( 1 ) )
          end do !-- iV
          !$OMP end parallel do
        case DEFAULT
          call show ( CoordinateSystem, 'CoordinateSystem' )
          call Show ( 'Coordinate System not implemented', CONSOLE % ERROR )
          call Show ( 'HomogeneousSphere_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetInteractions', 'subroutine', CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- CoordinateSystem

        !$OMP parallel do private ( iV )
        do iV = 1, nValues
          if ( R ( iV ) < 1.0_KDR ) then
            ED ( iV ) = EquilibriumDensity
            EO ( iV ) = EffectiveOpacity 
            TO ( iV ) = TransportOpacity 
          else 
            ED ( iV ) = 0.0_KDR
            EO ( iV ) = 0.0_KDR
            TO ( iV ) = 0.0_KDR 
          end if
        end do !-- iV
        !$OMP end parallel do
      end select !-- PS

  end subroutine SetInteractions


  subroutine SetRadiationKernel ( J, HX, HY, HZ )

    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    

    nValues = size ( J )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      J ( iV )   =  0.0_KDR
      HX ( iV )  =  0.0_KDR
      HY ( iV )  =  0.0_KDR
      HZ ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetRadiationKernel


end module HomogeneousSphere_Form
