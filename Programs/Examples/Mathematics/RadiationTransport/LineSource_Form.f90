module LineSource_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_ASC__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form

  implicit none
  private

  type, public, extends ( Integrator_C_Template ) :: LineSourceForm
    type ( Interactions_ASC_Form ), allocatable :: &
      Interactions_ASC
    type ( RadiationMoments_ASC_Form ), allocatable :: &
      Reference, &
      Difference
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type LineSourceForm

    private :: &
      SetRadiation, &
      SetReference

      private :: &
        SetRadiationKernel

      real ( KDR ), private :: &
      EquilibriumDensity, &
      EffectiveOpacity, &
      TransportOpacity

contains


  subroutine Initialize ( LS, Name )

    class ( LineSourceForm ), intent ( inout ) :: &
      LS
    character ( * ), intent ( in )  :: &
      Name

    class ( RadiationMomentsForm ), pointer :: &
      RM
    integer ( KDI ), dimension ( 3 ) :: &
      nCells
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate
    real ( KDR ) :: &
      MaxRadius
    character ( LDL ) :: &
      CoordinateSystem

    if ( LS % Type == '' ) &
      LS % Type = 'a LineSource'

    EquilibriumDensity = 0.0_KDR
    EffectiveOpacity   = 0.0_KDR
    TransportOpacity   = 0.0_KDR

    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: LS % PositionSpace )
    select type ( PS => LS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    MaxRadius = 1.5_KDR
    call PROGRAM_HEADER % GetParameter ( MaxRadius, 'MaxRadius' )

    CoordinateSystem = 'CYLINDRICAL'
    call PROGRAM_HEADER % GetParameter ( CoordinateSystem, 'CoordinateSystem' )

    select case ( CoordinateSystem )
      case ( 'CYLINDRICAL' )
        associate ( Pi => CONSTANT % PI )
          MinCoordinate = [ - MaxRadius, -1.0_KDR, 0.0_KDR ]
          MaxCoordinate = [ + MaxRadius, + 1.0_KDR, 2.0_KDR * Pi]
         ! MinCoordinate = [   0.0_KDR, 0.0_KDR, - MaxRadius ]
         ! MaxCoordinate = [ MaxRadius, 2.0_KDR * Pi, + MaxRadius]
        end associate !-- Pi
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 1 )
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
        
        nCells = [ 128, 128, 1 ]
        call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
      case ( 'CARTESIAN' )
       ! MinCoordinate = [ - MaxRadius, - MaxRadius, 0.0_KDR ]
       ! MaxCoordinate = [ + MaxRadius, + MaxRadius, 0.0_KDR ]
         MinCoordinate = [   0.0_KDR, 0.0_KDR, - MaxRadius ]
         MaxCoordinate = [ MaxRadius, 2.0_KDR * CONSTANT % PI, + MaxRadius]
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 1 )
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
         nCells = [ 128, 128, 1 ]
        call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
      case DEFAULT 
        call show ( CoordinateSystem, 'CoordinateSystem' )
        call Show ( 'Coordinate System not implemented', CONSOLE % ERROR )
        call Show ( 'LineSource_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    call PS % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells )

    !-- Geometry of PositionSpace

    allocate ( LS % Geometry_ASC )
    associate ( GA => LS % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Interactions
    allocate ( LS % Interactions_ASC )
    associate ( IA => LS % Interactions_ASC )
    call IA % Initialize &
               ( PS, EquilibriumDensityOption = EquilibriumDensity, &
                 EffectiveOpacityOption = EffectiveOpacity, &
                 TransportOpacityOption = TransportOpacity )

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: LS % Current_ASC )
    select type ( RMA => LS % Current_ASC )  !-- FluidAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )
    call RMA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_Form :: LS % Step )
    select type ( S => LS % Step )
    class is ( Step_RK2_C_Form )
    call S % Initialize ( Name )
    S % ApplyRelaxation  =>  ApplyRelaxation_Interactions
    S % ApplySources  => ApplySourcesCurvilinear_RadiationMoments
    end select !-- S

    !-- Diagnostics

    allocate ( LS % Reference )
    allocate ( LS % Difference )
    call LS % Reference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Reference' )
    call LS % Difference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Difference' )
    LS % SetReference => SetReference
    
    !-- Initial conditions
    
    RM => RMA % RadiationMoments ( )
    call SetRadiation &
         ( LS, RM, Time = 0.0_KDR, isInitial=.true. )
    
    !-- Initialize template
    
    call LS % InitializeTemplate_C ( Name )
    
    !-- Cleanup
           
    end select !-- RMA
    end associate !-- IA
    end select !-- PS
    nullify ( RM )
    
  end subroutine Initialize 

  subroutine ComputeError ( LS )

    class ( LineSourceForm ), intent ( in ) :: &
      LS
    
    real ( KDR ) :: &
      L1       
    class ( RadiationMomentsForm ), pointer :: &
      RM         
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( PS => LS % PositionSpace )
    class is ( Atlas_SC_Form ) 
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
    
    RM => LS % Difference % RadiationMoments ( )
    
    associate &
      ( Difference => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ) )
    call CO % Initialize ( PS % Communicator, [ 2 ], [ 2 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )
    end associate !-- Difference
    
    associate &
      ( DifferenceSum => CO % Incoming % Value ( 1 ), &
        nValues => CO % Incoming % Value ( 2 ) )
    L1 = DifferenceSum / nValues
    end associate

    call Show ( L1, '*** L1 error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end select !-- PS

  end subroutine ComputeError


  impure elemental subroutine Finalize ( LS )

    type ( LineSourceForm ), intent ( inout ) :: &
      LS

    if ( allocated ( LS % Difference ) ) &
      deallocate ( LS % Difference )
    if ( allocated ( LS % Reference ) ) &
      deallocate ( LS % Reference )

    call LS % FinalizeTemplate_C ( )

  end subroutine Finalize


  subroutine SetRadiation ( LS, RM, Time, isInitial )

    class ( LineSourceForm ), intent ( in ) :: &
      LS
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM
    real ( KDR ), intent ( in ) :: &
      Time
    logical ( KDL ), intent ( in ) :: &
      isInitial

    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      R
    character ( LDL ) :: &
      CoordinateSystem
    
    select type ( PS => LS % PositionSpace )
      class is ( Atlas_SC_Form )
      G => PS % Geometry ( )
      
      allocate ( R ( size ( G % Value ( :, G % CENTER ( 1 ) ) ) ) )

      associate &
        ( X =>  G % Value ( :, G % CENTER ( 1 ) ), &
          Z =>  G % Value ( :, G % CENTER ( 2 ) ) )

      R = sqrt ( X ** 2 + Z ** 2 )
        
      call SetRadiationKernel &
             ( X, &
               Time, &
               isInitial, &
               J  = RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
               HX = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 1 ) ), &
               HY = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 2 ) ), &
               HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_DENSITY_U ( 3 ) ) )
    

      call RM % ComputeFromPrimitive ( G )
   
      end associate !-- X, etc.
    end select    !-- PS
    
    nullify ( G )
    deallocate ( R )

  end subroutine SetRadiation


  subroutine SetReference ( LS )

    class ( IntegratorTemplate ), intent ( in ) :: &
      LS

    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( LS )
    class is ( LineSourceForm )

    select type ( RMA => LS % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- RMA

    RM_R => LS % Reference % RadiationMoments ( )
    call SetRadiation ( LS, RM_R, LS % Time, isInitial=.false. )

    RM_D => LS % Difference % RadiationMoments ( )
!    RM_D % Value  =  RM % Value  -  RM_R % Value
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- LS
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


  subroutine SetRadiationKernel ( R, T, isInitial, J, HX, HY, HZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R    
    real ( KDR ), intent ( in ) :: &
      T
    logical ( KDL ), intent ( in ) :: &
      isInitial
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      w
    

    nValues = size ( R )

    w = 0.03_KDR 

    if ( isInitial ) then
      !$OMP parallel do private ( iV )
      do iV = 1, nValues
        J ( iV )  =  max ( 1.0_KDR / ( 2.0_KDR * w ** 2 ) * &
                       exp ( - R ( iV ) ** 2 / ( 2 * w ** 2 ) ), &
                          10.0_KDR **(-4)  )

        HX ( iV )  =  0.0_KDR
        HY ( iV )  =  0.0_KDR
        HZ ( iV )  =  0.0_KDR

      end do !-- iV
      !$OMP end parallel do
    else 
      !$OMP parallel do private ( iV )
      do iV = 1, nValues
        J ( iV )  = 0.0_KDR
         if ( ( T ** 2 - R ( iV ) ** 2 ) > 0.0_KDR ) then
          J ( iV ) = 2.0_KDR / ( T * sqrt ( T ** 2 - R ( iV ) ** 2 ) )
        end if
      end do !--iV
      !$OMP end parallel do
   end if

  end subroutine SetRadiationKernel


end module LineSource_Form
