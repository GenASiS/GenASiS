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
    procedure, public, pass :: &
      ComputeError
    final :: &
      Finalize
  end type HomogeneousSphereForm

    private :: &
      SetRadiation, &
      SetReference

      private :: &
        SetRadiationKernel

      integer ( KDI ), private, parameter :: &
      iRADIUS_TS                     = 1, &  !-- must match the profile file columns
      iCOMOVING_ENERGY_DENSITY_TS    = 2, &
      iFLUX_FACTOR_TS                = 3, &
      iVARIABLE_EDDINGTON_FACTOR_TS  = 4
    integer ( KDI ), private, parameter :: &
      iCOMOVING_ENERGY_DENSITY_SI      = 1, & !-- spline interpolation
      iFLUX_FACTOR_SI                  = 2, & 
      iVARIABLE_EDDINGTON_FACTOR_SI    = 3
      
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
    class ( GeometryFlatForm ), pointer :: &
      G
    type ( SplineInterpolationForm ), dimension ( 3 ) :: &
      SI
    integer ( KDI ) :: &
      nCellsRadius, &
      iV
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

    EquilibriumDensity = 10.0_KDR
    EffectiveOpacity   = 250.0_KDR

    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )

    TransportOpacity   = EffectiveOpacity

     call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: HS % PositionSpace )
    select type ( PS => HS % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    select case ( PS % nDimensions )
      case ( 1 )
        CoordinateSystem = 'SPHERICAL'
      case ( 2 ) 
         CoordinateSystem = 'CYLINDRICAL'
      case ( 3 )
         CoordinateSystem = 'CARTESIAN'
      case DEFAULT
        call show ( PS % nDimensions, 'nDimensions' )
        call Show ( 'Dimensionality not supported', CONSOLE % ERROR )
        call Show ( 'HomogeneousSphere_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
    end select
    
    call PROGRAM_HEADER % GetParameter ( CoordinateSystem, 'CoordinateSystem' )

    MaxRadius = 5.0_KDR
    call PROGRAM_HEADER % GetParameter ( MaxRadius, 'MaxRadius' )

    nCellsRadius = 100
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
             TransportOpacity, ED, EO, TO )
    end associate !-- ED, etc.
    end select !-- I
    end select !-- IC

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: HS % Current_ASC )
    select type ( RMA => HS % Current_ASC )  !-- RadiationMoments Atlas
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
    
     !-- Diagnostics

    allocate ( HS % Reference )
    allocate ( HS % Difference )
    call HS % Reference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Reference' )
    call HS % Difference % Initialize &
           ( PS, 'GENERIC', NameOutputOption = 'Difference' )
    HS % SetReference => SetReference

    call PrepareInterpolation ( SI )
    RM => HS % Reference % RadiationMoments ( )
    G => PS % Geometry ( )
    select type ( Grid => PS % Chart )
    class is ( Chart_SL_Template )
    associate &
      ( J   => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
        FF  => RM % Value ( :, RM % FLUX_FACTOR ), &
        VEF => RM % Value ( :, RM % VARIABLE_EDDINGTON_FACTOR ), &
        R   => G % Value ( :, G % CENTER ( 1 ) ) )

      do iV = 1, size ( J )
        if ( .not. Grid % IsProperCell ( iV ) ) cycle
        call SI ( iCOMOVING_ENERGY_DENSITY_SI ) &
                % Evaluate ( R ( iV ), J ( iV ) )
        call SI ( iFLUX_FACTOR_SI ) &
                % Evaluate ( R ( iV ), FF ( iV ) ) 
        call SI ( iVARIABLE_EDDINGTON_FACTOR_SI ) &
                % Evaluate ( R ( iV ), VEF ( iV ) ) 

      end do
      
      end associate !-- J, etc. 
      end select !
    !-- Initial conditions
    
    RM => RMA % RadiationMoments ( )
    call SetRadiation ( HS, RM )
    
    !-- Initialize template
    
    call HS % InitializeTemplate_C ( Name, FinishTimeOption = 15.0_KDR )
    
    !-- Cleanup
           
    end select !-- RMA
    end associate !-- IA
    end select !-- PS
    nullify ( G, RM )
    
  end subroutine Initialize 


  subroutine ComputeError ( HS )

    class ( HomogeneousSphereForm ), intent ( in ) :: &
      HS
    
    real ( KDR ) :: &
      L1_J, &
      L1_FF, &
      L1_VEF
    class ( RadiationMomentsForm ), pointer :: &
      RM         
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( PS => HS % PositionSpace )
    class is ( Atlas_SC_Form ) 
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
    
    RM => HS % Difference % RadiationMoments ( )
    
    associate &
      ( Difference_J   => RM % Value ( :, RM % COMOVING_ENERGY_DENSITY ), &
        Difference_FF  => RM % Value ( :, RM % FLUX_FACTOR ), &
        Difference_VEF => RM % Value ( :, RM % VARIABLE_EDDINGTON_FACTOR ) )
    call CO % Initialize ( PS % Communicator, [ 4 ], [ 4 ] )
    CO % Outgoing % Value ( 1 ) = sum ( abs ( Difference_J ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 2 ) = sum ( abs ( Difference_FF ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 3 ) = sum ( abs ( Difference_VEF ), &
                                        mask = C % IsProperCell )
    CO % Outgoing % Value ( 4 ) = C % nProperCells
    call CO % Reduce ( REDUCTION % SUM )
    end associate !-- Difference_J, etc.
    
    associate &
      ( DifferenceSum_J   => CO % Incoming % Value ( 1 ), &
        DifferenceSum_FF  => CO % Incoming % Value ( 2 ), &
        DifferenceSum_VEF => CO % Incoming % Value ( 3 ), &
        nValues => CO % Incoming % Value ( 4 ) )
    L1_J   = DifferenceSum_J   / nValues
    L1_FF  = DifferenceSum_FF  / nValues
    L1_VEF = DifferenceSum_VEF / nValues
    end associate

    call Show ( L1_J, '*** L1_J error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_FF, '*** L1_FF error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )
    call Show ( L1_VEF, '*** L1_VEF error', nLeadingLinesOption = 2, &
                nTrailingLinesOption = 2 )

    end select !-- C
    end select !-- PS

  end subroutine ComputeError

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
    
    type ( SplineInterpolationForm ), dimension ( 1 ) :: &
      SI
    class ( GeometryFlatForm ), pointer :: &
      G
    integer :: &
      iV
    
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
      nullify ( G )

  end subroutine SetRadiation


  subroutine SetInteractions &
               ( HS, CoordinateSystem, EquilibriumDensity, EffectiveOpacity, &
                 TransportOpacity, ED, EO, TO )
    class ( HomogeneousSphereForm ), intent ( in ) :: &
      HS
    character ( LDL ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), intent ( in ) :: &
      EquilibriumDensity,&
      EffectiveOpacity, &
      TransportOpacity
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


  subroutine SetReference ( HS )

    class ( IntegratorTemplate ), intent ( in ) :: &
      HS

    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( HS )
    class is ( HomogeneousSphereForm )

    select type ( RMA => HS % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- RMA

    RM_R => HS % Reference % RadiationMoments ( )

    RM_D => HS % Difference % RadiationMoments ( )
!    RM_D % Value  =  RM % Value  -  RM_R % Value
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- HS
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


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


  subroutine PrepareInterpolation ( SI )

    type ( SplineInterpolationForm ), dimension ( 3 ), intent ( inout ) :: &
      SI

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Slope_J, &
      Slope_FF, &
      Slope_VEF
    real ( KDR ), dimension ( : ), allocatable :: &
      RC, &   !-- RadiusCenter
      dRC, &  !-- WidthCenter
      Radius, &
      ComovingEnergyDensity, &
      FluxFactor, &
      VariableEddingtonFactor
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Profile
    character ( LDF ) :: &
      Path, &
      Filename
    type ( TableStreamForm ) :: &
      TS

    Path = '../Parameters/'
    Filename = 'Oconnor_Abdikamalov.curve'
    
    call PROGRAM_HEADER % GetParameter &
           ( Filename, 'Filename' )

    call TS % Initialize &
           ( Filename, PROGRAM_HEADER % Communicator % Rank, &
             PathOption = Path )
    call TS % Read ( Profile, oRowOption = 1 )

    !-- Set "edge" values

    associate &
      (   R   => Profile ( :, iRADIUS_TS ), &                      !-- cell outter edge
          J   => Profile ( :, iCOMOVING_ENERGY_DENSITY_TS ), &     !-- cell center 
          FF  => Profile ( :, iFLUX_FACTOR_TS ), &                 !-- cell center
          VEF => Profile ( :, iVARIABLE_EDDINGTON_FACTOR_TS ), &   !-- cell center
          nProfile => size ( Profile, dim = 1 ) )

    allocate ( Radius ( nProfile + 1 ) )
    allocate ( dRC ( nProfile ) )
    allocate ( RC ( nProfile ) )
    Radius ( 1 )          =  0.0_KDR
    do iV = 2, nProfile + 1
      Radius  ( iV )  =  R ( iV - 1 )
      dRC ( iV - 1 )  =  Radius ( iV )  -  Radius ( iV - 1 )
      RC  ( iV - 1 )  =  Radius ( iV - 1 )  +  0.5_KDR * dRC ( iV - 1 )
    end do

    allocate ( ComovingEnergyDensity ( nProfile + 1 ) )
    allocate ( FluxFactor ( nProfile + 1 ) )
    allocate ( VariableEddingtonFactor ( nProfile + 1 ) )

    !-- First edge extrapolated
    Slope_J      =  ( J ( 2 )  -  J ( 1 ) )  &
                    /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_FF     =  ( FF ( 2 )  -  FF ( 1 ) )  &
                    /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
    Slope_VEF    =  ( VEF ( 2 )  -  VEF ( 1 ) )  &
                    /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )

    ComovingEnergyDensity ( 1 )    =  J ( 1 )  &
                                  +  Slope_J   * ( Radius ( 1 )  -  RC ( 1 ) )
    FluxFactor ( 1 )               =  FF ( 1 )  &
                                  +  Slope_FF   * ( Radius ( 1 )  -  RC ( 1 ) )
    VariableEddingtonFactor ( 1 )  =  VEF ( 1 )  &
                                  +  Slope_VEF   * ( Radius ( 1 )  -  RC ( 1 ) )

    do iV = 2, nProfile + 1

      if ( iV <= nProfile ) then
        Slope_J      =  ( J ( iV )  - J ( iV - 1 ) )  &
                        /  ( 0.5_KDR * ( dRC ( iV - 1 )  +  dRC ( iV ) ) )
        Slope_FF     =  ( FF ( iV )  -  FF ( iV - 1 ) )  &
                        /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
        Slope_VEF    =  ( VEF ( iV )  -  VEF ( iV - 1 ) )  &
                        /  ( 0.5_KDR * ( dRC ( 1 )  +  dRC ( 2 ) ) )
      else
        !-- Last edge extrapolated with same slope
      end if

      ComovingEnergyDensity ( iV )  &
        =  J   ( iV - 1 )  +  Slope_J * ( Radius ( iV )  -  RC ( iV - 1 ) )
      FluxFactor ( iV ) &
        =  FF ( iV - 1 ) +  Slope_FF   * ( Radius ( iV )  -  RC ( iV - 1 ) )
      VariableEddingtonFactor ( iV )  &
           =  VEF ( iV -1 ) +  Slope_VEF   * ( Radius ( iV )  -  RC ( iV - 1 ) )
    end do

    end associate !-- R, etc.

    call Show ( 'First few values' )
    call Show ( Profile ( 1 : 5, iRADIUS_TS ), 'RadiusTable' )
    call Show ( Radius ( 1 : 5 ), 'RadiusEdge' )
    call Show ( Profile ( 1 : 5, iCOMOVING_ENERGY_DENSITY_TS ), &
                'ComovingEnergyDensityTable' )
    call Show ( ComovingEnergyDensity ( 1 : 5 ), 'ComovingEnergyDensityEdge' )

    !-- SplineInterpolation initialization

    call SI ( iCOMOVING_ENERGY_DENSITY_SI ) % Initialize &
           ( Radius, ComovingEnergyDensity, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iFLUX_FACTOR_SI ) % Initialize &
           ( Radius, FluxFactor, VerbosityOption = CONSOLE % INFO_3 )
    call SI ( iVARIABLE_EDDINGTON_FACTOR_SI ) % Initialize &
           ( Radius, VariableEddingtonFactor, VerbosityOption = CONSOLE % INFO_3 )

  end subroutine PrepareInterpolation

end module HomogeneousSphere_Form
