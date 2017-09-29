module DynamicDiffusion_Form

  use Basics
  use Mathematics
  use Interactions_Template
  use Interactions_F__Form
  use Interactions_DD_G__Form
  use Interactions_ASC__Form
  use Interactions_CSL__Form
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use ApplyRelaxation_RM__Command

!-- Try long compared to diffusion time scale 
!-- i.e. large opacity corresponds to advection times much shorter than diffusion times

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: DynamicDiffusionForm
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
  end type DynamicDiffusionForm

    private :: &
      SetRadiation, &
      SetReference, &
      SetInteractions, &
      ApplySources_Radiation

      private :: &
        SetRadiationKernel, &
        SetReferenceKernel

      real ( KDR ), private :: &
        EquilibriumDensity, &
        EffectiveOpacity, &
        TransportOpacity, &
        X0, &
        Y0, &
        D0, &
        J0, &
        T0, &
        Delta, &
        Beta
      character ( LDF ) :: &
        InteractionsType
      class ( InteractionsTemplate ), pointer :: &
        Interactions => null ( )
      

contains


  subroutine Initialize ( DD, Name )

    class ( DynamicDiffusionForm ), intent ( inout ) :: &
      DD
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

    X0     = 1.0_KDR 
    Y0     = 1.0_KDR 
    D0     = 3.0e-3_KDR
    J0     = 1.0_KDR
    T0     = 5.0_KDR
    Delta  = 0.4_KDR

    Beta   = 0.0_KDR
    call PROGRAM_HEADER % GetParameter ( Beta, 'Beta' )

    EquilibriumDensity = 0.0_KDR
    EffectiveOpacity   = 0.0_KDR
    TransportOpacity   = D0 ** ( -1.0 ) / 3.0_KDR

    call PROGRAM_HEADER % GetParameter &
           ( EquilibriumDensity, 'EquilibriumDensity' )
    call PROGRAM_HEADER % GetParameter &
           ( EffectiveOpacity, 'EffectiveOpacity' )
    call PROGRAM_HEADER % GetParameter &
           ( TransportOpacity, 'TransportOpacity' )

    InteractionsType = 'DYNAMIC_DIFFUSION_GREY'
    call PROGRAM_HEADER % GetParameter ( InteractionsType, 'InteractionsType' )

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: DD % PositionSpace )
    select type ( PS => DD % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( 'PositionSpace', PROGRAM_HEADER % Communicator )
    
    select case ( PS % nDimensions )
      case ( 1 ) 
         CoordinateSystem = 'CARTESIAN'
      case ( 2 ) 
         CoordinateSystem = 'CARTESIAN'
      case DEFAULT
        call show ( PS % nDimensions, 'nDimensions' )
        call Show ( 'Dimensionality not supported', CONSOLE % ERROR )
        call Show ( 'DynamicDiffusion_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
    end select

    call PROGRAM_HEADER % GetParameter ( CoordinateSystem, 'CoordinateSystem' )

    select case ( CoordinateSystem )
      case ( 'CARTESIAN' )
        MinCoordinate = [ 0.0_KDR, 1.0_KDR, 0.0_KDR ]
        MaxCoordinate = [ 3.0_KDR, 1.0_KDR, 0.0_KDR ]
        call PS % SetBoundaryConditionsFace &
               ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 1 )
        nCells = [ 512, 1, 1 ]
        if ( PS % nDimensions > 1 ) then
          MinCoordinate = [ 0.0_KDR, 0.0_KDR, 0.0_KDR ]
          MaxCoordinate = [ 3.0_KDR, 2.0_KDR, 0.0_KDR ]
          call PS % SetBoundaryConditionsFace &
                 ( [ 'OUTFLOW', 'OUTFLOW' ], iDimension = 2 )
           nCells = [ 512, 256, 1 ]
        end if 
 
        call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )
      case DEFAULT 
        call show ( CoordinateSystem, 'CoordinateSystem' )
        call Show ( 'Coordinate System not implemented', CONSOLE % ERROR )
        call Show ( 'DynamicDiffusion_Form', 'module', CONSOLE % ERROR )
        call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
    end select !-- CoordinateSystem

    call PS % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells )

    !-- Geometry of PositionSpace

    allocate ( DD % Geometry_ASC )
    associate ( GA => DD % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Interactions
    allocate ( DD % Interactions_ASC )
    associate ( IA => DD % Interactions_ASC )

    call IA % Initialize ( PS, InteractionsType )

    !-- RadiationMoments ( Generic )

    allocate ( RadiationMoments_ASC_Form :: DD % Current_ASC )
    select type ( RMA => DD % Current_ASC )  !-- FluidAtlas
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize ( PS, 'GENERIC' )
    call RMA % SetInteractions ( IA )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: DD % Step )
    select type ( S => DD % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( RMA, Name )
    S % ApplyRelaxation % Pointer => ApplyRelaxation_RM
    S % ApplySources    % Pointer => ApplySources_Radiation
    end select !-- S

    !-- Diagnostics

    allocate ( DD % Reference )
    allocate ( DD % Difference )
    call DD % Reference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    call DD % Difference % Initialize &
           ( PS, 'GENERIC', NameShortOption = 'Difference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )

    DD % SetReference => SetReference
    
    !-- Initial conditions
    
    RM => RMA % RadiationMoments ( )
    call SetRadiation ( DD, RM )
    call SetInteractions ( DD )
    
    !-- Initialize template
    
    call DD % InitializeTemplate_C_PS ( Name, FinishTimeOption = 5.0_KDR )
    
    !-- Cleanup
           
    end select !-- RMA
    end associate !-- IA
    end select !-- PS
    nullify ( RM )
    
  end subroutine Initialize 


  subroutine ComputeError ( DD )

    class ( DynamicDiffusionForm ), intent ( in ) :: &
      DD
    
    real ( KDR ) :: &
      L1       
    class ( RadiationMomentsForm ), pointer :: &
      RM         
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( PS => DD % PositionSpace )
    class is ( Atlas_SC_Form ) 
    select type ( C => PS % Chart )
    class is ( Chart_SL_Template )
    
    RM => DD % Difference % RadiationMoments ( )
    
    associate &
      ( Difference => RM % Value ( :, RM % COMOVING_ENERGY ) )
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


  impure elemental subroutine Finalize ( DD )

    type ( DynamicDiffusionForm ), intent ( inout ) :: &
      DD

    if ( allocated ( DD % Difference ) ) &
      deallocate ( DD % Difference )
    if ( allocated ( DD % Reference ) ) &
      deallocate ( DD % Reference )

    call DD % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetRadiation ( DD, RM )

    class ( DynamicDiffusionForm ), intent ( in ) :: &
      DD
    class ( RadiationMomentsForm ), intent ( inout ) :: &
      RM

    class ( GeometryFlatForm ), pointer :: &
      G
    character ( LDL ) :: &
      CoordinateSystem
    
    select type ( PS => DD % PositionSpace )
      class is ( Atlas_SC_Form )
      G => PS % Geometry ( )

      associate &
        ( X =>  G % Value ( :, G % CENTER ( 1 ) ), &
          Y =>  G % Value ( :, G % CENTER ( 2 ) ) )
        
      call SetRadiationKernel &
             ( DD, X, Y, &
               J  = RM % Value ( :, RM % COMOVING_ENERGY ), &
               HX = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               HY = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               HZ = RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               VX = RM % Value ( :, RM % FLUID_VELOCITY_U ( 1 ) ), &
               VY = RM % Value ( :, RM % FLUID_VELOCITY_U ( 2 ) ), &
               VZ = RM % Value ( :, RM % FLUID_VELOCITY_U ( 3 ) ) )

      call RM % ComputeFromPrimitive ( G )
   
      end associate !-- X, etc.
    end select    !-- PS
    
    nullify ( G )

  end subroutine SetRadiation


  subroutine SetReference ( DD )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      DD

    real ( KDR ), allocatable, dimension ( : ) :: &
      R_sq
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( DD )
    class is ( DynamicDiffusionForm )

    select type ( RMA => DD % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RMA % RadiationMoments ( )
    end select !-- RMA

    select type ( PS => DD % PositionSpace )
      class is ( Atlas_SC_Form )
      G => PS % Geometry ( )
      
      allocate ( R_sq ( size ( G % Value ( :, G % CENTER ( 1 ) ) ) ) )

      associate &
        ( X  =>  G % Value ( :, G % CENTER ( 1 ) ), &
          Y  =>  G % Value ( :, G % CENTER ( 2 ) ), &
          nD => PS % nDimensions )

      R_sq =  ( X - X0 ) ** 2

      if ( nD > 1 ) R_sq = R_sq + ( Y - Y0 ) ** 2
      
      end associate !-- X, etc.
    end select !-- PS

    RM_R => DD % Reference % RadiationMoments ( )

    associate ( J => RM_R % Value ( :, RM_R % COMOVING_ENERGY ) )

    call SetReferenceKernel ( R_sq, DD % Time, J )

    end associate !-- E

    RM_D => DD % Difference % RadiationMoments ( )
!    RM_D % Value  =  RM % Value  -  RM_R % Value
    call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- DD
    nullify ( G, RM, RM_R, RM_D )
    deallocate ( R_sq )

  end subroutine SetReference


  subroutine SetInteractions ( DD )
    type ( DynamicDiffusionForm ), intent ( inout ) :: &
      DD           
    
    class ( RadiationMomentsForm ), pointer :: &
      RM
    class ( InteractionsTemplate), pointer :: &
      I
    class ( GeometryFlatForm ), pointer :: &
      G

    associate ( IA => DD % Interactions_ASC )

    I => IA % Interactions ( )

    select type ( I )
      type is ( Interactions_F_Form )
         call Show ( 'Interactions_F_Form' )
        associate &
          ( Xi_J  => I % Value ( :, I % EMISSIVITY_J ), &
            Chi_J => I % Value ( :, I % OPACITY_J ), &
            Chi_H => I % Value ( :, I % OPACITY_H ) )

        call SetFixedInteractions &
               ( DD, EquilibriumDensity, &
                 EffectiveOpacity, TransportOpacity, &
                 Xi_J, Chi_J, Chi_H )
        end associate !-- Xi_i, etc.
      type is ( Interactions_DD_G_Form )
        call Show ( '>>>Interactions_DD_G_Form' )
        select type ( RMA => DD % Current_ASC )
          class is ( RadiationMoments_ASC_Form )
          RM => RMA % RadiationMoments ( )
        select type ( PS => DD % PositionSpace )
               class is ( Atlas_SC_Form )
                 G => PS % Geometry ( )
        call I % Set ( RM, G, DD % Time, D0, Delta, X0, Y0 )
        end select !-- G

        call I % Compute ( RM )

        !-- Module variable for accessibility in ApplySources_Fluid below
        Interactions => I
        
        end select !-- RMA
      class default
      call Show ( 'Interactions type not recognized', CONSOLE % ERROR )
      call Show ( 'DynamicDiffusion_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetInteractions', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- I

    end associate !-- IA

  end subroutine SetInteractions


  subroutine ApplySources_Radiation &
               ( S, Sources_RM, Increment, Radiation, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Radiation
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    !-- No sources applied here; just an occasion to compute interactions
    !   to be used in relaxation.

    if ( associated ( Interactions ) ) then
      call Interactions % Compute ( Radiation )
    end if

  end subroutine ApplySources_Radiation


  subroutine SetFixedInteractions &
               ( DD, EquilibriumDensity, &
                 EffectiveOpacity, TransportOpacity, &
                 Xi_J, Chi_J, Chi_H )

    class ( DynamicDiffusionForm ), intent ( in ) :: &
      DD
    real ( KDR ), intent ( in ) :: &
      EquilibriumDensity,&
      EffectiveOpacity, &
      TransportOpacity
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ) :: &
      R_sq
    integer ( KDI ) :: &
     iV, & !-- iValue
     nValues

    select type ( PS => DD % PositionSpace )
      class is ( Atlas_SC_Form )
        G => PS % Geometry ( )
      
       nValues = size ( Xi_J )

      associate &
        ( X  =>  G % Value ( :, G % CENTER ( 1 ) ), &
          Y  =>  G % Value ( :, G % CENTER ( 2 ) ) )

      ! !$OMP parallel do private ( iV )
      ! do iV = 1, nValues
      !   Xi_J  ( iV ) = EffectiveOpacity * EquilibriumDensity
      !   Chi_J ( iV ) = EffectiveOpacity 
      !   Chi_H ( iV ) = TransportOpacity 
      ! end do !-- iV
      ! !$OMP end parallel do

      !$OMP parallel do private ( iV )
      do iV = 1, nValues
        R_sq         =  ( X ( iV ) - X0 ) ** 2 + ( Y ( iV ) - Y0 ) ** 2
        Xi_J  ( iV ) = EffectiveOpacity * EquilibriumDensity
        Chi_J ( iV ) = EffectiveOpacity 
        Chi_H ( iV ) = 1.0_KDR / ( 3.0_KDR * D0 ) * exp ( - R_sq / Delta ** 2 )
      end do !-- iV
      !$OMP end parallel do

      end associate !-- X, etc.
   end select !-- PS

  end subroutine SetFixedInteractions


  subroutine SetRadiationKernel ( DD, X, Y, J, HX, HY, HZ, VX, VY, VZ )

    class ( DynamicDiffusionForm ), intent ( in ) :: &
      DD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, &
      Y
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J, &
      HX, HY, HZ, &
      VX, VY, VZ

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      R_sq

    nValues = size ( X )

    select type ( PS => DD % PositionSpace )
      class is ( Atlas_SC_Form )
    associate ( nD => PS % nDimensions ) 

      !$OMP parallel do private ( iV )
      do iV = 1, nValues
        R_sq       =  ( X ( iV ) - X0 ) ** 2 + ( Y ( iV ) - Y0 ) ** 2
        J ( iV )   = max ( J0 * exp ( - R_sq / ( 4.0_KDR * D0 * T0 ) ), &
                           sqrt ( tiny ( 0.0_KDR ) ) )
        HX ( iV )  = J0 * ( X ( iV ) - X0 ) / ( 2.0_KDR * T0 ) &
                           * exp ( - R_sq / ( 4.0_KDR * D0  * T0 ) )
        HY ( iV )  = J0 * ( Y ( iV ) - Y0 ) / ( 2.0_KDR * T0 ) &
                           * exp ( - R_sq / ( 4.0_KDR * D0 * T0 ) )
        HZ ( iV )  = 0.0_KDR
        VX ( iV )  = Beta * CONSTANT % SPEED_OF_LIGHT
        VY ( iV )  = 0.0_KDR
        VZ ( iV )  = 0.0_KDR
      end do !-- iV
      !$OMP end parallel do
    
    end associate !-- nD
    end select !-- PS

  end subroutine SetRadiationKernel


  subroutine SetReferenceKernel ( R_sq, T, J )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R_sq
    real ( KDR ), intent ( in ) :: &
      T
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      J

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      J0, &
      T0

    nValues = size ( R_sq )

    J0 = 1.0_KDR
    T0 = 5.0_KDR

   
    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      J ( iV )   = J0 * ( T0 / ( T0 + T ) ) &
                   * exp ( - R_sq ( iV ) / ( 4.0_KDR * D0 * ( T0 + T ) ) )
    end do !-- iV
    !$OMP end parallel do

  end subroutine SetReferenceKernel


end module DynamicDiffusion_Form
