module StreamingProfile_Form

  !-- Muller et al. 2010 Test 3.1

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use RadiationMoments_ASC__Form
  use ApplyCurvilinear_RM__Command

  implicit none
  private

  type, public, extends ( Integrator_C_PS_Template ) :: StreamingProfileForm
     type ( RadiationMoments_ASC_Form ), allocatable :: &
      Reference
     real ( KDR ) :: &
       InnerRadius, &
       OutterRadius, &
       InflectionRadius
  contains
    procedure, public, pass :: &
      Initialize 
    final :: &
      Finalize
  end type StreamingProfileForm

    private :: &
      SetRadiation, &
      SetInteractions, &
      SetReference

    real ( KDR ), dimension ( 3 ), private :: &
      MinCoordinate, &
      MaxCoordinate
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

contains

  subroutine Initialize ( SP, Name )
  
    class ( StreamingProfileForm ), intent ( inout ) :: &
      SP
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    integer ( KDI ), dimension ( 3 ) :: &
      nCells, &
      nGhostLayers
    real ( KDR ) :: &
      FinishTime, &
      dR
    real ( KDR ), dimension ( 3 ) :: &
      Ratio
    type ( MeasuredValueForm ) :: &
      TimeUnit, &
      MassDensityUnit, &
      EnergyDensityUnit, &
      MassUnit, &
      EnergyUnit, &
      MomentumUnit
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit, &
      VelocityUnit, &
      MomentumDensity_U_Unit, &
      MomentumDensity_D_Unit

    associate &
      ( IR     => SP % InnerRadius, &
        IFR    => SP % InflectionRadius, &
        OR     => SP % OutterRadius )

    IR  = 4.0_KDR
    IFR = 150.0_KDR
    OR  = 500.0_KDR 
    
    dR  = OR - IR

    call PROGRAM_HEADER % GetParameter ( IR,     'InnerRadius' )
    call PROGRAM_HEADER % GetParameter ( IFR,    'InflectionRadius' )
    call PROGRAM_HEADER % GetParameter ( OR,     'OutterRadius' )

    IR  = IR   * UNIT % KILOMETER
    IFR = IFR  * UNIT % KILOMETER
    OR  = OR  * UNIT % KILOMETER

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: SP % PositionSpace )
    select type ( PS => SP % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize &
         ( 'PositionSpace', PROGRAM_HEADER % Communicator )!, &
           !iDimensionalityOption = 1_KDI )

    CoordinateSystem = 'SPHERICAL'
    call PROGRAM_HEADER % GetParameter ( CoordinateSystem, 'CoordinateSystem' )

   ! do iD = 1, PS % nDimensions
      call PS % SetBoundaryConditionsFace &
             ( [ 'INFLOW ', 'OUTFLOW' ], iDimension = 1_KDI )
   ! end do !-- iD

    CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

    MinCoordinate = [ IR, 0.0_KDR, 0.0_KDR ]
    MaxCoordinate = [ OR, 0.0_KDR, 0.0_KDR ]

    nCells = [ 4000, 1, 1 ]
    nGhostLayers = [ 2, 1, 1 ]
    call PROGRAM_HEADER % GetParameter ( nCells, 'nCells' )

    dR = dR / nCells ( 1 ) 
    Spacing = [ 'PROPORTIONAL', 'EQUAL       ', 'EQUAL       ' ]
    Ratio   = [ dR, 0.0_KDR, 0.0_KDR ]

    call PS % CreateChart &
           ( CoordinateSystemOption = CoordinateSystem, &
             CoordinateUnitOption = CoordinateUnit, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayers )!, &
!             SpacingOption = Spacing, &
!             RatioOption = Ratio )

    call PS % SetGeometry ( )
    
    !-- Prepare Units
    
    TimeUnit = UNIT % SECOND
    
    VelocityUnit           =  UNIT % SPEED_OF_LIGHT
   ! MassDensityUnit        =  UNIT % MASS_DENSITY_CGS
    EnergyDensityUnit      =  UNIT % KILOMETER ** ( - 2 )
    MomentumDensity_U_Unit =  EnergyDensityUnit 
    MomentumDensity_D_Unit =  EnergyDensityUnit 
    !MassUnit               =  MassDensityUnit  *  CoordinateUnit ( 1 )
    !EnergyUnit             =  MassUnit  *  UNIT % SPEED_OF_LIGHT ** 2
    !MomentumUnit           =  MassUnit  *  UNIT % SPEED_OF_LIGHT

    

    !-- RadiationMoments

    allocate ( RadiationMoments_ASC_Form :: SP % Current_ASC )
    select type ( RMA => SP % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    call RMA % Initialize &
           ( PS, 'GENERIC', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             UseLimiterOption = .true., LimiterParameterOption = 2.0_KDR )

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: SP % Step )
    select type ( S => SP % Step )
    class is ( Step_RK2_C_ASC_Form )

    call S % Initialize ( SP, RMA, Name )

    S % ApplySources % Pointer =>  ApplyCurvilinear_RM

    end select !-- S

    !-- Initial conditions

    call SetRadiation ( SP )

    !-- Diagnostics

    allocate ( SP % Reference )
    call SP % Reference % Initialize &
           ( PS, 'GENERIC', &
             MomentumDensity_U_UnitOption = MomentumDensity_U_Unit, &
             MomentumDensity_D_UnitOption = MomentumDensity_D_Unit, &
             EnergyDensityUnitOption = EnergyDensityUnit, &
             NameShortOption = 'Reference', &
             AllocateSourcesOption = .false., &
             IgnorabilityOption = CONSOLE % INFO_2 )
    
    call SetReference ( SP )

    !-- Initialize template

    FinishTime  =  1.0e-2_KDR  *  UNIT % SECOND

    call SP % InitializeTemplate_C_PS &
           ( Name, TimeUnitOption = TimeUnit, FinishTimeOption = FinishTime )

    !-- Cleanup

    end select !-- RMA
    end select !-- PS
    end associate !-- IR, etc.

  end subroutine Initialize

  
  impure elemental subroutine Finalize ( SP )
    
    type ( StreamingProfileForm ), intent ( inout ) :: &
      SP

    call SP % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetRadiation ( SP )

    type ( StreamingProfileForm ), intent ( inout ) :: &
      SP
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( RadiationMomentsForm ), pointer :: &
      RM

    select type ( RA => SP % Current_ASC )
    class is ( RadiationMoments_ASC_Form )
    RM => RA % RadiationMoments ( )

    select type ( PS => SP % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    associate &
      (   J => RM % Value ( :, RM % COMOVING_ENERGY ), &
        H_1 => RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
        V_1 => RM % Value ( :, RM % FLUID_VELOCITY_U ( 1 ) ), &
         FF => RM % Value ( :, RM % FLUX_FACTOR ), &
         SF => RM % Value ( :, RM % STRESS_FACTOR ), &
          R =>  G % Value ( :, G % CENTER_U ( 1 ) ) )

    J    =  0.0_KDR
    H_1  =  0.0_KDR
    V_1  =  0.0_KDR
    
    where ( R < MinCoordinate ( 1 ) )
       J   = 1.0_KDR / ( 4.0_KDR * CONSTANT % PI * ( R ) ** 2 )
       H_1 = 1.0_KDR / ( 4.0_KDR * CONSTANT % PI * ( R ) ** 2 )
       FF  = 1.0_KDR
       SF  = 1.0_KDR
   end where

    where ( R >= SP % InflectionRadius * 0.9_KDR &
           .and. R < SP % InflectionRadius ) 
      V_1 = - 0.2_KDR * CONSTANT % SPEED_OF_LIGHT * &
              ( R - SP % InflectionRadius * 0.9_KDR ) &
              / ( SP % InflectionRadius * 0.1_KDR )
    end where

    where ( R >= SP % InflectionRadius )
      V_1 = - 0.2_KDR * CONSTANT % SPEED_OF_LIGHT &
              * ( SP % InflectionRadius / R ) ** 2.0_KDR
    end where

    call RM % ComputeFromPrimitive ( G )

    end associate !-- J, etc.
    end select !-- PS
    end select !-- RA
    nullify ( RM, G )

  end subroutine SetRadiation


  subroutine SetInteractions ( SP, Xi_J, Chi_J, Chi_H )

    class ( StreamingProfileForm ), intent ( in ) :: &
      SP
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H

    Xi_J  = 0.0_KDR
    Chi_J = 0.0_KDR
    Chi_H = 0.0_KDR

  end subroutine SetInteractions

  
  subroutine SetReference ( SP )

    class ( IntegratorTemplate ), intent ( inout ) :: &
      SP

    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      R
    class ( RadiationMomentsForm ), pointer :: &
      RM, &
      RM_R, &  !-- RM_Reference
      RM_D     !-- RM_Difference

    select type ( SP )
    class is ( StreamingProfileForm )
       
    select type ( RMA => SP % Current_ASC )  !-- RadiationMoments Atlas
    class is ( RadiationMoments_ASC_Form )
   
    RM => RMA % RadiationMoments ( )
    RM_R => SP % Reference % RadiationMoments ( )

    select type ( PS => SP % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )
    
    associate &
      (   J => RM_R % Value ( :, RM % COMOVING_ENERGY ), &
        H_1 => RM_R % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
        V_1 => RM % Value ( :, RM % FLUID_VELOCITY_U ( 1 ) ), &
        R   => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Pi  => CONSTANT % PI, &
        c   => CONSTANT % SPEED_OF_LIGHT )

    H_1 =  c * 1.0_KDR / ( 1.0_KDR + 2.0_KDR * V_1 ) &
               / ( 4.0_KDR * Pi * R ** 2.0 )
    
    end associate !-- J, etc.

   ! call MultiplyAdd ( RM % Value, RM_R % Value, -1.0_KDR, RM_D % Value )

    end select !-- G
    end select !-- RMA
    end select !-- HS
    nullify ( RM, RM_R, RM_D )

  end subroutine SetReference


end Module StreamingProfile_Form
