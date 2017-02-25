module RiemannProblem_Form

  use Basics
  use Mathematics
  use Fluid_P_P__Form
  use Fluid_ASC__Form
  
  implicit none
  private
  
  type, public, extends ( Integrator_C_PS_Template ) :: RiemannProblemForm
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize  
  end type RiemannProblemForm

    private :: &
      SetFluid

    real ( KDR ), private :: &
      Density_L, Density_R, &  !-- Left and Right
      Pressure_L, Pressure_R, &
      Energy_L, Energy_R, &
      Speed_L, Speed_R, &
      AdiabaticIndex, &
      SinTheta, CosTheta, &
      SinPhi, CosPhi
    real ( KDR ), dimension ( 3 ) :: &
      DP_1, DP_2, DP_3  !-- DiscontinuityPoint_1, etc.
    type ( MeasuredValueForm ), private :: &
      DensityUnit, &
      EnergyUnit, &
      SpeedUnit

contains


  subroutine Initialize ( RP, Name )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    real ( KDR ), dimension ( 3 ) :: &
      Normal, &
      UnitNormal

    if ( RP % Type == '' ) &
      RP % Type = 'a RiemannProblem' 

    !-- PositionSpace

    allocate ( Atlas_SC_Form :: RP % PositionSpace )
    select type ( PS => RP % PositionSpace )
    class is ( Atlas_SC_Form )
    call PS % Initialize ( Name, PROGRAM_HEADER % Communicator )

    do iD = 1, PS % nDimensions
      call PS % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iDimension = iD )
    end do !-- iD

    call PS % CreateChart ( )

    !-- Geometry of PositionSpace

    allocate ( RP % Geometry_ASC )
    associate ( GA => RP % Geometry_ASC )  !-- GeometryAtlas
    call GA % Initialize ( PS )
    call PS % SetGeometry ( GA )
    end associate !-- GA

    !-- Fluid

    allocate ( Fluid_ASC_Form :: RP % Current_ASC )
    select type ( FA => RP % Current_ASC )
    class is ( Fluid_ASC_Form )
    call FA % Initialize ( PS, 'POLYTROPIC' )
    end select !-- FA

    !-- Step

    allocate ( Step_RK2_C_ASC_Form :: RP % Step )
    select type ( S => RP % Step )
    class is ( Step_RK2_C_ASC_Form )
    call S % Initialize ( Name )
    end select !-- S

    !-- Problem definition

    !-- Left and Right states

    Density_L      = 1.0_KDR
    Pressure_L     = 1.0_KDR
    Speed_L        = 0.0_KDR
    Density_R      = 0.125_KDR
    Pressure_R     = 0.1_KDR
    Speed_R        = 0.0_KDR
    AdiabaticIndex = 1.4_KDR

    DensityUnit = UNIT % IDENTITY
    EnergyUnit  = UNIT % IDENTITY
    SpeedUnit   = PS % Chart % CoordinateUnit ( 1 ) / RP % TimeUnit
 
    call PROGRAM_HEADER % GetParameter &
           ( Density_L, 'DensityLeft', InputUnitOption = DensityUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Pressure_L, 'PressureLeft', InputUnitOption = EnergyUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Speed_L, 'SpeedLeft', InputUnitOption = SpeedUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Density_R, 'DensityRight', InputUnitOption = DensityUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Pressure_R, 'PressureRight', InputUnitOption = EnergyUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Speed_R, 'SpeedRight', InputUnitOption = SpeedUnit )
    call PROGRAM_HEADER % GetParameter ( AdiabaticIndex, 'AdiabaticIndex' )

    Energy_L = Pressure_L / ( AdiabaticIndex - 1.0_KDR )
    Energy_R = Pressure_R / ( AdiabaticIndex - 1.0_KDR )

    !-- Three points define the plane of discontinuity

    DP_1 = [ 0.5_KDR, 0.0_KDR, 0.0_KDR ]
    DP_2 = [ 0.0_KDR, 0.5_KDR, 0.0_KDR ]
    DP_3 = [ 0.0_KDR, 0.0_KDR, 0.5_KDR ]

    if ( PS % nDimensions < 3 ) DP_3 ( 3 ) = 0.1 * sqrt ( huge ( 1.0_KDR ) )
    if ( PS % nDimensions < 2 ) DP_2 ( 2 ) = 0.1 * sqrt ( huge ( 1.0_KDR ) )

    call PROGRAM_HEADER % GetParameter &
           ( DP_1 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_1' )
    if ( PS % nDimensions > 1 ) &
      call PROGRAM_HEADER % GetParameter &
             ( DP_2 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_2' )
    if ( PS % nDimensions > 2 ) &
      call PROGRAM_HEADER % GetParameter &
             ( DP_3 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_3' )

    !-- Normal vector ( DP_2 - DP_1 ) x ( DP_3 - DP_1 )

    Normal ( 1 ) &
      = DP_3 ( 2 ) * ( DP_1 ( 3 ) - DP_2 ( 3 ) ) &
          + DP_1 ( 2 ) * (   DP_2 ( 3 ) - DP_3 ( 3 ) ) &
          + DP_2 ( 2 ) * ( - DP_1 ( 3 ) + DP_3 ( 3 ) )
    Normal ( 2 ) &
      = DP_3 ( 1 ) * ( - DP_1 ( 3 ) + DP_2 ( 3 ) ) &
          + DP_2 ( 1 ) * (   DP_1 ( 3 ) - DP_3 ( 3 ) ) &
          + DP_1 ( 1 ) * ( - DP_2 ( 3 ) + DP_3 ( 3 ) )
    Normal ( 3 ) &
      = DP_3 ( 1 ) * ( DP_1 ( 2 ) - DP_2 ( 2 )) &
          + DP_1 ( 1 ) * (   DP_2 ( 2 ) - DP_3 ( 2 ) ) &
          + DP_2 ( 1 ) * ( - DP_1 ( 2 ) + DP_3 ( 2 ) )
    Normal = Normal / maxval ( Normal ) !-- to avoid overflow in the next line
    UnitNormal = Normal / sqrt ( dot_product ( Normal, Normal ) )
    call Show ( UnitNormal, 'UnitNormal', CONSOLE % INFO_3 )
  
    CosTheta = dot_product ( UnitNormal, [ 0.0_KDR, 0.0_KDR, 1.0_KDR ] )
    SinTheta = sqrt ( 1.0_KDR - CosTheta ** 2 )
    if ( SinTheta /= 0.0_KDR ) then
      CosPhi &
        = dot_product ( UnitNormal, [ 1.0_KDR, 0.0_KDR, 0.0_KDR ] ) &
            / SinTheta
      SinPhi &
        = dot_product ( UnitNormal, [ 0.0_KDR, 1.0_KDR, 0.0_KDR ] ) &
            / SinTheta
    else
      CosPhi = 1.0_KDR
      SinPhi = 0.0_KDR
    end if
    call Show ( 'Angles defining normal vector', CONSOLE % INFO_3 )
    call Show ( SinTheta, 'SinTheta', CONSOLE % INFO_3 )
    call Show ( CosTheta, 'CosTheta', CONSOLE % INFO_3 )
    call Show ( SinPhi, 'SinPhi', CONSOLE % INFO_3 )
    call Show ( CosPhi, 'CosPhi', CONSOLE % INFO_3 )

    !-- Set fluid and initialize Integrator template

    call SetFluid ( RP )
    call RP % InitializeTemplate_C_PS ( Name )
    
    !-- Cleanup

    end select !-- PS

  end subroutine Initialize

  
  subroutine Finalize ( RP )

    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    call RP % FinalizeTemplate_C_PS ( )

  end subroutine Finalize


  subroutine SetFluid ( RP )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_P_Form ), pointer :: &
      F

    select type ( PS => RP % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( FA => RP % Current_ASC )
    class is ( Fluid_ASC_Form )

    F => FA % Fluid_P_P ( )

    call F % SetAdiabaticIndex ( AdiabaticIndex )

    !-- Translate to origin, rotate normal to xz plane and then to x axis

    associate &
      ( X  => G % Value ( :, G % CENTER ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER ( 3 ) ), &
        N  => F % Value ( :, F % COMOVING_DENSITY ), &
        VX => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        VY => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        VZ => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
        E  => F % Value ( :, F % INTERNAL_ENERGY ) )

    where ( SinTheta * CosPhi * ( X  - DP_1 ( 1 ) ) &
              + SinTheta * SinPhi * ( Y - DP_1 ( 2 ) ) &
              + CosTheta * ( Z - DP_1 ( 3 ) ) <= 1.e-10_KDR )     
      N  = Density_L
      E  = Energy_L 
      VX = Speed_L * SinTheta * CosPhi
      VY = Speed_L * SinTheta * SinPhi
      VZ = Speed_L * CosTheta
    elsewhere
      N  = Density_R  
      E  = Energy_R 
      VX = Speed_R * SinTheta * CosPhi
      VY = Speed_R * SinTheta * SinPhi
      VZ = Speed_R * CosTheta
    end where

    call F % ComputeFromPrimitive ( G )

    end associate !-- X, etc.
    end select !-- FA
    end select !-- PS
    nullify ( F, G )
    
  end subroutine SetFluid

  
end module RiemannProblem_Form
