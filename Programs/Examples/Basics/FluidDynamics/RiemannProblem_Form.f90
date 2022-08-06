module RiemannProblem_Form

  use Basics
  use PolytropicFluid_Form
  use ConservationLawEvolution_Template

  implicit none
  private

  type, public, extends ( ConservationLawEvolutionTemplate ) :: &
    RiemannProblemForm
  type ( QuantityForm ), private :: &
    DensityUnit, &
    EnergyUnit, &
    SpeedUnit
  contains
    procedure, private, pass :: &
      Initialize_RP
    generic, public :: &
      Initialize => Initialize_RP
    final :: &
      Finalize
  end type RiemannProblemForm
  
    private :: &
      SetInitial

contains


  subroutine Initialize_RP ( RP )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    type ( StorageForm ) :: &
      Primitive
    type ( QuantityForm ) :: &
      DensityUnit,  &
      EnergyUnit, &
      SpeedUnit

    RP % Type = 'a RiemannProblem' 

    call RP % Initialize &
           ( PROGRAM_HEADER % Communicator, &
             BoundaryConditionOption = 'REFLECTING'  )

    associate ( DM => RP % DistributedMesh )

    allocate ( PolytropicFluidForm :: RP % ConservedFields )
    select type ( PF => RP % ConservedFields )
    type is ( PolytropicFluidForm )

    call PF % Initialize ( DM, NameOption = 'PolytropicFluid' )
    
    if ( RP % UseDevice ) &
      call PF % AllocateDevice ( )
    
    call Primitive % Initialize ( PF, iaSelectedOption = PF % iaPrimitive )
    call DM % SetGhostExchange ( Primitive )
    
    RP % DensityUnit  = UNIT % IDENTITY
    RP % EnergyUnit   = UNIT % IDENTITY
    RP % SpeedUnit    = DM % CoordinateUnit ( 1 ) / RP % TimeUnit

    call SetInitial ( PF, RP )
    
    call PF % ComputeAuxiliaryFromPressure ( PF % Value )
    call PF % ComputeConserved ( PF % Value, UseDeviceOption = .false. )
    
    call PF % SetOutputPolytropic &
           ( VelocityUnitOption = spread ( RP % SpeedUnit, 1, 3 ), &
             DensityUnitOption = RP % DensityUnit, &
             EnergyUnitOption = RP % EnergyUnit )
    
    end select    !-- PF
    end associate !-- DM

  end subroutine Initialize_RP 
  
  
  subroutine Finalize ( RP )

    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    if ( allocated ( RP % ConservedFields ) ) &
      deallocate ( RP % ConservedFields )

  end subroutine Finalize
  
  
  subroutine SetInitial ( PF, RP )
  
    class ( PolytropicFluidForm ), intent ( inout ) :: &
      PF
    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Density_L, Density_R, &  !-- Left and Right
      Pressure_L, Pressure_R, &
      Speed_L, Speed_R, &
      AdiabaticIndex, &
      SinTheta, CosTheta, &
      SinPhi, CosPhi
    real ( KDR ), dimension ( 3 ) :: &
      DP_1, DP_2, DP_3, &  !-- DiscontinuityPoint_1, etc.
      Normal, &
      UnitNormal
      
    associate ( DM => RP % DistributedMesh )

    !-- Left and Right states

    Density_L      = 1.0_KDR
    Pressure_L     = 1.0_KDR
    Speed_L        = 0.0_KDR
    Density_R      = 0.125_KDR
    Pressure_R     = 0.1_KDR
    Speed_R        = 0.0_KDR
    AdiabaticIndex = 1.4_KDR

    call PROGRAM_HEADER % GetParameter &
           ( Density_L, 'DensityLeft', InputUnitOption = RP % DensityUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Pressure_L, 'PressureLeft', InputUnitOption = RP % EnergyUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Speed_L, 'SpeedLeft', InputUnitOption = RP % SpeedUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Density_R, 'DensityRight', InputUnitOption = RP % DensityUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Pressure_R, 'PressureRight', InputUnitOption = RP % EnergyUnit )
    call PROGRAM_HEADER % GetParameter &
           ( Speed_R, 'SpeedRight', InputUnitOption = RP % SpeedUnit )
    call PROGRAM_HEADER % GetParameter ( AdiabaticIndex, 'AdiabaticIndex' )
    
    !-- Three points define the plane of discontinuity

    DP_1 = [ 0.500001_KDR, 0.0_KDR, 0.0_KDR ]
    DP_2 = [ 0.0_KDR, 0.500001_KDR, 0.0_KDR ]
    DP_3 = [ 0.0_KDR, 0.0_KDR, 0.500001_KDR ]

    if ( DM % nDimensions < 3 ) DP_3 ( 3 ) = 0.1 * sqrt ( huge ( 1.0_KDR ) )
    if ( DM % nDimensions < 2 ) DP_2 ( 2 ) = 0.1 * sqrt ( huge ( 1.0_KDR ) )

    call PROGRAM_HEADER % GetParameter &
           ( DP_1 ( 1 : DM % nDimensions ), 'DiscontinuityPoint_1' )
    if ( DM % nDimensions > 1 ) &
      call PROGRAM_HEADER % GetParameter &
             ( DP_2 ( 1 : DM % nDimensions ), 'DiscontinuityPoint_2' )
    if ( DM % nDimensions > 2 ) &
      call PROGRAM_HEADER % GetParameter &
             ( DP_3 ( 1 : DM % nDimensions ), 'DiscontinuityPoint_3' )

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

    !-- Translate to origin, rotate normal to xz plane and then to x axis

    associate &
      ( X     => DM % Position % Value ( :, 1 ), &
        Y     => DM % Position % Value ( :, 2 ), &
        Z     => DM % Position % Value ( :, 3 ), &
        VX    => PF % Value ( :, PF % VELOCITY ( 1 ) ), &
        VY    => PF % Value ( :, PF % VELOCITY ( 2 ) ), &
        VZ    => PF % Value ( :, PF % VELOCITY ( 3 ) ), &
        N     => PF % Value ( :, PF % COMOVING_DENSITY ), &
        P     => PF % Value ( :, PF % PRESSURE ), &
        Gamma => PF % Value ( :, PF % ADIABATIC_INDEX ) )

    !-- FIXME: This where block is not reliable with NAG on my laptop with 
    !          many cores
    ! where ( SinTheta * CosPhi * ( X  - DP_1 ( 1 ) ) &
    !           + SinTheta * SinPhi * ( Y - DP_1 ( 2 ) ) &
    !           + CosTheta * ( Z - DP_1 ( 3 ) ) < 0.0_KDR )     
    !   N  = Density_L
    !   P  = Pressure_L 
    !   VX = Speed_L * SinTheta * CosPhi
    !   VY = Speed_L * SinTheta * SinPhi
    !   VZ = Speed_L * CosTheta
    ! elsewhere
    !   N  = Density_R  
    !   P  = Pressure_R 
    !   VX = Speed_R * SinTheta * CosPhi
    !   VY = Speed_R * SinTheta * SinPhi
    !   VZ = Speed_R * CosTheta
    ! end where
    do iV = 1, PF % nValues
      if ( SinTheta * CosPhi * ( X ( iV ) - DP_1 ( 1 ) ) &
                + SinTheta * SinPhi * ( Y ( iV ) - DP_1 ( 2 ) ) &
                + CosTheta * ( Z ( iV ) - DP_1 ( 3 ) ) < 0.0_KDR ) &
      then
        N  ( iV ) = Density_L
        P  ( iV ) = Pressure_L 
        VX ( iV ) = Speed_L * SinTheta * CosPhi
        VY ( iV ) = Speed_L * SinTheta * SinPhi
        VZ ( iV ) = Speed_L * CosTheta
      else
        N  ( iV ) = Density_R  
        P  ( iV ) = Pressure_R 
        VX ( iV ) = Speed_R * SinTheta * CosPhi
        VY ( iV ) = Speed_R * SinTheta * SinPhi
        VZ ( iV ) = Speed_R * CosTheta
      end if 
    end do
    Gamma = AdiabaticIndex
    
    end associate !-- X, etc.
    
    end associate !-- DM

  end subroutine SetInitial


end module RiemannProblem_Form
