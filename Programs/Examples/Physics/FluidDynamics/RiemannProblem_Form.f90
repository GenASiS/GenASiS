module RiemannProblem_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( UniverseTemplate ) :: RiemannProblemForm
    real ( KDR ) :: &
      Density_L, Density_R, &  !-- Left and Right
      Pressure_L, Pressure_R, &
      Energy_L, Energy_R, &
      Speed_L, Speed_R, &
      AdiabaticIndex, &
      SinTheta, CosTheta, &
      SinPhi, CosPhi
    real ( KDR ), dimension ( 3 ) :: &
      DP_1, DP_2, DP_3  !-- DiscontinuityPoint_1, etc.
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type RiemannProblemForm

    private :: &
      SetFluid

contains


  subroutine Initialize ( RP, Name )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    character ( * ), intent ( in ) :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension
    real ( KDR ), dimension ( 3 ) :: &
      Normal, &
      UnitNormal
    type ( Character_1D_Form ), dimension ( 3 ) :: &
      BoundaryConditionsFace


    if ( RP % Type == '' ) &
      RP % Type = 'a RiemannProblem'

    call RP % InitializeTemplate ( Name )


    !-- Integrator

    allocate ( FluidBoxForm :: RP % Integrator )
    select type ( FB => RP % Integrator )
    type is ( FluidBoxForm )

    associate ( BCF => BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )     
    end do
    call FB % Initialize &
           ( Name, FluidType = 'IDEAL', GeometryType = 'GALILEAN', &
             BoundaryConditionsFaceOption = BCF )
    ! FB % SetReference => SetReference
    end associate !-- BCF

    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )

    select type ( FA => FB % Current_ASC )
    class is ( Fluid_ASC_Form )


    !-- Parameters

    RP % Density_L      = 1.0_KDR
    RP % Pressure_L     = 1.0_KDR
    RP % Speed_L        = 0.0_KDR
    RP % Density_R      = 0.125_KDR
    RP % Pressure_R     = 0.1_KDR
    RP % Speed_R        = 0.0_KDR
    RP % AdiabaticIndex = 1.4_KDR

    call PROGRAM_HEADER % GetParameter ( RP % Density_L, 'DensityLeft' )
    call PROGRAM_HEADER % GetParameter ( RP % Pressure_L, 'PressureLeft' )
    call PROGRAM_HEADER % GetParameter ( RP % Speed_L, 'SpeedLeft' )
    call PROGRAM_HEADER % GetParameter ( RP % Density_R, 'DensityRight' )
    call PROGRAM_HEADER % GetParameter ( RP % Pressure_R, 'PressureRight' )
    call PROGRAM_HEADER % GetParameter ( RP % Speed_R, 'SpeedRight' )
    call PROGRAM_HEADER % GetParameter ( RP % AdiabaticIndex, 'AdiabaticIndex' )

    RP % Energy_L = RP % Pressure_L / ( RP % AdiabaticIndex - 1.0_KDR )
    RP % Energy_R = RP % Pressure_R / ( RP % AdiabaticIndex - 1.0_KDR )


    !-- Three points define the plane of discontinuity

    RP % DP_1 = [ 0.5_KDR, 0.0_KDR, 0.0_KDR ]
    RP % DP_2 = [ 0.0_KDR, 0.5_KDR, 0.0_KDR ]
    RP % DP_3 = [ 0.0_KDR, 0.0_KDR, 0.5_KDR ]

    if ( PS % nDimensions < 3 ) RP % DP_3 ( 3 ) &
      = 0.1 * sqrt ( huge ( 1.0_KDR ) )
    if ( PS % nDimensions < 2 ) RP % DP_2 ( 2 ) &
      = 0.1 * sqrt ( huge ( 1.0_KDR ) )

    call PROGRAM_HEADER % GetParameter &
           ( RP % DP_1 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_1' )
    if ( PS % nDimensions > 1 ) &
      call PROGRAM_HEADER % GetParameter &
             ( RP % DP_2 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_2' )
    if ( PS % nDimensions > 2 ) &
      call PROGRAM_HEADER % GetParameter &
             ( RP % DP_3 ( 1 : PS % nDimensions ), 'DiscontinuityPoint_3' )


    !-- Normal vector ( DP_2 - DP_1 ) x ( DP_3 - DP_1 )

    Normal ( 1 ) &
      = RP % DP_3 ( 2 ) * ( RP % DP_1 ( 3 ) - RP % DP_2 ( 3 ) ) &
          + RP % DP_1 ( 2 ) * (   RP % DP_2 ( 3 ) - RP % DP_3 ( 3 ) ) &
          + RP % DP_2 ( 2 ) * ( - RP % DP_1 ( 3 ) + RP % DP_3 ( 3 ) )
    Normal ( 2 ) &
      = RP % DP_3 ( 1 ) * ( - RP % DP_1 ( 3 ) + RP % DP_2 ( 3 ) ) &
          + RP % DP_2 ( 1 ) * (   RP % DP_1 ( 3 ) - RP % DP_3 ( 3 ) ) &
          + RP % DP_1 ( 1 ) * ( - RP % DP_2 ( 3 ) + RP % DP_3 ( 3 ) )
    Normal ( 3 ) &
      = RP % DP_3 ( 1 ) * ( RP % DP_1 ( 2 ) - RP % DP_2 ( 2 )) &
          + RP % DP_1 ( 1 ) * (   RP % DP_2 ( 2 ) - RP % DP_3 ( 2 ) ) &
          + RP % DP_2 ( 1 ) * ( - RP % DP_1 ( 2 ) + RP % DP_3 ( 2 ) )
    Normal = Normal / maxval ( Normal ) !-- to avoid overflow in the next line
    UnitNormal = Normal / sqrt ( dot_product ( Normal, Normal ) )
    call Show ( UnitNormal, 'UnitNormal', CONSOLE % INFO_3 )
  
    RP % CosTheta = dot_product ( UnitNormal, [ 0.0_KDR, 0.0_KDR, 1.0_KDR ] )
    RP % SinTheta = sqrt ( 1.0_KDR - RP % CosTheta ** 2 )
    if ( RP % SinTheta /= 0.0_KDR ) then
      RP % CosPhi &
        = dot_product ( UnitNormal, [ 1.0_KDR, 0.0_KDR, 0.0_KDR ] ) &
            / RP % SinTheta
      RP % SinPhi &
        = dot_product ( UnitNormal, [ 0.0_KDR, 1.0_KDR, 0.0_KDR ] ) &
            / RP % SinTheta
    else
      RP % CosPhi = 1.0_KDR
      RP % SinPhi = 0.0_KDR
    end if
    call Show ( 'Angles defining normal vector', CONSOLE % INFO_3 )
    call Show ( RP % SinTheta, 'SinTheta', CONSOLE % INFO_3 )
    call Show ( RP % CosTheta, 'CosTheta', CONSOLE % INFO_3 )
    call Show ( RP % SinPhi, 'SinPhi', CONSOLE % INFO_3 )
    call Show ( RP % CosPhi, 'CosPhi', CONSOLE % INFO_3 )


    !-- Initial conditions

    call SetFluid ( RP )


    !-- Cleanup

    end select !-- FA
    end select !-- PS
    end select !-- FB

  end subroutine Initialize


  subroutine Finalize ( RP )

    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    call RP % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine SetFluid ( RP )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( FB => RP % Integrator )
    class is ( FluidBoxForm )

    select type ( PS => FB % PositionSpace )
    class is ( Atlas_SC_Form )
    G => PS % Geometry ( )

    select type ( FA => FB % Current_ASC )
    class is ( Fluid_ASC_Form )

    F => FA % Fluid_P_I ( )

    call F % SetAdiabaticIndex ( RP % AdiabaticIndex )

    !-- Translate to origin, rotate normal to xz plane and then to x axis

    associate &
      ( X  => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Y  => G % Value ( :, G % CENTER_U ( 2 ) ), &
        Z  => G % Value ( :, G % CENTER_U ( 3 ) ), &
        N  => F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
        VX => F % Value ( :, F % VELOCITY_U ( 1 ) ), &
        VY => F % Value ( :, F % VELOCITY_U ( 2 ) ), &
        VZ => F % Value ( :, F % VELOCITY_U ( 3 ) ), &
        E  => F % Value ( :, F % INTERNAL_ENERGY ) )

    where ( RP % SinTheta  *  RP % CosPhi  *  ( X  -  RP % DP_1 ( 1 ) ) &
              + RP % SinTheta  *  RP % SinPhi *  ( Y  -  RP % DP_1 ( 2 ) ) &
              + RP % CosTheta  *  ( Z  -  RP % DP_1 ( 3 ) ) <= 1.e-10_KDR )     
      N   =  RP % Density_L
      E   =  RP % Energy_L 
      VX  =  RP % Speed_L  *  RP % SinTheta  *  RP % CosPhi
      VY  =  RP % Speed_L  *  RP % SinTheta  *  RP % SinPhi
      VZ  =  RP % Speed_L  *  RP % CosTheta
    elsewhere
      N   =  RP % Density_R  
      E   =  RP % Energy_R 
      VX  =  RP % Speed_R  *  RP % SinTheta  *  RP % CosPhi
      VY  =  RP % Speed_R  *  RP % SinTheta  *  RP % SinPhi
      VZ  =  RP % Speed_R  *  RP % CosTheta
    end where

    call F % ComputeFromPrimitive ( G )

    end associate !-- X, etc.
    end select !-- FA
    end select !-- PS
    end select !-- FB
    nullify ( F, G )
    
  end subroutine SetFluid

  
end module RiemannProblem_Form
