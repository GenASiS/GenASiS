module RiemannProblem_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( FluidBoxForm ) :: RiemannProblemForm
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
    procedure, private, pass :: &
      Initialize_RP
    generic, public :: &
      Initialize => Initialize_RP
    final :: &
      Finalize
  end type RiemannProblemForm

    private :: &
      InitializeFluidBox, &
      SetProblem

    private :: &
      SetFluid

        private :: &
          SetFluidKernel
    
contains


  subroutine Initialize_RP ( RP, Name )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    character ( * ), intent ( in ) :: &
      Name

    if ( RP % Type == '' ) &
      RP % Type = 'a RiemannProblem'

    call InitializeFluidBox ( RP, Name )
    call SetProblem ( RP )

  end subroutine Initialize_RP


  subroutine Finalize ( RP )

    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP

  end subroutine Finalize


  subroutine InitializeFluidBox ( RP, Name )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD  !-- iDimension

    allocate ( RP % BoundaryConditionsFace ( 3 ) )
    associate ( BCF => RP % BoundaryConditionsFace )
    do iD = 1, 3
      call BCF ( iD ) % Initialize ( [ 'REFLECTING', 'REFLECTING' ] )     
    end do
    end associate !-- BCF
    
    call RP % Initialize &
           ( FluidType = 'IDEAL', GeometryType = 'GALILEAN', Name = Name, &
             nCellsOption = [ 128, 128, 128 ] )

  end subroutine InitializeFluidBox


  subroutine SetProblem ( RP )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP

    real ( KDR ), dimension ( 3 ) :: &
      Normal, &
      UnitNormal
    class ( GeometryFlatForm ), pointer :: &
      G
    class ( Fluid_P_I_Form ), pointer :: &
      F

    select type ( I => RP % Integrator )
    class is ( Integrator_C_PS_Form )

    select type ( FA => I % Current_ASC )
    class is ( Fluid_ASC_Form )

    select type ( PS => I % PositionSpace )
    class is ( Atlas_SC_Form )

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
    call PROGRAM_HEADER % GetParameter ( RP % AdiabaticIndex, &
                                         'AdiabaticIndex' )

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

    G => PS % Geometry ( )
    F => FA % Fluid_P_I ( )
    call SetFluid ( RP, F, G )

    end select !-- PS
    end select !-- FA
    end select !-- I
    nullify ( G, F )

  end subroutine SetProblem


  subroutine SetFluid ( RP, F, G )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F
    class ( GeometryFlatForm ), intent ( in ) :: &
      G

    call F % SetAdiabaticIndex ( RP % AdiabaticIndex )

    call SetFluidKernel &
           (      X = G % Value ( :, G % CENTER_U ( 1 ) ), &
                  Y = G % Value ( :, G % CENTER_U ( 2 ) ), &
                  Z = G % Value ( :, G % CENTER_U ( 3 ) ), &
               DP_1 = RP % DP_1, &
             sTheta = RP % SinTheta, &
             cTheta = RP % CosTheta, &
               sPhi = RP % SinPhi, &
               cPhi = RP % CosPhi, &
                N_L = RP % Density_L, &
                N_R = RP % Density_R, &
                E_L = RP % Energy_L, &
                E_R = RP % Energy_R, &
                V_L = RP % Speed_L, &
                V_R = RP % Speed_R, &
                  N = F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
                  E = F % Value ( :, F % INTERNAL_ENERGY ), &
                 VX = F % Value ( :, F % VELOCITY_U ( 1 ) ), &
                 VY = F % Value ( :, F % VELOCITY_U ( 2 ) ), &
                 VZ = F % Value ( :, F % VELOCITY_U ( 3 ) ) )

  end subroutine SetFluid


  subroutine SetFluidKernel &
               ( X, Y, Z, DP_1, sTheta, cTheta, sPhi, cPhi, N_L, N_R, &
                 E_L, E_R, V_L, V_R, N, E, VX, VY, VZ )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      X, Y, Z
    real ( KDR ), dimension ( 3 ), intent ( in ) :: &
      DP_1
    real ( KDR ), intent ( in ) :: &
      sTheta, cTheta, &
      sPhi, cPhi, &
      N_L, N_R, &
      E_L, E_R, &
      V_L, V_R
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      N, &
      E, &
      VX, VY, VZ

!     !-- Translate to origin, rotate normal to xz plane and then to x axis

    where (     sTheta  *  cPhi  *  ( X  -  DP_1 ( 1 ) ) &
              + sTheta  *  sPhi  *  ( Y  -  DP_1 ( 2 ) ) &
              + cTheta  *           ( Z  -  DP_1 ( 3 ) ) &
           <=  1.e-10_KDR )     
      N   =  N_L
      E   =  E_L 
      VX  =  V_L  *  sTheta  *  cPhi
      VY  =  V_L  *  sTheta  *  sPhi
      VZ  =  V_L  *  cTheta
    elsewhere
      N   =  N_R  
      E   =  E_R 
      VX  =  V_R  *  sTheta  *  cPhi
      VY  =  V_R  *  sTheta  *  sPhi
      VZ  =  V_R  *  cTheta
    end where

  end subroutine SetFluidKernel


end module RiemannProblem_Form
