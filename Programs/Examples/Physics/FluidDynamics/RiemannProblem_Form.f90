module RiemannProblem_Form

  use GenASiS

  implicit none
  private

  type, public, extends ( Universe_F_B_Form ) :: RiemannProblemForm
    real ( KDR ) :: &
      Density_L, Density_R, &  !-- Left and Right
      Pressure_L, Pressure_R, &
      Energy_L, Energy_R, &
      Speed_L, Speed_R, &
      AdiabaticIndex, &
      SinTheta, CosTheta, &
      SinPhi, CosPhi
    real ( KDR ), dimension ( 3 ) :: &
      DP_1, DP_2, DP_3, &  !-- DiscontinuityPoint_1, etc.
      UnitNormal
  contains
    procedure, private, pass :: &
      Initialize_H
    final :: &
      Finalize
    procedure, public, pass :: &
      ShowParameters
  end type RiemannProblemForm

    private :: &
      InitializeUniverse, &
      SetInitial

      private :: &
        SetFluid

        private :: &
          SetFluidKernel
    
contains


  subroutine Initialize_H ( U, NameOption )

    class ( RiemannProblemForm ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in ), optional  :: &
      NameOption

    character ( LDL ) :: &
      Name

    if ( U % Type  ==  '' ) &
      U % Type  =  'a RiemannProblem'

    Name  =  'RiemannProblem'
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call InitializeUniverse ( U, Name )

  end subroutine Initialize_H


  subroutine Finalize ( RP )

    type ( RiemannProblemForm ), intent ( inout ) :: &
      RP

  end subroutine Finalize


  subroutine ShowParameters ( U )

    class ( RiemannProblemForm ), intent ( in ) :: &
      U

    call U % Universe_F_B_Form % ShowParameters ( )

    call Show ( U % Density_L,      'Density_L' )
    call Show ( U % Pressure_L,     'Pressure_L' )
    call Show ( U % Speed_L,        'Speed_L' )
    call Show ( U % Density_R,      'Density_R' )
    call Show ( U % Pressure_R,     'Pressure_R' )
    call Show ( U % Speed_R,        'Speed_R' )
    call Show ( U % AdiabaticIndex, 'AdiabaticIndex' )
    call Show ( U % DP_1,           'DiscontinuityPoint_1' )
    call Show ( U % DP_2,           'DiscontinuityPoint_2' )
    call Show ( U % DP_3,           'DiscontinuityPoint_3' )
    call Show ( U % UnitNormal,     'UnitNormal' )
    call Show ( U % SinTheta,       'SinTheta' )
    call Show ( U % CosTheta,       'CosTheta' )
    call Show ( U % SinPhi,         'SinPhi' )
    call Show ( U % CosPhi,         'CosPhi' )

  end subroutine ShowParameters


  subroutine InitializeUniverse ( RP, Name )

    class ( RiemannProblemForm ), intent ( inout ), target :: &
      RP
    character ( * ), intent ( in )  :: &
      Name

    integer ( KDI ) :: &
      iD

    call RP % Initialize &
           ( FluidType = 'IDEAL', &
             GravitationType = 'GALILEO', &
             NameOption = Name, &
             nCellsOption = [ 128, 128, 128 ] )

    select type ( I  =>  RP % Integrator )
      class is ( Integrator_CS_Form )
    associate &
      ( F  =>  I % CurrentSet_X )
    do iD  =  1, 3
      call F % SetBoundaryConditionsFace &
             ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = iD )
    end do !-- iD
    end associate !-- F
    end select !-- I
             
    RP % Integrator % SetInitial  =>  SetInitial
    RP % Integrator % System      =>  RP

  end subroutine InitializeUniverse


  subroutine SetInitial ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    real ( KDR ), dimension ( 3 ) :: &
      Normal

    select type ( RP  =>  I % System )
      class is ( RiemannProblemForm )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( F  =>  I % CurrentSet_X )
      class is ( Fluid_P_I_Form )
    select type ( A  =>  F % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

    RP % Density_L       =  1.0_KDR
    RP % Pressure_L      =  1.0_KDR
    RP % Speed_L         =  0.0_KDR
    RP % Density_R       =  0.125_KDR
    RP % Pressure_R      =  0.1_KDR
    RP % Speed_R         =  0.0_KDR
    RP % AdiabaticIndex  =  1.4_KDR

    call PROGRAM_HEADER % GetParameter ( RP % Density_L,  'Density_L' )
    call PROGRAM_HEADER % GetParameter ( RP % Pressure_L, 'Pressure_L' )
    call PROGRAM_HEADER % GetParameter ( RP % Speed_L,    'Speed_L' )
    call PROGRAM_HEADER % GetParameter ( RP % Density_R,  'Density_R' )
    call PROGRAM_HEADER % GetParameter ( RP % Pressure_R, 'Pressure_R' )
    call PROGRAM_HEADER % GetParameter ( RP % Speed_R,    'Speed_R' )
    call PROGRAM_HEADER % GetParameter ( RP % AdiabaticIndex, &
                                         'AdiabaticIndex' )

    RP % Energy_L  =  RP % Pressure_L  /  ( RP % AdiabaticIndex - 1.0_KDR )
    RP % Energy_R  =  RP % Pressure_R  /  ( RP % AdiabaticIndex - 1.0_KDR )


    !-- Three points define the plane of discontinuity

    RP % DP_1  =  [ 0.5_KDR, 0.0_KDR, 0.0_KDR ]
    RP % DP_2  =  [ 0.0_KDR, 0.5_KDR, 0.0_KDR ]
    RP % DP_3  =  [ 0.0_KDR, 0.0_KDR, 0.5_KDR ]

    if ( C % nDimensions  <  3 )  &
      RP % DP_3 ( 3 )  =  0.1 * sqrt ( huge ( 1.0_KDR ) )
    if ( C % nDimensions  <  2 )  &
      RP % DP_2 ( 2 )  =  0.1 * sqrt ( huge ( 1.0_KDR ) )

    call PROGRAM_HEADER % GetParameter &
           ( RP % DP_1 ( 1 : C % nDimensions ), 'DiscontinuityPoint_1' )
    if ( C % nDimensions  >  1 ) &
      call PROGRAM_HEADER % GetParameter &
             ( RP % DP_2 ( 1 : C % nDimensions ), 'DiscontinuityPoint_2' )
    if ( C % nDimensions  >  2 ) &
      call PROGRAM_HEADER % GetParameter &
             ( RP % DP_3 ( 1 : C % nDimensions ), 'DiscontinuityPoint_3' )


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
    Normal  =  Normal / maxval ( Normal ) !-- to avoid overflow in the next line
    RP % UnitNormal  =  Normal / sqrt ( dot_product ( Normal, Normal ) )
  
    RP % CosTheta  =  &
      dot_product ( RP % UnitNormal, [ 0.0_KDR, 0.0_KDR, 1.0_KDR ] )
    RP % SinTheta  =  &
      sqrt ( 1.0_KDR  -  RP % CosTheta ** 2 )
    if ( RP % SinTheta /= 0.0_KDR ) then
      RP % CosPhi &
        =  dot_product ( RP % UnitNormal, [ 1.0_KDR, 0.0_KDR, 0.0_KDR ] ) &
             / RP % SinTheta
      RP % SinPhi &
        =  dot_product ( RP % UnitNormal, [ 0.0_KDR, 1.0_KDR, 0.0_KDR ] ) &
             / RP % SinTheta
    else
      RP % CosPhi  =  1.0_KDR
      RP % SinPhi  =  0.0_KDR
    end if
    call SetFluid ( RP, F )

    end associate !-- C
    end select !-- A
    end select !-- F
    end select !-- I
    end select !-- RP

  end subroutine SetInitial


  subroutine SetFluid ( RP, F )

    class ( RiemannProblemForm ), intent ( inout ) :: &
      RP
    class ( Fluid_P_I_Form ), intent ( inout ) :: &
      F

    call F % SetAdiabaticIndex ( RP % AdiabaticIndex )

    associate &
      ( G  =>  F % Geometry )
    associate &
      ( FV  =>  F % Storage_GS % Value, &
        GV  =>  G % Storage_GS % Value )

    call SetFluidKernel &
           (      X = GV ( :, G % CENTER_U ( 1 ) ), &
                  Y = GV ( :, G % CENTER_U ( 2 ) ), &
                  Z = GV ( :, G % CENTER_U ( 3 ) ), &
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
                  N = FV ( :, F % BARYON_DENSITY_C ), &
                  E = FV ( :, F % ENERGY_DENSITY_C ), &
                 VX = FV ( :, F % VELOCITY_U_1 ), &
                 VY = FV ( :, F % VELOCITY_U_2 ), &
                 VZ = FV ( :, F % VELOCITY_U_3 ) )

    end associate !-- FV, etc.
    end associate !-- G

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
