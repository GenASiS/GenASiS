module Universe_F_CC__Form

  !-- Universe_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Universe_F_C__Form

  implicit none
  private

  type, public, extends ( Universe_F_C_Form ) :: Universe_F_CC_Form
    real ( KDR ) :: &
      VelocityMax, &
      Radius_V_Max, &
      Baryons_V_Max, &
      Mass_V_Max, &
      BaryonDensity_V_Max, &
      MassDensity_V_Max, &
      BaryonDensity_C, &
      MassDensity_C, &
      Temperature_C, &
      EntropyPerBaryon_C, &
      ElectronFraction_C
  contains
    procedure, private, pass :: &
      Initialize_F_CC
    generic, public :: &
      Initialize => Initialize_F_CC
    final :: &
      Finalize
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, private, pass :: &
      InitializeAtlas
    procedure, public, pass :: &
      Compute_dT_G_CGS
  end type Universe_F_CC_Form

    private :: &
      Set_T_CheckpointInterval, &
      Compute_dT_Local

      private :: &
        Compute_dT_G_CGS_Kernel

    interface
    
      module subroutine Compute_dT_G_CGS_Kernel &
               ( dT, ProperCell, GradPhi_1, GradPhi_2, GradPhi_3, &
                 M_UU_11, M_UU_22, M_UU_33, dX_1, dX_2, dX_3, &
                 nDimensions, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), intent ( inout ) :: &
          dT
        logical ( KDL ), dimension ( : ), intent ( in ) :: &
          ProperCell
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          GradPhi_1, GradPhi_2, GradPhi_3, &
          M_UU_11, M_UU_22, M_UU_33, &
          dX_1, dX_2, dX_3
        integer ( KDI ), intent ( in ) :: &
          nDimensions
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine Compute_dT_G_CGS_Kernel

    end interface

contains


  subroutine Initialize_F_CC &
               ( U, FluidType, GravitationType, NameOption, &
                 DimensionlessOption, FinishTimeOption, RadiusMaxOption, &
                 RadiusCoreOption, RadialRatioOption, GravityFactorOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      DimensionlessOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption, &
      GravityFactorOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_CC'

    call U % Initialize_F_C &
           ( FluidType, GravitationType, &
             NameOption = NameOption, &
             DimensionlessOption = DimensionlessOption, &
             FinishTimeOption = FinishTimeOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             GravityFactorOption = GravityFactorOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nWriteOption = nWriteOption ) 

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
    I % Compute_dT_Local          =>  Compute_dT_Local
    end associate !-- I

  end subroutine Initialize_F_CC


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_CC_Form ), intent ( inout ) :: &
      U

  end subroutine Finalize


  subroutine SetBoundaryConditions ( U )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    
    select type ( I  =>  U % Integrator )
      class is ( Integrator_CS_Form )

    associate &
      ( F  =>  I % CurrentSet_X )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'OUTFLOW   ' ], iC = 1, iD = 1 )
    call F % SetBoundaryConditionsFace &
           ( [ 'REFLECTING', 'REFLECTING' ], iC = 1, iD = 2 )
    call F % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = 3 )
    end associate !-- F

    end select !-- I

  end subroutine SetBoundaryConditions


  subroutine InitializeAtlas &
               ( U, RadiusMaxOption, RadiusCoreOption, RadiusExcisionOption, &
                 RadialRatioOption, nCellsPolarOption )

      class ( Universe_F_CC_Form ), intent ( inout ) :: &
        U
      real ( KDR ), intent ( in ), optional :: &
        RadiusMaxOption, &
        RadiusCoreOption, &
        RadiusExcisionOption, &
        RadialRatioOption
      integer ( KDI ), intent ( in ), optional :: &
        nCellsPolarOption

    real ( KDR ) :: &
      RadiusMax, &
      RadiusCore, &
      RadialRatio

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_CC_Form :: I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_CC_Form )

    if ( U % Dimensionless ) then

      RadiusMax  =  10.0_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax  =  RadiusMaxOption

      RadiusCore  =  10.0_KDR / 8.0_KDR
      if ( present ( RadiusCoreOption ) ) &
        RadiusCore  =  RadiusCoreOption

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusCore = RadiusCore, &
               CommunicatorOption = PROGRAM_HEADER % Communicator, &
               NameOption = 'PositionSpace', &
               nCellsPolarOption = nCellsPolarOption )

    else

      RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadialRatio  =  5.9_KDR

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusCore = RadiusCore, &
               CommunicatorOption = PROGRAM_HEADER % Communicator, &
               NameOption = 'PositionSpace', &
               CoordinateUnitOption = U % Units_F ( 1 ) % Coordinate_PS, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end if !-- Dimensionless

    !-- Azimuthal average
    if ( PS % Chart_GS_CC % nDimensions  >  2 ) then
      allocate ( Atlas_SCG_CC_Form :: U % PositionSpace_AA )
      select type ( PS_SA  =>  U % PositionSpace_AA )
        class is ( Atlas_SCG_CC_Form )
      call PS_SA % Initialize ( PS, nDimensions = 2 )
      end select !-- PS_SA
    end if !-- nDimensions

    !-- Spherical average
    if ( PS % Chart_GS_CC % nDimensions  >  1 ) then
      allocate ( Atlas_SCG_CC_Form :: U % PositionSpace_SA )
      select type ( PS_SA  =>  U % PositionSpace_SA )
        class is ( Atlas_SCG_CC_Form )
      call PS_SA % Initialize ( PS, nDimensions = 1 )
      end select !-- PS_SA
    end if !-- nDimensions

    end select !-- PS
    end associate !-- I

  end subroutine InitializeAtlas


  subroutine Compute_dT_G_CGS ( U, dT, iC, T_Option )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    real ( KDR ), intent ( inout ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( G  =>  U % Integrator % Geometry_X )
      class is ( Gravitation_N_H_Form )
    select type ( A  =>  G % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C   =>  A % Chart_GS, &
        GV  =>  G % Storage_GS % Value )

    call Compute_dT_G_CGS_Kernel &
           ( dT, C % ProperCell, &
             GV ( :, G % POTENTIAL_GRADIENT_D_1 ), &
             GV ( :, G % POTENTIAL_GRADIENT_D_2 ), &
             GV ( :, G % POTENTIAL_GRADIENT_D_3 ), &
             GV ( :, G % METRIC_F_UU_11 ), &
             GV ( :, G % METRIC_F_UU_22 ), &
             GV ( :, G % METRIC_F_UU_33 ), &
             GV ( :, G % WIDTH_U_1 ), &
             GV ( :, G % WIDTH_U_2 ), &
             GV ( :, G % WIDTH_U_3 ), &
             C % nDimensions, &
             UseDeviceOption = G % DeviceMemory )

    end associate !-- C, etc.

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Universe_F_CC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Compute_dT_G_CGS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A

    end select !-- G

    dT  =  U % GravityFactor  *  dT
    
  end subroutine Compute_dT_G_CGS


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iR, &  !-- iRadius
      iC, &  !-- iCell
      oC, &  !-- oCell
      oI, &  !-- oIncoming
      nF
    real ( KDR ) :: &
      Constant_G, &
      T_V, &
      T_N
    real ( KDR ), dimension ( : ), allocatable :: &
       R, &
      dV, &
       N, &
       V, &
       T, &
       S, &
       Y
    real ( KDR ), dimension ( : ), pointer :: &
      T_P, &
      S_P, &
      Y_P
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO
    class ( Atlas_H_Form ), pointer :: &
      A_SA
    class ( FieldSetForm ), pointer :: &
      G_SA, &
      F_SA      

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )
    select type ( I )
      class is ( Integrator_CS_Form )

    if ( allocated ( U % PositionSpace_SA ) ) then
      A_SA  =>  U % PositionSpace_SA
      G_SA  =>  U % SA_Gravitation % FieldSet_SA
      F_SA  =>  U % SA_Fluid % FieldSet_SA
    else !-- 1D
      A_SA  =>  I % X
      G_SA  =>  I % Geometry_X
      F_SA  =>  I % CurrentSet_X
    end if

    select type ( F_SA )
      class is ( Fluid_D_Form )
    select type ( G_SA )
      class is ( Gravitation_G_Form )
    select type ( A_SA )
      class is ( Atlas_SCG_Form )
    associate &
      ( C_SA    =>  A_SA % Chart_GS, &
        G_SA_V  =>  G_SA % Storage_GS % Value, &
        F_SA_V  =>  F_SA % Storage_GS % Value, &
        M_B     =>  F_SA % BaryonMass, &
        N_Min   =>  F_SA % BaryonDensityMin )
    associate &
      ( nGL   =>  C_SA % nGhostLayers ( 1 ), &
        nC    =>  C_SA % nCells ( 1 ), &
        nCB   =>  C_SA % nCellsBrick ( 1 ), &
        nP    =>  C_SA % Communicator % Size, &
         R_P  =>  G_SA_V ( :, G_SA % CENTER_U_1 ), &
        dV_P  =>  G_SA_V ( :, G_SA % VOLUME ), &
         N_P  =>  F_SA_V ( :, F_SA % BARYON_DENSITY_C ), &
         V_P  =>  F_SA_V ( :, F_SA % VELOCITY_U_1 ) )

    nF   =   4
    T_P  =>  null ( )
    S_P  =>  null ( )
    Y_P  =>  null ( )

    select type ( F_SA )
    class is ( Fluid_P_Form )
      nF   =   nF + 2
      T_P  =>  F_SA_V ( :, F_SA % TEMPERATURE )
      S_P  =>  F_SA_V ( :, F_SA % ENTROPY_PER_BARYON )
    end select !-- F_SA

    select type ( F_SA )
    class is ( Fluid_P_HN_Form )
      nF   =   nF + 1
      Y_P  =>  F_SA_V ( :, F_SA % ELECTRON_FRACTION )
    end select !-- F_SA

    !-- Gather spherically averaged density and velocity
    !   (assume decomposition in spherical shells)

    allocate ( CO )
    call CO % Initialize &
           ( C_SA % Communicator, &
             nOutgoing = [ nF * nCB ], nIncoming = [ nF * nC ] )

    Outgoing_2D ( 1 : nCB, 1 : nF )  =>  CO % Outgoing % Value
    Outgoing_2D ( 1 : nCB, 1 )  =   R_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 2 )  =  dV_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 3 )  =   N_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 4 )  =   V_P ( nGL + 1 : nGL + nCB )
    if ( associated ( T_P ) ) &
      Outgoing_2D ( 1 : nCB, 5 )  =  T_P ( nGL + 1 : nGL + nCB )
    if ( associated ( S_P ) ) &
      Outgoing_2D ( 1 : nCB, 6 )  =  S_P ( nGL + 1 : nGL + nCB )
    if ( associated ( Y_P ) ) &
      Outgoing_2D ( 1 : nCB, 7 )  =  Y_P ( nGL + 1 : nGL + nCB )

    call CO % Gather ( )

    allocate ( R ( nC ), dV ( nC ), N ( nC ), V ( nC ) )
    if ( associated ( T_P ) ) allocate ( T ( nC ) )
    if ( associated ( S_P ) ) allocate ( S ( nC ) )
    if ( associated ( Y_P ) ) allocate ( Y ( nC ) )
    do iP  =  0,  nP - 1
      oC  =  iP * nCB
      oI  =  oC * nF
      Incoming_2D ( 1 : nCB, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCB * nF ) 
       R ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 1 )
      dV ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 2 )
       N ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 3 )
       V ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 4 )
      if ( allocated ( T ) ) &
        T ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 5 )
      if ( allocated ( S ) ) &
        S ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 6 )
      if ( allocated ( Y ) ) &
        Y ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 7 )
    end do !-- iP

    associate &
      ( Units_F    =>  U % Units_F ( 1 ), &
          V_Max    =>  U % VelocityMax, &
          R_V_Max  =>  U % Radius_V_Max, &
          B_V_Max  =>  U % Baryons_V_Max, &
          M_V_Max  =>  U % Mass_V_Max, &
          N_V_Max  =>  U % BaryonDensity_V_Max, &
        Rho_V_Max  =>  U % MassDensity_V_Max, &
          N_C      =>  U % BaryonDensity_C, &
        Rho_C      =>  U % MassDensity_C, &
          T_C      =>  U % Temperature_C, &
          S_C      =>  U % EntropyPerBaryon_C, &
          Y_C      =>  U % ElectronFraction_C )

    !-- VelocityMax

    V_Max  =  maxval ( abs ( V ) )

    iR  =  nC
    do iC  =  nC, 1, -1
      if ( N ( iC )  >  1.01_KDR  *  N_Min ) then
        if ( abs ( V ( iC ) )  ==  V_Max ) then
          iR  =  iC
          exit
        end if
      end if
    end do !-- iC

      R_V_Max  =  R ( iR )
      B_V_Max  =  sum ( N ( : iR )  *  dV ( : iR ) )
      M_V_Max  =  M_B * B_V_Max
      N_V_Max  =  B_V_Max  /  sum ( dV ( : iR ) )
    Rho_V_Max  =  M_B * N_V_Max

    !-- Center

      N_C  =  N ( 1 )
    Rho_C  =  M_B  *  N_C
    if ( allocated ( S ) )  S_C  =  S ( 1 )
    if ( allocated ( T ) )  T_C  =  T ( 1 )
    if ( allocated ( Y ) )  Y_C  =  Y ( 1 )

    !-- Time scales and CheckpointTimeInterval

    T_V  =  R_V_Max  /  max ( V_Max, sqrt ( tiny ( 0.0_KDR ) ) )

    if ( U % Dimensionless ) then
      Constant_G  =  1.0_KDR
    else
      Constant_G  =  CONSTANT % GRAVITATIONAL
    end if

    T_N  =  ( Constant_G * M_B * min ( N_V_Max, N_C ) ) ** ( -0.5_KDR )

    I % T_CheckpointInterval  =  min ( T_V, T_N )  /  I % nWrite

    !-- Display

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( V_Max, Units_F % Velocity_U ( 1 ), &
                'VelocityMax', I % IGNORABILITY )
    call Show ( R_V_Max, Units_F % Coordinate_PS ( 1 ), &
                'Radius_V_Max', I % IGNORABILITY )
    call Show ( B_V_Max, Units_F % Number, &
                'Baryons_V_Max', I % IGNORABILITY )
    call Show ( M_V_Max, Units_F % Mass, &
                'Mass_V_Max', I % IGNORABILITY )
    call Show ( N_V_Max, Units_F % NumberDensity, &
                'BaryonDensity_V_Max', I % IGNORABILITY )
    call Show ( Rho_V_Max, Units_F % MassDensity, &
                'MassDensity_V_Max', I % IGNORABILITY )
    call Show ( N_C, Units_F % NumberDensity, &
                'BaryonDensity_C', I % IGNORABILITY )
    call Show ( Rho_C, Units_F % MassDensity, &
                'MassDensity_C', I % IGNORABILITY )
    if ( allocated ( T ) ) &
      call Show ( T_C, Units_F % Temperature, &
                  'Temperature_C', I % IGNORABILITY )
    if ( allocated ( S ) ) &
      call Show ( S_C, Units_F % EnergyDensity  /  Units_F % NumberDensity  &
                       /  Units_F % Temperature, &
                  'EntropyPerBaryon_C', I % IGNORABILITY )
    if ( allocated ( Y ) ) &
      call Show ( Y_C, &
                  'ElectronFraction_C', I % IGNORABILITY )
    call Show ( T_V, I % Unit_T, &
                'T_Velocity', I % IGNORABILITY )
    call Show ( T_N, I % Unit_T, &
                'T_Density', I % IGNORABILITY )

    !-- Cleanup

    end associate !-- Units_F, etc.
    end associate !-- nGL, etc.
    end associate !-- C_SA
    end select !-- A_SA
    end select !-- G_SA
    end select !-- F_SA
    end select !-- I
    end select !-- U

  end subroutine Set_T_CheckpointInterval


  subroutine Compute_dT_Local ( I, dT_Candidate, iC, T_Option )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      dT_Candidate
    integer ( KDI ), intent ( in ) :: &
      iC
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
    class is ( Universe_F_C_Form )
      if ( U % Coarsen ) then
        call U % Compute_dT_CS_CGS_C ( dT_Candidate ( 1 ), iC, T_Option )
      else !-- .not. Coarsen
        select type ( I )
        class is ( Integrator_CS_Form )
          call I % Compute_dT_CS_CGS ( dT_Candidate ( 1 ), iC, T_Option )
        end select !-- I
      end if !-- Coarsen
    end select !-- U

    select type ( U  =>  I % System )
    class is ( Universe_F_CC_Form )
      call U % Compute_dT_G_CGS ( dT_Candidate ( 2 ), iC, T_Option )
    end select !-- U

  end subroutine Compute_dT_Local


end module Universe_F_CC__Form
