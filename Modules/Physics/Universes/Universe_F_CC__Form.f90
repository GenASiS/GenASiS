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
               ( U, FluidType, GravitationType, NameOption, FinishTimeOption, &
                 RadiusMaxOption, RadiusCoreOption, RadialRatioOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), intent ( in ), optional :: &
      FinishTimeOption, &
      RadiusMaxOption, &
      RadiusCoreOption, &
      RadialRatioOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nWriteOption

    if ( U % Type == '' ) &
      U % Type = 'a Universe_F_CC'

    call U % Initialize_F_C &
           ( FluidType, GravitationType, &
             NameOption = NameOption, &
             FinishTimeOption = FinishTimeOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nWriteOption = nWriteOption ) 

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

    allocate ( U % Units_F ( 1 ) )

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
             NameOption = 'PositionSpace' )

    allocate ( Atlas_SCG_CC_Form :: U % PositionSpace_SA )
    select type ( PS_SA  =>  U % PositionSpace_SA )
      class is ( Atlas_SCG_CC_Form )
    call PS_SA % Initialize ( PS )
    end select !-- PS_SA

    ! if ( FC % Dimensionless ) then

    !   call PS % CreateChart_CC ( )

    ! else

    !   RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
    !   RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
    !   RadialRatio  =  2.4_KDR

    !   call PS % CreateChart_CC &
    !          ( CoordinateUnitOption = FC % Units % Coordinate_PS, &
    !            RadiusCoreOption = RadiusCore, &
    !            RadiusMaxOption = RadiusMax, &
    !            RadialRatioOption = RadialRatio, &
    !            nCellsPolarOption = nCellsPolarOption )

    !   FC % RadiusPolarMomentum  =  8.0_KDR  *  UNIT % KILOMETER
    !   call PROGRAM_HEADER % GetParameter &
    !          ( FC % RadiusPolarMomentum, 'RadiusPolarMomentum' )

    ! end if !-- Dimensionless

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

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    integer ( KDI ) :: &
      iP, &  !-- iProcess
      iR, &  !-- iRadius
      iC, &  !-- iCell
      oC, &  !-- oCell
      oI     !-- oIncoming
    real ( KDR ) :: &
      Constant_G, &
      D_Min, &
      V_Max, &
      R_V_Max, &
      B_V_Max, &
      N_V_Max, &
      N_Max, &
      T_V, &
      T_N
    real ( KDR ), dimension ( : ), allocatable :: &
       R, &
      dV, &
       N, &
       V
    real ( KDR ), dimension ( :, : ), pointer :: &
      Outgoing_2D, &
      Incoming_2D
    type ( CollectiveOperation_R_Form ), allocatable :: &
      CO

    select type ( U  =>  I % System )
      class is ( Universe_F_C_Form )
    select type ( F_SA  =>  U % SA_Fluid % FieldSet_SA )
      class is ( Fluid_D_Form )
    select type ( G_SA  =>  U % SA_Gravitation % FieldSet_SA )
      class is ( Gravitation_G_Form )
    select type ( A_SA  =>  F_SA % Atlas )
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
        nF    =>  4, &
         R_P  =>  G_SA_V ( :, G_SA % CENTER_U_1 ), &
        dV_P  =>  G_SA_V ( :, G_SA % VOLUME ), &
         N_P  =>  F_SA_V ( :, F_SA % BARYON_DENSITY_C ), &
         V_P  =>  F_SA_V ( :, F_SA % VELOCITY_U_1 ) )

    !-- Gather spherically averaged density and velocity
    !   (assume decomposition in spherical shells)

    allocate ( CO )
    call CO % Initialize &
           ( C_SA % Communicator, &
             nOutgoing = [ nF * nCB ], nIncoming = [ nF * nC ] )

    oC  =  C_SA % nGhostLayers ( 1 )
    Outgoing_2D ( 1 : nCB, 1 : nF )  =>  CO % Outgoing % Value
    Outgoing_2D ( 1 : nCB, 1 )      =    R_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 2 )      =   dV_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 3 )      =    N_P ( nGL + 1 : nGL + nCB )
    Outgoing_2D ( 1 : nCB, 4 )      =    V_P ( nGL + 1 : nGL + nCB )

    call CO % Gather ( )

    allocate ( R ( nC ), dV ( nC ), N ( nC ), V ( nC ) )
    do iP  =  0,  nP - 1
      oC  =  iP * nCB
      oI  =  oC * nF
      Incoming_2D ( 1 : nCB, 1 : nF )  &
        =>  CO % Incoming % Value ( oI + 1 : oI + nCB * nF ) 
       R ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 1 )
      dV ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 2 )
       N ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 3 )
       V ( oC + 1 : oC + nCB )  =  Incoming_2D ( 1 : nCB, 4 )
    end do !-- iP

    !-- Velocity

    V_Max  =  maxval ( abs ( V ) )

    iR  =  nC
    do iC  =  nC, 1, -1
      if ( N ( iC )  >  1.01_KDR  *  N_Min ) then
        if ( V ( iC )  ==  V_Max ) then
          iR  =  iC
          exit
        end if
      end if
    end do !-- iC

    R_V_Max  =  R ( iR )

    !-- Density

    B_V_Max  =  sum ( N ( : iR )  *  dV ( : iR ) )

    N_V_Max  =  B_V_Max  /  sum ( dV ( : iR ) )

    N_Max  =  maxval ( N )

    !-- Time scales and CheckpointTimeInterval

    T_V  =  R_V_Max  /  max ( V_Max, sqrt ( tiny ( 0.0_KDR ) ) )

    !-- FIXME: nontrivial units
    Constant_G  =  1.0_KDR

    T_N  =  ( Constant_G * M_B * min ( N_V_Max, N_Max ) ) ** ( -0.5_KDR )

    I % T_CheckpointInterval  =  min ( T_V, T_N )  /  I % nWrite

    !-- Display

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( V_Max, F_SA % Unit ( F_SA % VELOCITY_U_1, 1 ), &
                'VelocityMax', I % IGNORABILITY )
    call Show ( R_V_Max, G_SA % Unit ( G_SA % CENTER_U_1, 1 ), &
                'RadiusVelocityMax', I % IGNORABILITY )
    call Show ( B_V_Max, &
                'BaryonsVelocityMax', I % IGNORABILITY )
    call Show ( N_V_Max, F_SA % Unit ( F_SA % BARYON_DENSITY_C, 1 ), &
                'DensityVelocityMax', I % IGNORABILITY )
    call Show ( N_Max, F_SA % Unit ( F_SA % BARYON_DENSITY_C, 1 ), &
                'DensityMax', I % IGNORABILITY )
    ! call Show ( NumberDensityAve, FC % Units % NumberDensity, &
    !             'NumberDensityAve', I % IGNORABILITY )
    ! call Show ( NumberDensityMax, FC % Units % NumberDensity, &
    !             'NumberDensityMax', I % IGNORABILITY )
    ! call Show ( BaryonMass * NumberDensityAve, FC % Units % MassDensity, &
    !             'MassDensityAve', I % IGNORABILITY )
    ! call Show ( BaryonMass * NumberDensityMax, FC % Units % MassDensity, &
    !             'MassDensityMax', I % IGNORABILITY )
    call Show ( T_V, I % Unit_T, &
                'T_Velocity', I % IGNORABILITY )
    call Show ( T_N, I % Unit_T, &
                'T_Density', I % IGNORABILITY )

    !-- Cleanup

    deallocate ( V, N, dV, R )
    deallocate ( CO )

    end associate !-- nC, etc.
    end associate !-- C_SA
    end select !-- A_SA
    end select !-- G_SA
    end select !-- F_SA
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

    select type ( I )
      class is ( Integrator_CS_Form )
    call I % Compute_dT_CGS ( dT_Candidate ( 1 ), iC, T_Option )
    end select !-- I

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )
    call U % Compute_dT_G_CGS ( dT_Candidate ( 2 ), iC, T_Option )
    end select !-- U

  end subroutine Compute_dT_Local


end module Universe_F_CC__Form
