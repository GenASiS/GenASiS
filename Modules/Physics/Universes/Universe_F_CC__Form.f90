module Universe_F_CC__Form

  !-- Universe_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use Measures_F_CC__Form
  use Series_F_CC__Form
  use Universe_F_C__Form

  implicit none
  private

  type, public, extends ( Universe_F_C_Form ) :: Universe_F_CC_Form
    class ( Measures_F_CC_Form ), allocatable :: &
      Measures
  contains
    procedure, public, pass :: &
      Initialize_F_CC
    generic, public :: &
      Initialize => Initialize_F_CC
    final :: &
      Finalize
    procedure, public, pass :: &
      SetBoundaryConditions
    procedure, public, pass :: &
      SetMeasures
    procedure, public, pass :: &
      InitializeAtlas
    procedure, public, pass :: &
      Compute_dT_G_CGS
    procedure, public, nopass :: &
      InitializeSeries_CC      
    procedure, public, nopass :: &
      Set_T_CheckpointInterval_F_CC
  end type Universe_F_CC_Form

    private :: &
      Analyze_F_CC, &
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
               ( U, FluidType, GravitationType, Name, &
                 UnitsTypeOption, FinishTimeOption, RadiusMaxOption, &
                 RadiusCoreOption, RadialRatioOption, GravityFactorOption, &
                 nCellsPolarOption, nWriteOption )

    class ( Universe_F_CC_Form ), intent ( inout ), target :: &
      U
    character ( * ), intent ( in )  :: &
      FluidType, &
      GravitationType, &
      Name
    character ( * ), intent ( in ), optional :: &
      UnitsTypeOption
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
           ( FluidType, GravitationType, Name, &
             UnitsTypeOption = UnitsTypeOption, &
             FinishTimeOption = FinishTimeOption, &
             RadiusMaxOption = RadiusMaxOption, &
             RadiusCoreOption = RadiusCoreOption, &
             RadialRatioOption = RadialRatioOption, &
             GravityFactorOption = GravityFactorOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nWriteOption = nWriteOption ) 

    call U % SetMeasures ( )

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Compute_dT_Local          =>  Compute_dT_Local
    I % InitializeSeries          =>  InitializeSeries_CC
    I % Analyze                   =>  Analyze_F_CC
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval_F_CC
    end associate !-- I

  end subroutine Initialize_F_CC


  impure elemental subroutine Finalize ( U )

    type ( Universe_F_CC_Form ), intent ( inout ) :: &
      U

    if ( allocated ( U % Measures ) ) &
      deallocate ( U % Measures )

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


  subroutine SetMeasures ( U )

    class ( Universe_F_CC_Form ), intent ( inout ), target :: &
      U

    class ( Atlas_H_Form ), pointer :: &
      A_SA, &
      A
    class ( FieldSet_BM_Form ), pointer :: &
      G_SA, &
      F_SA, &
      F

    if ( allocated ( U % PositionSpace_SA ) ) then
      A_SA  =>  U % PositionSpace_SA
      G_SA  =>  U % SA_Gravitation % FieldSet_SA
      F_SA  =>  U % SA_Fluid % FieldSet_SA
      F     =>  U % SA_Fluid % FieldSet
    else !-- 1D
      select type ( I  =>  U % Integrator )
        class is ( Integrator_CS_Form )
      A_SA  =>  I % X
      G_SA  =>  I % Geometry_X
      F_SA  =>  I % CurrentSet_X
      F     =>  I % CurrentSet_X
      end select !-- I
    end if

    allocate ( U % Measures )
    associate ( M  =>  U % Measures )
    call M % Initialize ( F, F_SA, G_SA, A_SA, Units_F = U % Units_F ( 1 ) )
    end associate !-- M

  end subroutine SetMeasures


  subroutine InitializeAtlas &
               ( U, CommunicatorOption, RadiusMaxOption, RadiusCoreOption, &
                 RadiusExcisionOption, RadialRatioOption, nCellsPolarOption )

    class ( Universe_F_CC_Form ), intent ( inout ) :: &
      U
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
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
    type ( QuantityForm ), dimension ( 3 ) :: &
      CoordinateUnit
    type ( CommunicatorForm ), pointer :: &
      Communicator

    if ( present ( CommunicatorOption ) ) then
      Communicator  =>  CommunicatorOption
    else
      Communicator  =>  U % Communicator
    end if

    associate ( I  =>  U % Integrator )

    allocate ( Atlas_SCG_CC_Form :: I % X )
    select type ( PS  =>  I % X )
      class is ( Atlas_SCG_CC_Form )

    select case ( trim ( U % UnitsType ) )
    case ( 'ASTROPHYSICS' )

      CoordinateUnit  =  [ UNIT % KILOMETER, UNIT % RADIAN, UNIT % RADIAN ]

      RadiusMax    =  1.0e4_KDR  *  UNIT % KILOMETER
      RadiusCore   =   16.0_KDR  *  UNIT % KILOMETER
!      RadiusCore   =   32.0_KDR  *  UNIT % KILOMETER

!-- 40 nCellsCore for 128 nCellsPolar, polar/radial aspect ratio close to 1
      RadialRatio  =  2.4_KDR  !-- RadiusCore = 16.0 km

!-- 100 nCellsCore for 128 nCellsPolar, polar/radial aspect ratio close to 2.5  
!      RadialRatio  =  5.9_KDR  !-- RadiusCore = 16.0 km
!      RadialRatio  =  5.3_KDR  !-- RadiusCore = 32.0 km

      if ( present ( RadiusMaxOption ) ) &
        RadiusMax  =  RadiusMaxOption
      if ( present ( RadiusCoreOption ) ) &
        RadiusCore  =  RadiusCoreOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio  =  RadialRatioOption

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusCore = RadiusCore, &
               CommunicatorOption = Communicator, &
               NameOption = 'PositionSpace', &
               DeviceMemoryOption = U % DeviceMemory, &
               CoordinateUnitOption = CoordinateUnit, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    case default

      RadiusMax   =  10.0_KDR
      RadiusCore  =  10.0_KDR / 8.0_KDR
      RadialRatio =  2.45_KDR
      if ( present ( RadiusMaxOption ) ) &
        RadiusMax  =  RadiusMaxOption
      if ( present ( RadiusCoreOption ) ) &
        RadiusCore  =  RadiusCoreOption
      if ( present ( RadialRatioOption ) ) &
        RadialRatio  =  RadialRatioOption

      call PS % Initialize &
             ( RadiusMax = RadiusMax, &
               RadiusCore = RadiusCore, &
               CommunicatorOption = Communicator, &
               NameOption = 'PositionSpace', &
               DeviceMemoryOption = U % DeviceMemory, &
               RadialRatioOption = RadialRatio, &
               nCellsPolarOption = nCellsPolarOption )

    end select !-- UnitsType

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

  end subroutine Compute_dT_G_CGS


  subroutine InitializeSeries_CC ( I )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I

    allocate ( Series_F_CC_Form :: I % Series )

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )
    select type ( I )
      class is ( Integrator_CS_Form )
    select type ( S  =>  I % Series )
      class is ( Series_F_CC_Form )
    call S % Initialize &
      ( U % Measures, I % CurrentSet_X, I % GridImageStream, I % dT_Label, &
        I % Unit_T, I % dT_Candidate, I % T, I % Communicator % Rank, &
        I % nWrite, I % iCycle )
    end select !-- S
    end select !-- I
    end select !-- U

  end subroutine InitializeSeries_CC


  subroutine Set_T_CheckpointInterval_F_CC ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    real ( KDR ) :: &
      Constant_G, &
      T_V_Max, &
      T_V_Shock, &
      T_N

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )
    associate &
      ( M  =>  U % Measures )

    associate &
      (     V_Max  =>  M % VelocityMax, &
          R_V_Max  =>  M % Radius_V_Max, &
        Rho_V_Max  =>  M % MassDensity_V_Max, &
        Rho_C      =>  M % MassDensity_C, &
          R_Shock  =>  M % RadiusShock )

    !-- Time scales and CheckpointTimeInterval

    T_V_Max    =  R_V_Max  /  max ( V_Max, sqrt ( tiny ( 0.0_KDR ) ) )

    T_V_Shock  =  R_Shock  &
                  /  max ( V_Max, sqrt ( tiny ( 0.0_KDR ) ) )
    !-- Note V_Max is still used even though it may not correspond to the 
    !   velocity near R_Shock when this differs from R_V_Max.

    if ( trim ( U % UnitsType )  ==  '' ) then
      Constant_G  =  1.0_KDR
    else
      Constant_G  =  CONSTANT % GRAVITATIONAL
    end if

    T_N  =  ( Constant_G * min ( Rho_V_Max, Rho_C ) ) ** ( -0.5_KDR )

    if ( R_Shock  >  R_V_Max ) then
      I % T_CheckpointInterval  =  T_V_Shock  /  I % nWrite
    else
      I % T_CheckpointInterval  =  min ( T_V_Max, T_N )  /  I % nWrite
    end if

    !-- Display

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( T_N, I % Unit_T, 'T_Density',  I % IGNORABILITY )
    call Show ( T_V_Max, I % Unit_T, 'T_Velocity_Max', I % IGNORABILITY )
    call Show ( T_V_Shock, I % Unit_T, 'T_Velocity_Shock', I % IGNORABILITY )

    !-- Cleanup

    end associate !-- V_Max, etc.
    end associate !-- M
    end select !-- U

  end subroutine Set_T_CheckpointInterval_F_CC


  subroutine Analyze_F_CC ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )

    call U % Average_F_C ( )

    associate ( M  =>  U % Measures )
    call M % Compute ( )
    end associate !-- M

    end select !-- U

    select type ( I )
      class is ( Integrator_CS_Form )
    call I % Analyze_CS ( I, Ignorability, T_Option )
    end select !-- I

  end subroutine Analyze_F_CC


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
      class is ( Universe_F_CC_Form )
    select type ( I )
      class is ( Integrator_CS_Form )
    associate &
      ( dT_1  =>  dT_Candidate ( 1 ), &
        dT_2  =>  dT_Candidate ( 2 ) )

    !-- Advection step

    if ( U % Coarsen ) then
      call U % Compute_dT_CS_CGS_C &
             ( I % EigenspeedSet_X, dT_1, U % Coarsening_F, iC, T_Option )
    else !-- .not. Coarsen
      call I % Compute_dT_CS_CGS &
             ( I % EigenspeedSet_X, dT_1, iC, T_Option )
    end if !-- Coarsen
    dT_1  =  I % CourantFactor  *  dT_1
    
    !-- Gravity step

    call U % Compute_dT_G_CGS ( dT_2, iC, T_Option )
    dT_2  =  U % GravityFactor  *  dT_2    

    end associate !-- dT_1, etc.
    end select !-- I
    end select !-- U

  end subroutine Compute_dT_Local


end module Universe_F_CC__Form
