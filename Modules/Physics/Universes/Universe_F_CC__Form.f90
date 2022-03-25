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
    procedure, public, nopass :: &
      Analyze_F_CC
  end type Universe_F_CC_Form

    private :: &
      Set_T_CheckpointInterval, &
      Compute_dT_Local, &
      InitializeSeries

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

    class ( Universe_F_CC_Form ), intent ( inout ), target :: &
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

    class ( Atlas_H_Form ), pointer :: &
      A_SA
    class ( FieldSetForm ), pointer :: &
      G_SA, &
      F_SA      

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

    !-- Measures

    if ( allocated ( U % PositionSpace_SA ) ) then
      A_SA  =>  U % PositionSpace_SA
      G_SA  =>  U % SA_Gravitation % FieldSet_SA
      F_SA  =>  U % SA_Fluid % FieldSet_SA
    else !-- 1D
      select type ( I  =>  U % Integrator )
        class is ( Integrator_CS_Form )
      A_SA  =>  I % X
      G_SA  =>  I % Geometry_X
      F_SA  =>  I % CurrentSet_X
      end select !-- I
    end if

    allocate ( U % Measures )
    associate ( M  =>  U % Measures )
    call M % Initialize ( F_SA, G_SA, A_SA, Units_F = U % Units_F ( 1 ) )
    end associate !-- M

    !-- Integrator methods

    associate ( I  =>  U % Integrator )
    I % Compute_dT_Local          =>  Compute_dT_Local
    I % InitializeSeries          =>  InitializeSeries
    I % Analyze                   =>  Analyze_F_CC
    I % Set_T_CheckpointInterval  =>  Set_T_CheckpointInterval
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


  subroutine Analyze_F_CC ( I, Ignorability, T_Option )

    class ( Integrator_H_Form ), intent ( inout ) :: &
      I
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )

    call U % Analyze_F_C ( I, Ignorability, T_Option )

    associate ( M  =>  U % Measures )
    call M % Compute ( )
    end associate !-- M

    end select !-- U

  end subroutine Analyze_F_CC


  subroutine Set_T_CheckpointInterval ( I )

    class ( Integrator_H_Form ), intent ( inout ), target :: &
      I

    real ( KDR ) :: &
      Constant_G, &
      T_V, &
      T_N

    select type ( U  =>  I % System )
      class is ( Universe_F_CC_Form )
    associate &
      ( M  =>  U % Measures )

    associate &
      (     V_Max  =>  M % VelocityMax, &
          R_V_Max  =>  M % Radius_V_Max, &
        Rho_V_Max  =>  M % MassDensity_V_Max, &
        Rho_C      =>  M % MassDensity_C )

    !-- Time scales and CheckpointTimeInterval

    T_V  =  R_V_Max  /  max ( V_Max, sqrt ( tiny ( 0.0_KDR ) ) )

    if ( U % Dimensionless ) then
      Constant_G  =  1.0_KDR
    else
      Constant_G  =  CONSTANT % GRAVITATIONAL
    end if

    T_N  =  ( Constant_G * min ( Rho_V_Max, Rho_C ) ) ** ( -0.5_KDR )

    I % T_CheckpointInterval  =  min ( T_V, T_N )  /  I % nWrite

    !-- Display

    call Show ( 'Time Scales', I % IGNORABILITY )
    call Show ( T_V, I % Unit_T, 'T_Velocity', I % IGNORABILITY )
    call Show ( T_N, I % Unit_T, 'T_Density',  I % IGNORABILITY )

    !-- Cleanup

    end associate !-- V_Max, etc.
    end associate !-- M
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


  subroutine InitializeSeries ( I )

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

  end subroutine InitializeSeries


end module Universe_F_CC__Form
