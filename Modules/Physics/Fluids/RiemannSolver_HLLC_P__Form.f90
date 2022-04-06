module RiemannSolver_HLLC_P__Form

  !-- RiemannSolver_HartenLaxVanLeerContact_Perfect__Form

  use Basics
  use Mathematics
  use Fluid_P__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_SOLVER_SPEEDS_HLLC  =  1

  type, public, extends ( RiemannSolver_HLL_Form ) :: RiemannSolver_HLLC_P_Form
    integer ( KDI ) :: &
      N_SOLVER_SPEEDS_HLLC = N_SOLVER_SPEEDS_HLLC
    integer ( KDI ) :: &
      ALPHA_CENTER_U = 0
    type ( FieldSetForm ), allocatable :: &
      Metric_I, &
      CurrentSet_ICL, CurrentSet_ICR
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeCenterStates
  end type RiemannSolver_HLLC_P_Form

    private :: &
      ComputeCenterStates_P

    private :: &
      ComputeCenterSpeedKernel

    interface 

      module subroutine ComputeCenterSpeedKernel &
               ( F_D_IL, F_D_IR, F_S_IL, F_S_IR, M_IL, M_IR, D_IL, D_IR, &
                 S_IL, S_IR, AP_I, AM_I, M_UU, AC_I, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          F_D_IL, F_D_IR, &
          F_S_IL, F_S_IR, &
          M_IL, M_IR, &
          D_IL, D_IR, &
          S_IL, S_IR, &
          AP_I, &
          AM_I, &
          M_UU
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          AC_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeCenterSpeedKernel

      module subroutine ComputeCenterStatesKernel &
               ( M_IL, M_IR, D_IL, D_IR, &
                 V_1_IL, V_2_IL, V_3_IL, V_D_IL, &
                 V_1_IR, V_2_IR, V_3_IR, V_D_IR, &
                 S_1_IL, S_2_IL, S_3_IL, S_D_IL, &
                 S_1_IR, S_2_IR, S_3_IR, S_D_IR, &
                 G_IL, G_IR, P_IL, P_IR, &
                 AP_I, AM_I, AC_I, M_DD_11, M_DD_22, M_DD_33, &
                 M_ICL, M_ICR, D_ICL, D_ICR, &
                 V_1_ICL, V_2_ICL, V_3_ICL, V_D_ICL, &
                 V_1_ICR, V_2_ICR, V_3_ICR, V_D_ICR, &
                 S_1_ICL, S_2_ICL, S_3_ICL, S_D_ICL, &
                 S_1_ICR, S_2_ICR, S_3_ICR, S_D_ICR, &
                 G_ICL, G_ICR, P_ICL, P_ICR, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_IL, M_IR, &
          D_IL, D_IR, &
          V_1_IL, V_2_IL, V_3_IL, V_D_IL, &
          V_1_IR, V_2_IR, V_3_IR, V_D_IR, &
          S_1_IL, S_2_IL, S_3_IL, S_D_IL, &
          S_1_IR, S_2_IR, S_3_IR, S_D_IR, &
          G_IL, G_IR, &
          P_IL, P_IR, &
          AP_I, &
          AM_I, &
          AC_I, &
          M_DD_11, M_DD_22, M_DD_33
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          M_ICL, M_ICR, &
          D_ICL, D_ICR, &
          V_1_ICL, V_2_ICL, V_3_ICL, V_D_ICL, &
          V_1_ICR, V_2_ICR, V_3_ICR, V_D_ICR, &
          S_1_ICL, S_2_ICL, S_3_ICL, S_D_ICL, &
          S_1_ICR, S_2_ICR, S_3_ICR, S_D_ICR, &
          G_ICL, G_ICR, &
          P_ICL, P_ICR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeCenterStatesKernel

    end interface


contains


  subroutine InitializeAllocate_RS &
               ( RS, CS, FieldOption, ReconstructedSetOption, PrefixOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      ReconstructedSetOption, &
      PrefixOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      oF, &  !-- oField
      nFields
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RS % Type  ==  '' ) &
      RS % Type  =  'a RiemannSolver_HLLC_P' 
    
    associate ( nB  =>  CS % nBalanced )

    !-- Field indices

    oF  =  nB  +  RS % N_SOLVER_SPEEDS_HLL

    nFields  =  oF  +  RS % N_SOLVER_SPEEDS_HLLC
    if ( present ( nFieldsOption ) ) &
      nFields  =  nFieldsOption

    RS % ALPHA_CENTER_U   =  oF  +  1

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( oF + 1 : oF + RS % N_SOLVER_SPEEDS_HLLC ) &
      =  [ 'AlphaCenter_U' ]
          
    !-- FieldSet

    call RS % RiemannSolver_HLL_Form % Initialize &
           ( CS, &
             FieldOption = Field, &
             ReconstructedSetOption = ReconstructedSetOption, &
             PrefixOption = PrefixOption, &
             nFieldsOption = nFields )

    end associate !-- nB

    !-- Reconstructed metric

    allocate ( RS % Metric_I )
    associate ( M_I  =>  RS % Metric_I )
    call M_I % Initialize &
           ( CS % Atlas, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = 6 )
    end associate !-- M_I

    !-- Center states

    allocate &
      ( RS % CurrentSet_ICL, &
        RS % CurrentSet_ICR )
    associate &
      ( CS_ICL  =>  RS % CurrentSet_ICL, &
        CS_ICR  =>  RS % CurrentSet_ICR )
    call CS_ICL % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Field, &
             NameOption = trim ( CS % Name ) // '_ICL', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             UnitOption = CS % Unit, &
             nFieldsOption = CS % nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    call CS_ICR % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Field, &
             NameOption = trim ( CS % Name ) // '_ICR', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             UnitOption = CS % Unit, &
             nFieldsOption = CS % nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- CS_ICL, etc.

  end subroutine InitializeAllocate_RS


  subroutine Compute ( RS, DP, iC, iD, T_Option )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iDensity, &
      iMomentum
    real ( KDR ), dimension ( : ), pointer :: &
      M_UU

    call RS % RiemannSolver_HLL_Form % Compute ( DP, iC, iD, T_Option )

    select type ( CS  =>  RS % CurrentSet )
      class is ( Fluid_P_Form )
    associate &
      (  G      =>  CS % Geometry, &
        CS_ICL  =>  RS % CurrentSet_ICL, &
        CS_ICR  =>  RS % CurrentSet_ICR, &
         M_I    =>  RS % Metric_I)
    
    call G % ComputeReconstruction ( M_I, iC, iD )

    associate &
      ( RSV      =>  RS % Storage ( iC ) % Value, &
        FS_IL_V  =>  RS % FluxSet_IL % Storage ( iC ) % Value, &
        FS_IR_V  =>  RS % FluxSet_IR % Storage ( iC ) % Value, &
        CS_IL_V  =>  RS % CurrentSet_IL % Storage ( iC ) % Value, &
        CS_IR_V  =>  RS % CurrentSet_IR % Storage ( iC ) % Value )
    associate &
      ( M_DD_11  =>  M_I % Storage ( iC ) % Value ( :, 1 ), &
        M_DD_22  =>  M_I % Storage ( iC ) % Value ( :, 2 ), &
        M_DD_33  =>  M_I % Storage ( iC ) % Value ( :, 3 ), &
        M_UU_11  =>  M_I % Storage ( iC ) % Value ( :, 4 ), &
        M_UU_22  =>  M_I % Storage ( iC ) % Value ( :, 5 ), &
        M_UU_33  =>  M_I % Storage ( iC ) % Value ( :, 6 ) )

    select case ( iD )
    case ( 1 )
      M_UU  =>  M_UU_11
    case ( 2 ) 
      M_UU  =>  M_UU_22
    case ( 3 )
      M_UU  =>  M_UU_33
    end select

    call Search ( CS % iaBalanced, CS % BARYON_DENSITY_C,          iDensity )
    call Search ( CS % iaBalanced, CS % MOMENTUM_DENSITY_D ( iD ), iMomentum )

    call ComputeCenterSpeedKernel &
           ( F_D_IL = FS_IL_V ( :, iDensity ), &
             F_D_IR = FS_IR_V ( :, iDensity ), &
             F_S_IL = FS_IL_V ( :, iMomentum ), &
             F_S_IR = FS_IR_V ( :, iMomentum ), &
               M_IL = CS_IL_V ( :, CS % BARYON_MASS ), &
               M_IR = CS_IR_V ( :, CS % BARYON_MASS ), &
               D_IL = CS_IL_V ( :, CS % BARYON_DENSITY_C ), &
               D_IR = CS_IR_V ( :, CS % BARYON_DENSITY_C ), &
               S_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               S_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
              AP_I  = RSV     ( :, RS % ALPHA_PLUS_U ), &
              AM_I  = RSV     ( :, RS % ALPHA_MINUS_U ), &
               M_UU = M_UU, &
              AC_I  = RSV     ( :, RS % ALPHA_CENTER_U ), &
             UseDeviceOption = RS % DeviceMemory )

    call RS % ComputeCenterStates ( iC, iD )

    end associate !-- M_DD_11, etc.
    end associate !-- RSV, etc.
    end associate !-- G, etc.
    end select !-- CS

  end subroutine Compute


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS

    if ( allocated ( RS % CurrentSet_ICR ) ) &
      deallocate ( RS % CurrentSet_ICR )
    if ( allocated ( RS % CurrentSet_ICL ) ) &
      deallocate ( RS % CurrentSet_ICL )
    if ( allocated ( RS % Metric_I ) ) &
      deallocate ( RS % Metric_I )

  end subroutine Finalize


  subroutine ComputeCenterStates ( RS, iC, iD )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions

    call ComputeCenterStates_P ( RS, iC, iD )

  end subroutine ComputeCenterStates


  subroutine ComputeCenterStates_P ( RS, iC, iD )

    class ( RiemannSolver_HLLC_P_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions

    select type ( CS  =>  RS % CurrentSet )
      class is ( Fluid_P_Form )
    associate &
      ( CS_ICL  =>  RS % CurrentSet_ICL, &
        CS_ICR  =>  RS % CurrentSet_ICR, &
         M_I    =>  RS % Metric_I)
    associate &
      ( RSV       =>  RS % Storage ( iC ) % Value, &
        CS_ICL_V  =>  RS % CurrentSet_ICL % Storage ( iC ) % Value, &
        CS_ICR_V  =>  RS % CurrentSet_ICR % Storage ( iC ) % Value, &
        CS_IL_V   =>  RS % CurrentSet_IL % Storage ( iC ) % Value, &
        CS_IR_V   =>  RS % CurrentSet_IR % Storage ( iC ) % Value )
    associate &
      ( M_DD_11  =>  M_I % Storage ( iC ) % Value ( :, 1 ), &
        M_DD_22  =>  M_I % Storage ( iC ) % Value ( :, 2 ), &
        M_DD_33  =>  M_I % Storage ( iC ) % Value ( :, 3 ), &
        M_UU_11  =>  M_I % Storage ( iC ) % Value ( :, 4 ), &
        M_UU_22  =>  M_I % Storage ( iC ) % Value ( :, 5 ), &
        M_UU_33  =>  M_I % Storage ( iC ) % Value ( :, 6 ) )

      call ComputeCenterStatesKernel &
             ( M_IL   = CS_IL_V ( :, CS % BARYON_MASS ), &
               M_IR   = CS_IR_V ( :, CS % BARYON_MASS ), &
               D_IL   = CS_IL_V ( :, CS % BARYON_DENSITY_B ), &
               D_IR   = CS_IR_V ( :, CS % BARYON_DENSITY_B ), &
               V_1_IL = CS_IL_V ( :, CS % VELOCITY_U_1 ), &
               V_2_IL = CS_IL_V ( :, CS % VELOCITY_U_2 ), &
               V_3_IL = CS_IL_V ( :, CS % VELOCITY_U_3 ), &
               V_D_IL = CS_IL_V ( :, CS % VELOCITY_U ( iD ) ), &
               V_1_IR = CS_IR_V ( :, CS % VELOCITY_U_1 ), &
               V_2_IR = CS_IR_V ( :, CS % VELOCITY_U_2 ), &
               V_3_IR = CS_IR_V ( :, CS % VELOCITY_U_3 ), &
               V_D_IR = CS_IR_V ( :, CS % VELOCITY_U ( iD ) ), &
               S_1_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_D_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               S_1_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_D_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               G_IL   = CS_IL_V ( :, CS % ENERGY_DENSITY_B ), &
               G_IR   = CS_IR_V ( :, CS % ENERGY_DENSITY_B ), &
               P_IL   = CS_IL_V ( :, CS % PRESSURE ), &
               P_IR   = CS_IR_V ( :, CS % PRESSURE ), &
               AP_I   = RSV ( :, RS % ALPHA_PLUS_U ), &
               AM_I   = RSV ( :, RS % ALPHA_MINUS_U ), &
               AC_I   = RSV ( :, RS % ALPHA_CENTER_U ), &
               M_DD_11 = M_DD_11, &
               M_DD_22 = M_DD_22, &
               M_DD_33 = M_DD_33, &
               M_ICL   = CS_ICL_V ( :, CS % BARYON_MASS ), &
               M_ICR   = CS_ICR_V ( :, CS % BARYON_MASS ), &
               D_ICL   = CS_ICL_V ( :, CS % BARYON_DENSITY_B ), &
               D_ICR   = CS_ICR_V ( :, CS % BARYON_DENSITY_B ), &
               V_1_ICL = CS_ICL_V ( :, CS % VELOCITY_U_1 ), &
               V_2_ICL = CS_ICL_V ( :, CS % VELOCITY_U_2 ), &
               V_3_ICL = CS_ICL_V ( :, CS % VELOCITY_U_3 ), &
               V_D_ICL = CS_ICL_V ( :, CS % VELOCITY_U ( iD ) ), &
               V_1_ICR = CS_ICR_V ( :, CS % VELOCITY_U_1 ), &
               V_2_ICR = CS_ICR_V ( :, CS % VELOCITY_U_2 ), &
               V_3_ICR = CS_ICR_V ( :, CS % VELOCITY_U_3 ), &
               V_D_ICR = CS_ICR_V ( :, CS % VELOCITY_U ( iD ) ), &
               S_1_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_D_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               S_1_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_D_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D ( iD ) ), &
               G_ICL   = CS_ICL_V ( :, CS % ENERGY_DENSITY_B ), &
               G_ICR   = CS_ICR_V ( :, CS % ENERGY_DENSITY_B ), &
               P_ICL   = CS_ICL_V ( :, CS % PRESSURE ), &
               P_ICR   = CS_ICR_V ( :, CS % PRESSURE ), &
               UseDeviceOption = RS % DeviceMemory )

    end associate !-- M_DD_11, etc.
    end associate !-- RSV, etc.
    end associate !-- CS_ICL, etc.
    end select !-- CS

  end subroutine ComputeCenterStates_P


end module RiemannSolver_HLLC_P__Form
