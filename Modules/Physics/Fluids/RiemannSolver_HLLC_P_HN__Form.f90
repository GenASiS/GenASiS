module RiemannSolver_HLLC_P_HN__Form

  !-- RiemannSolver_HartenLaxVanLeerContact_Perfect_HeavyNucleus__Form

  use Basics
  use Mathematics
  use Fluid_P_HN__Form
  use RiemannSolver_HLLC_P__Form

  implicit none
  private

  type, public, extends ( RiemannSolver_HLLC_P_Form ) :: &
    RiemannSolver_HLLC_P_HN_Form
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    final :: &
      Finalize
    procedure, public, pass :: &
      ComputeCenterStates
  end type RiemannSolver_HLLC_P_HN_Form

    private :: &
      ComputeCenterStates_P_HN

    private :: &
      ComputeCenterStatesKernel

    interface 

      module subroutine ComputeCenterStatesKernel &
               ( M_IL, M_IR, D_IL, D_IR, &
                 V_1_IL, V_2_IL, V_3_IL, V_1_IR, V_2_IR, V_3_IR, &
                 S_1_IL, S_2_IL, S_3_IL, S_1_IR, S_2_IR, S_3_IR, &
                 G_IL, G_IR, P_IL, P_IR, DE_IL, DE_IR, &
                 AP_I, AM_I, AC_I, M_DD_11, M_DD_22, M_DD_33, iD, &
                 M_ICL, M_ICR, D_ICL, D_ICR, &
                 V_1_ICL, V_2_ICL, V_3_ICL, V_1_ICR, V_2_ICR, V_3_ICR, &
                 S_1_ICL, S_2_ICL, S_3_ICL, S_1_ICR, S_2_ICR, S_3_ICR, &
                 G_ICL, G_ICR, P_ICL, P_ICR, DE_ICL, DE_ICR, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          M_IL, M_IR, &
          D_IL, D_IR, &
          V_1_IL, V_2_IL, V_3_IL, V_1_IR, V_2_IR, V_3_IR, &
          S_1_IL, S_2_IL, S_3_IL, S_1_IR, S_2_IR, S_3_IR, &
          G_IL, G_IR, &
          P_IL, P_IR, &
          DE_IL, DE_IR, &
          AP_I, &
          AM_I, &
          AC_I, &
          M_DD_11, M_DD_22, M_DD_33
        integer ( KDI ), intent ( in ) :: &
          iD
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          M_ICL, M_ICR, &
          D_ICL, D_ICR, &
          V_1_ICL, V_2_ICL, V_3_ICL, V_1_ICR, V_2_ICR, V_3_ICR, &
          S_1_ICL, S_2_ICL, S_3_ICL, S_1_ICR, S_2_ICR, S_3_ICR, &
          G_ICL, G_ICR, &
          P_ICL, P_ICR, &
          DE_ICL, DE_ICR
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeCenterStatesKernel

    end interface


contains


  subroutine InitializeAllocate_RS &
               ( RS, CS, FieldOption, ReconstructedSetOption, PrefixOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLLC_P_HN_Form ), intent ( inout ) :: &
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

    if ( RS % Type  ==  '' ) &
      RS % Type  =  'a RiemannSolver_HLLC_P_HN' 
    
    call RS % RiemannSolver_HLLC_P_Form % Initialize &
           ( CS, FieldOption, ReconstructedSetOption, PrefixOption, &
             nFieldsOption )

  end subroutine InitializeAllocate_RS


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLLC_P_HN_Form ), intent ( inout ) :: &
      RS

  end subroutine Finalize


  subroutine ComputeCenterStates ( RS, iC, iD )

    class ( RiemannSolver_HLLC_P_HN_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions

    call ComputeCenterStates_P_HN ( RS, iC, iD )

  end subroutine ComputeCenterStates


  subroutine ComputeCenterStates_P_HN ( RS, iC, iD )

    class ( RiemannSolver_HLLC_P_HN_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions

    select type ( CS  =>  RS % CurrentSet )
      class is ( Fluid_P_HN_Form )
    associate &
      ( M_I  =>  RS % Metric_I )
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
               V_1_IR = CS_IR_V ( :, CS % VELOCITY_U_1 ), &
               V_2_IR = CS_IR_V ( :, CS % VELOCITY_U_2 ), &
               V_3_IR = CS_IR_V ( :, CS % VELOCITY_U_3 ), &
               S_1_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_IL = CS_IL_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_1_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_IR = CS_IR_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               G_IL   = CS_IL_V ( :, CS % ENERGY_DENSITY_B ), &
               G_IR   = CS_IR_V ( :, CS % ENERGY_DENSITY_B ), &
               P_IL   = CS_IL_V ( :, CS % PRESSURE ), &
               P_IR   = CS_IR_V ( :, CS % PRESSURE ), &
               DE_IL  = CS_IL_V ( :, CS % ELECTRON_DENSITY_B ), &
               DE_IR  = CS_IR_V ( :, CS % ELECTRON_DENSITY_B ), &
               AP_I   = RSV ( :, RS % ALPHA_PLUS_U ), &
               AM_I   = RSV ( :, RS % ALPHA_MINUS_U ), &
               AC_I   = RSV ( :, RS % ALPHA_CENTER_U ), &
               M_DD_11 = M_DD_11, &
               M_DD_22 = M_DD_22, &
               M_DD_33 = M_DD_33, &
               iD = iD, &
               M_ICL   = CS_ICL_V ( :, CS % BARYON_MASS ), &
               M_ICR   = CS_ICR_V ( :, CS % BARYON_MASS ), &
               D_ICL   = CS_ICL_V ( :, CS % BARYON_DENSITY_B ), &
               D_ICR   = CS_ICR_V ( :, CS % BARYON_DENSITY_B ), &
               V_1_ICL = CS_ICL_V ( :, CS % VELOCITY_U_1 ), &
               V_2_ICL = CS_ICL_V ( :, CS % VELOCITY_U_2 ), &
               V_3_ICL = CS_ICL_V ( :, CS % VELOCITY_U_3 ), &
               V_1_ICR = CS_ICR_V ( :, CS % VELOCITY_U_1 ), &
               V_2_ICR = CS_ICR_V ( :, CS % VELOCITY_U_2 ), &
               V_3_ICR = CS_ICR_V ( :, CS % VELOCITY_U_3 ), &
               S_1_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_ICL = CS_ICL_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               S_1_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_1 ), &
               S_2_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_2 ), &
               S_3_ICR = CS_ICR_V ( :, CS % MOMENTUM_DENSITY_D_3 ), &
               G_ICL   = CS_ICL_V ( :, CS % ENERGY_DENSITY_B ), &
               G_ICR   = CS_ICR_V ( :, CS % ENERGY_DENSITY_B ), &
               P_ICL   = CS_ICL_V ( :, CS % PRESSURE ), &
               P_ICR   = CS_ICR_V ( :, CS % PRESSURE ), &
               DE_ICL  = CS_ICL_V ( :, CS % ELECTRON_DENSITY_B ), &
               DE_ICR  = CS_ICR_V ( :, CS % ELECTRON_DENSITY_B ), &
               UseDeviceOption = RS % DeviceMemory )

    end associate !-- M_DD_11, etc.
    end associate !-- RSV, etc.
    end associate !-- M_I
    end select !-- CS

  end subroutine ComputeCenterStates_P_HN


end module RiemannSolver_HLLC_P_HN__Form
