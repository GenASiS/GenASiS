module Step_RK_NM_G_1D_C__Form

  !-- Step_RungeKutta_NeutrinoMoments_Grey_1D_Collected__Form

  use Basics
  use Mathematics
  use Gravitations
  use Fluids
  use NeutrinoMoments_G__Form
  use Interactions_NM_G__Form

  implicit none
  private

  type, public, extends ( Step_RK_CS_1D_C_CS_Form ) :: Step_RK_NM_G_1D_C_Form
    real ( KDR ), dimension ( : ), allocatable :: &
      Residual_J_Eq_E,  Residual_N_Eq_E, &
      Residual_J_Eq_EB, Residual_N_Eq_EB, &
      Residual_J_Eq_X,  Residual_N_Eq_X
  contains
    procedure, private, pass :: &
      Initialize_CS_1D_C_CS
    final :: &
      Finalize
    procedure, public, pass :: &
      SolveUpdateImplicit
  end type Step_RK_NM_G_1D_C_Form

    private :: &
      SetBalancedIndices, &
      SetStoragePointers_F, &
      SetStoragePointers_R, &
      SetFieldPointers_FS_B, &
      SetFieldPointers_F, &
      SetFieldPointers_R, &
      SetFieldPointers_I, &
      SolveKernel, &
      SolveKernelDevice

    interface

    module subroutine SolveKernel &
               ( F_V, &
                 Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                 Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
                 Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
                 Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E, &
                 Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                 Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
                 Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
                 Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB, &
                 Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                 Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
                 Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
                 Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X, &
                 J_E, H_E_1, H_E_2, H_E_3, N_E, &
                 E_E, S_E_1, S_E_2, S_E_3, D_E, &
                 J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
                 T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E, &
                 J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
                 E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
                 J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
                 T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB, &
                 J_X, H_X_1, H_X_2, H_X_3, N_X, &
                 E_X, S_X_1, S_X_2, S_X_3, D_X, &
                 J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
                 T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X, &
                 E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
                 N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
                 M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, &
                 Z_F, A_F, Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F, &
                 Error, nIterations, Omega, Residual, &
                 ProperCell, &
                 ApplyImplicit_F, &
                 EOS, &
                 E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
                 E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
                 E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
                 E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
                 M_DD_11, M_DD_22, M_DD_33, &
                 M_UU_11, M_UU_22, M_UU_33, &
                 T_L_N, T_L_T, T_Ye, &
                 AA, Tol, dT, Rho_DB, &
                 M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
                 E_Shift, &
                 ia_F_I, ia_F_O, ia_E, &
                 mRI, mII, iC, iSolve, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
                 KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
                 Res_J_Eq_E,  Res_N_Eq_E, &
                 Res_J_Eq_EB, Res_N_Eq_EB, &
                 Res_J_Eq_X,  Res_N_Eq_X )
      implicit none
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        F_V
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_E,  Xi_H_E,  Xi_N_E,  Chi_J_E,  Chi_H_E,  Chi_N_E, &
        Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
        Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
        Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
        Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
        Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
        Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_X,  Xi_H_X,  Xi_N_X,  Chi_J_X,  Chi_H_X,  Chi_N_X, &
        Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
        Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
        Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_E, H_E_1, H_E_2, H_E_3, N_E, &
        E_E, S_E_1, S_E_2, S_E_3, D_E, &
        J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
        T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
        E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
        J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
        T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_X, H_X_1, H_X_2, H_X_3, N_X, &
        E_X, S_X_1, S_X_2, S_X_3, D_X, &
        J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
        T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
        N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
        M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
        Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Error, nIterations, Omega, Residual
      logical ( KDL ), dimension ( : ), intent ( in ) :: &
        ProperCell
      logical ( KDL ), intent ( in ) :: &
        ApplyImplicit_F
      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
        EOS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
        E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
        E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
        E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M_DD_11, M_DD_22, M_DD_33, &
        M_UU_11, M_UU_22, M_UU_33, &
        T_L_N, T_L_T, T_Ye
      real ( KDR ), intent ( in ) :: &
        AA, Tol, dT, &
        Rho_DB, &
        M_Ref, N_Min, E_Min, &
        T_Min, Y_Min, Y_Safe, &
        E_Shift
      integer ( KDI ), dimension ( : ), intent ( in ) :: &
        ia_F_I, ia_F_O, ia_E
      integer ( KDI ), intent ( in ) :: &
        mRI, &
        mII, &
        iC, &
        iSolve
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
        KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
        KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
        KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        Res_J_Eq_E,  Res_N_Eq_E, &
        Res_J_Eq_EB, Res_N_Eq_EB, &
        Res_J_Eq_X,  Res_N_Eq_X
    end subroutine SolveKernel

    
    module subroutine SolveKernelDevice &
               ( F_V, &
                 Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                 Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
                 Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
                 Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E, &
                 Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                 Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
                 Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
                 Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB, &
                 Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                 Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
                 Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
                 Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X, &
                 J_E, H_E_1, H_E_2, H_E_3, N_E, &
                 E_E, S_E_1, S_E_2, S_E_3, D_E, &
                 J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
                 T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E, &
                 J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
                 E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
                 J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
                 T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB, &
                 J_X, H_X_1, H_X_2, H_X_3, N_X, &
                 E_X, S_X_1, S_X_2, S_X_3, D_X, &
                 J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
                 T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X, &
                 E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
                 N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
                 M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, &
                 Z_F, A_F, Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F, &
                 Error, nIterations, Omega, Residual, &
                 ProperCell, &
                 ApplyImplicit_F, &
                 EOS, &
                 E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
                 E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
                 E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
                 E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
                 M_DD_11, M_DD_22, M_DD_33, &
                 M_UU_11, M_UU_22, M_UU_33, &
                 T_L_N, T_L_T, T_Ye, &
                 AA, Tol, dT, Rho_DB, &
                 M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, &
                 E_Shift, &
                 ia_F_I, ia_F_O, ia_E, &
                 mRI, mII, iC, iSolve, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
                 KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
                 Res_J_Eq_E,  Res_N_Eq_E, &
                 Res_J_Eq_EB, Res_N_Eq_EB, &
                 Res_J_Eq_X,  Res_N_Eq_X )
      implicit none
      real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
        F_V
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_E,  Xi_H_E,  Xi_N_E,  Chi_J_E,  Chi_H_E,  Chi_N_E, &
        Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
        Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
        Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
        Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
        Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
        Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Xi_J_X,  Xi_H_X,  Xi_N_X,  Chi_J_X,  Chi_H_X,  Chi_N_X, &
        Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
        Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
        Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_E, H_E_1, H_E_2, H_E_3, N_E, &
        E_E, S_E_1, S_E_2, S_E_3, D_E, &
        J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
        T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
        E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
        J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
        T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        J_X, H_X_1, H_X_2, H_X_3, N_X, &
        E_X, S_X_1, S_X_2, S_X_3, D_X, &
        J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
        T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
        N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
        M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
        Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F
      real ( KDR ), dimension ( : ), intent ( inout ) :: &
        Error, nIterations, Omega, Residual
      logical ( KDL ), dimension ( : ), intent ( in ) :: &
        ProperCell
      logical ( KDL ), intent ( in ) :: &
        ApplyImplicit_F
      real ( KDR ), dimension ( :, :, :, : ), intent ( in ) :: &
        EOS
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
        E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
        E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
        E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0
      real ( KDR ), dimension ( : ), intent ( in ) :: &
        M_DD_11, M_DD_22, M_DD_33, &
        M_UU_11, M_UU_22, M_UU_33, &
        T_L_N, T_L_T, T_Ye
      real ( KDR ), intent ( in ) :: &
        AA, Tol, dT, &
        Rho_DB, &
        M_Ref, N_Min, E_Min, &
        T_Min, Y_Min, Y_Safe, &
        E_Shift
      integer ( KDI ), dimension ( : ), intent ( in ) :: &
        ia_F_I, ia_F_O, ia_E
      integer ( KDI ), intent ( in ) :: &
        mRI, &
        mII, &
        iC, &
        iSolve
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
        KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
        KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
        KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D
      real ( KDR ), dimension ( : ), intent ( out ) :: &
        Res_J_Eq_E,  Res_N_Eq_E, &
        Res_J_Eq_EB, Res_N_Eq_EB, &
        Res_J_Eq_X,  Res_N_Eq_X
    end subroutine SolveKernelDevice

  end interface


contains


  subroutine Initialize_CS_1D_C_CS &
               ( S, CS_1D, CS, NameOption, ImplicitExplicitOption, &
                 ComputeExplicit_CS_Option, ComputeExplicit_CS_1D_Option, &
                 ComputeImplicit_CS_Option, ComputeImplicit_CS_1D_Option, &
                 nStagesOption )

    class ( Step_RK_NM_G_1D_C_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), dimension ( : ), intent ( in ), target :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption, &
      ComputeExplicit_CS_Option, ComputeExplicit_CS_1D_Option, &
      ComputeImplicit_CS_Option, ComputeImplicit_CS_1D_Option
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_NM_G_1D_C'

    call S % Step_RK_CS_1D_C_CS_Form % Initialize &
           ( CS_1D, CS, NameOption, ImplicitExplicitOption, &
             ComputeExplicit_CS_Option, ComputeExplicit_CS_1D_Option, &
             ComputeImplicit_CS_Option, ComputeImplicit_CS_1D_Option, &
             nStagesOption )

    associate ( mII  =>  S % MaxImplicitIterations )

    allocate ( S % Residual_J_Eq_E ( mII ) )
    allocate ( S % Residual_N_Eq_E ( mII ) )

    allocate ( S % Residual_J_Eq_EB ( mII ) )
    allocate ( S % Residual_N_Eq_EB ( mII ) )

    allocate ( S % Residual_J_Eq_X ( mII ) )
    allocate ( S % Residual_N_Eq_X ( mII ) )

    end associate !-- nNM, etc.

  end subroutine Initialize_CS_1D_C_CS


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_NM_G_1D_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Residual_N_Eq_X ) ) &
      deallocate ( S % Residual_N_Eq_X )
    if ( allocated ( S % Residual_J_Eq_X ) ) &
      deallocate ( S % Residual_J_Eq_X )

    if ( allocated ( S % Residual_N_Eq_EB ) ) &
      deallocate ( S % Residual_N_Eq_EB )
    if ( allocated ( S % Residual_J_Eq_EB ) ) &
      deallocate ( S % Residual_J_Eq_EB )

    if ( allocated ( S % Residual_N_Eq_E ) ) &
      deallocate ( S % Residual_N_Eq_E )
    if ( allocated ( S % Residual_J_Eq_E ) ) &
      deallocate ( S % Residual_J_Eq_E )

  end subroutine Finalize


  subroutine SolveUpdateImplicit  ( S, ApplyImplicit_CS, T, dT, iS )

    class ( Step_RK_NM_G_1D_C_Form ), intent ( inout ), target :: &
      S
    logical ( KDL ), intent ( in ) :: &
      ApplyImplicit_CS
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iR, &  !-- iRadiation
      iC     !-- iChart
    integer ( KDI ) :: &
      iEnergy_R, iNumber_R, &
      iEnergy_F, iNumber_F
    integer ( KDI ), dimension ( 3 ) :: &
      iMomentum_R, iMomentum_F
    !-- Field pointers
    real ( KDR ), dimension ( : ), pointer :: &
      Error, nIterations, Omega, Residual
    real ( KDR ), dimension ( : ), pointer :: &
      KK_F_E,   KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
      KK_E_E,   KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
      KK_EB_E,  KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
      KK_X_E,   KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D 
    real ( KDR ), dimension ( : ), pointer :: &
      E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
      E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
      E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
      E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0
    real ( KDR ), dimension ( : ), pointer :: &
      M_DD_11, M_DD_22, M_DD_33, &
      M_UU_11, M_UU_22, M_UU_33
    real ( KDR ), dimension ( : ), pointer :: &
      E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
      M_F, N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
      T_F, P_F, SB_F, SS_F, X_AA_F, X_p_F, X_n_F, X_A_F, &
      Z_F, A_F, Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F
    real ( KDR ), dimension ( : ), pointer :: &
      J_E, H_E_1, H_E_2, H_E_3, N_E, &
      E_E, S_E_1, S_E_2, S_E_3, D_E, &
      J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
      T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E
    real ( KDR ), dimension ( : ), pointer :: &
      J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
      E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
      J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
      T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB
    real ( KDR ), dimension ( : ), pointer :: &
      J_X, H_X_1, H_X_2, H_X_3, N_X, &
      E_X, S_X_1, S_X_2, S_X_3, D_X, &
      J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
      T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X
    real ( KDR ), dimension ( : ), pointer :: &
      Xi_J_E,  Xi_H_E,  Xi_N_E,  &
      Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
      Chi_J_E,  Chi_H_E,  Chi_N_E, &
      Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
      Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E, &
      Xi_J_EB,  Xi_H_EB,  Xi_N_EB,  &
      Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
      Chi_J_EB,  Chi_H_EB,  Chi_N_EB, &
      Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
      Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB, &
      Xi_J_X,  Xi_H_X,  Xi_N_X,  &
      Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
      Chi_J_X,  Chi_H_X,  Chi_N_X, &
      Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
      Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X
    !-- Storage % Value pointers
    real ( KDR ), dimension ( :, : ), pointer :: &
      ID_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      KK_F_V, KK_E_V, KK_EB_V, KK_X_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      Y_I_F_V, Y_I_E_V, Y_I_EB_V, Y_I_X_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      G_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      F_V
    real ( KDR ), dimension ( :, : ), pointer :: &
      R_E_V, R_EB_V, R_X_V, &
      I_E_V, I_EB_V, I_X_V
    !-- FieldSet pointers
    class ( ImplicitDiagnosticsForm ), pointer :: &
      ID
    class ( FieldSet_BM_Form ), pointer :: &
      KK_F, KK_E, KK_EB, KK_X
    class ( FieldSet_BM_Form ), pointer :: &
      Y_I_F, Y_I_E, Y_I_EB, Y_I_X
    class ( Gravitation_N_SG_Form ), pointer :: &
      G_N
    class ( Fluid_P_HN_Form ), pointer :: &
      F_HN
    class ( NeutrinoMoments_G_Form ), pointer :: &
      R_E, R_EB, R_X
    class ( Interactions_NM_G_Form ), pointer :: &
      I_E, I_EB, I_X

integer ( KDI ) :: &
  iV

    call Show ( 'SolveUpdateImplicit', CONSOLE % INFO_5 )
    call Show ( S % Name, 'Step', CONSOLE % INFO_5 )

    associate &
      (  S_R  =>  S % Step_CS_1D ( : ), &
         S_F  =>  S % Step_CS, &
         Res_J_Eq_E   =>  S % Residual_J_Eq_E, &
         Res_N_Eq_E   =>  S % Residual_N_Eq_E, &
         Res_J_Eq_EB  =>  S % Residual_J_Eq_EB, &
         Res_N_Eq_EB  =>  S % Residual_N_Eq_EB, &
         Res_J_Eq_X   =>  S % Residual_J_Eq_X, &
         Res_N_Eq_X   =>  S % Residual_N_Eq_X, &
         AA   =>  S % AA ( iS ) % Value ( iS ), &
         Tol  =>  S % ImplicitTolerance, &
        mII   =>  S % MaxImplicitIterations, &
        mRI   =>  S % MaxRelaxationIterations, & 
        nR    =>  S % nCurrentSets_1D )

    !-- FieldSet pointers

    ID  =>  S % ImplicitDiagnostics ( iS )

    select type ( F  =>  S_F % CurrentSet )
    class is ( Fluid_P_HN_Form )
      F_HN   =>  F 
      Y_I_F  =>  S_F % Intermediate
       KK_F  =>  S_F % SlopeStageImplicit ( iS ) % Element
    end select !-- F

    select type ( G  =>  F_HN % Geometry )
    class is ( Gravitation_N_SG_Form )
      G_N  =>  G
    end select !-- G

    do iR  =  1,  nR
      select type ( R  =>  S_R ( iR ) % CurrentSet )
        class is ( NeutrinoMoments_G_Form )
      select type ( I  =>  R % Interactions )
        class is ( Interactions_NM_G_Form )
      select case ( trim ( R % RadiationType ) )
      case ( 'NEUTRINOS_E' )
          I_E  =>  I
          R_E  =>  R
        Y_I_E  =>  S_R ( iR ) % Intermediate
         KK_E  =>  S_R ( iR ) % SlopeStageImplicit ( iS ) % Element
      case ( 'NEUTRINOS_EB' )
          I_EB  =>  I
          R_EB  =>  R
        Y_I_EB  =>  S_R ( iR ) % Intermediate
         KK_EB  =>  S_R ( iR ) % SlopeStageImplicit ( iS ) % Element
      case ( 'NEUTRINOS_HL' )
          I_X  =>  I
          R_X  =>  R
        Y_I_X  =>  S_R ( iR ) % Intermediate
         KK_X  =>  S_R ( iR ) % SlopeStageImplicit ( iS ) % Element
      case default
        call Show ( 'RadiationType not recognized', CONSOLE % ERROR )
        call Show ( R % RadiationType, 'RadiationType', &
                    CONSOLE % ERROR )
        call Show ( 'Step_RK_NM_G_1D_C__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- RadiationType
      end select !-- I
      end select !-- R
    end do !-- iR

    call SetBalancedIndices &
           ( R_E, F_HN, iMomentum_R, iMomentum_F, iEnergy_R, iEnergy_F, &
             iNumber_R, iNumber_F )

    do iC  =  1,  S % Atlas % nCharts
      select type ( C  =>  S % Atlas % Chart ( iC ) % Element )
        class is ( Chart_GS_Form )

      !-- Storage % Value pointers

      ID_V  =>   ID % Storage ( iC ) % Value

       G_V  =>  G_N % Storage ( iC ) % Value

      call SetStoragePointers_F &
             ( F_HN, Y_I_F,   KK_F,  iC, &
               F_V,  Y_I_F_V, KK_F_V )
      call SetStoragePointers_R &
             ( I_E,   R_E,   Y_I_E,   KK_E,  iC, &
               I_E_V, R_E_V, Y_I_E_V, KK_E_V )
      call SetStoragePointers_R &
             ( I_EB,   R_EB,   Y_I_EB,   KK_EB,  iC, &
               I_EB_V, R_EB_V, Y_I_EB_V, KK_EB_V )
      call SetStoragePointers_R &
             ( I_X,   R_X,   Y_I_X,   KK_X,  iC, &
               I_X_V, R_X_V, Y_I_X_V, KK_X_V )

      !-- Field pointers

            Error  =>  ID_V ( :, ID % ERROR )
      nIterations  =>  ID_V ( :, ID % N_ITERATIONS )
            Omega  =>  ID_V ( :, ID % RELAXATION )
         Residual  =>  ID_V ( :, ID % RESIDUAL )

      M_DD_11  =>  G_V ( :, G_N % METRIC_F_DD_11 )
      M_DD_22  =>  G_V ( :, G_N % METRIC_F_DD_22 )
      M_DD_33  =>  G_V ( :, G_N % METRIC_F_DD_33 )
      M_UU_11  =>  G_V ( :, G_N % METRIC_F_UU_11 )
      M_UU_22  =>  G_V ( :, G_N % METRIC_F_UU_22 )
      M_UU_33  =>  G_V ( :, G_N % METRIC_F_UU_33 )

      call SetFieldPointers_FS_B &
             ( KK_F_V, iMomentum_F, iEnergy_F, iNumber_F, &
               KK_F_S_1, KK_F_S_2, KK_F_S_3, KK_F_E, KK_F_D )
      call SetFieldPointers_FS_B &
             ( KK_E_V, iMomentum_R, iEnergy_R, iNumber_R, &
               KK_E_S_1, KK_E_S_2, KK_E_S_3, KK_E_E, KK_E_D )
      call SetFieldPointers_FS_B &
             ( KK_EB_V, iMomentum_R, iEnergy_R, iNumber_R, &
               KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_E, KK_EB_D )
      call SetFieldPointers_FS_B &
             ( KK_X_V, iMomentum_R, iEnergy_R, iNumber_R, &
               KK_X_S_1, KK_X_S_2, KK_X_S_3, KK_X_E, KK_X_D )

      call SetFieldPointers_FS_B &
             ( Y_I_F_V, iMomentum_F, iEnergy_F, iNumber_F, &
               S_F_1_0, S_F_2_0, S_F_3_0, E_F_0, D_F_0 )
      call SetFieldPointers_FS_B &
             ( Y_I_E_V, iMomentum_R, iEnergy_R, iNumber_R, &
               S_E_1_0, S_E_2_0, S_E_3_0, E_E_0, D_E_0 )
      call SetFieldPointers_FS_B &
             ( Y_I_EB_V, iMomentum_R, iEnergy_R, iNumber_R, &
               S_EB_1_0, S_EB_2_0, S_EB_3_0, E_EB_0, D_EB_0 )
      call SetFieldPointers_FS_B &
             ( Y_I_X_V, iMomentum_R, iEnergy_R, iNumber_R, &
               S_X_1_0, S_X_2_0, S_X_3_0, E_X_0, D_X_0 )

      call SetFieldPointers_F &
             ( F_HN, F_V, E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
               M_F, N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
               T_F, P_F, SB_F, SS_F, X_AA_F, X_p_F, X_n_F, X_A_F, Z_F, A_F, &
               Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F )

      call SetFieldPointers_R &
             ( R_E, R_E_V, &
               J_E, H_E_1, H_E_2, H_E_3, N_E, &
               E_E, S_E_1, S_E_2, S_E_3, D_E, &
               J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
               T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E )
      call SetFieldPointers_R &
             ( R_EB, R_EB_V, &
               J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
               E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
               J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
               T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB )
      call SetFieldPointers_R &
             ( R_X, R_X_V, &
               J_X, H_X_1, H_X_2, H_X_3, N_X, &
               E_X, S_X_1, S_X_2, S_X_3, D_X, &
               J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
               T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X )

      call SetFieldPointers_I &
             ( I_E, I_E_V, &
               Xi_J_E, Xi_H_E, Xi_N_E, &
               Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
               Chi_J_E, Chi_H_E, Chi_N_E, &
               Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
               Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E )
      call SetFieldPointers_I &
             ( I_EB, I_EB_V, &
               Xi_J_EB, Xi_H_EB, Xi_N_EB, &
               Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
               Chi_J_EB, Chi_H_EB, Chi_N_EB, &
               Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
               Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB )
      call SetFieldPointers_I &
             ( I_X, I_X_V, &
               Xi_J_X, Xi_H_X, Xi_N_X, &
               Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
               Chi_J_X, Chi_H_X, Chi_N_X, &
               Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
               Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X )

      associate &
        ( Rho_DB  => I_E % DensityDetailedBalance, &
          EOS     => F_HN % EOS % Table, &
          M_Ref   => F_HN % BaryonMass, &   
          N_Min   => F_HN % BaryonDensityMin, &
          E_Min   => F_HN % EnergyDensityMin, &
          T_Min   => F_HN % TemperatureMin, &
          Y_Min   => F_HN % ElectronFractionMin, &
          Y_Safe  => F_HN % ElectronFractionSafe, &
          T_L_N   => F_HN % EOS % LogDensity, &
          T_L_T   => F_HN % EOS % LogTemperature, &
          T_Ye    => F_HN % EOS % ElectronFraction, &
          E_Shift => F_HN % EOS % EnergyShift, &
          ia_F_I  => [ F_HN % BARYON_DENSITY_C, &
                       F_HN % TEMPERATURE, F_HN % ELECTRON_FRACTION ], &
          ia_F_O  => F_HN % EOS % iaFluidOutput, &
          ia_E    => F_HN % EOS % iaSelected, &
          iSolve  => F_HN % ENERGY_DENSITY_C )

      if ( F_HN % DeviceMemory ) then
        call SolveKernelDevice &
               ( F_V, &
                 Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                 Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
                 Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
                 Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E, &
                 Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                 Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
                 Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
                 Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB, &
                 Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                 Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
                 Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
                 Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X, &
                 J_E, H_E_1, H_E_2, H_E_3, N_E, &
                 E_E, S_E_1, S_E_2, S_E_3, D_E, &
                 J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
                 T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E, &
                 J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
                 E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
                 J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
                 T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB, &
                 J_X, H_X_1, H_X_2, H_X_3, N_X, &
                 E_X, S_X_1, S_X_2, S_X_3, D_X, &
                 J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
                 T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X, &
                 E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
                 N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
                 M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, &
                 Z_F, A_F, Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F, &
                 Error, nIterations, Omega, Residual, &
                 C % ProperCell, &
                 ApplyImplicit_CS, &
                 EOS, &
                 E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
                 E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
                 E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
                 E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
                 M_DD_11, M_DD_22, M_DD_33, &
                 M_UU_11, M_UU_22, M_UU_33, &
                 T_L_N, T_L_T, T_Ye, & 
                 AA, Tol, dT, Rho_DB, &
                 M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, E_Shift, &
                 ia_F_I, ia_F_O, ia_E, &
                 mRI, mII, iC, iSolve, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
                 KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
                 Res_J_Eq_E,  Res_N_Eq_E, &
                 Res_J_Eq_EB, Res_N_Eq_EB, &
                 Res_J_Eq_X,  Res_N_Eq_X )
      else 
        call SolveKernel &
               ( F_V, &
                 Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                 Xi_J_EA_N_E, Xi_J_EA_A_E, Xi_J_P_EP_E, Xi_J_S_EP_E, &
                 Chi_J_EA_N_E, Chi_J_EA_A_E, Chi_J_P_EP_E, Chi_J_S_EP_E, &
                 Chi_H_S_N_E, Chi_H_S_A_E, Chi_H_S_EP_E, &
                 Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                 Xi_J_EA_N_EB, Xi_J_EA_A_EB, Xi_J_P_EP_EB, Xi_J_S_EP_EB, &
                 Chi_J_EA_N_EB, Chi_J_EA_A_EB, Chi_J_P_EP_EB, Chi_J_S_EP_EB, &
                 Chi_H_S_N_EB, Chi_H_S_A_EB, Chi_H_S_EP_EB, &
                 Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                 Xi_J_EA_N_X, Xi_J_EA_A_X, Xi_J_P_EP_X, Xi_J_S_EP_X, &
                 Chi_J_EA_N_X, Chi_J_EA_A_X, Chi_J_P_EP_X, Chi_J_S_EP_X, &
                 Chi_H_S_N_X, Chi_H_S_A_X, Chi_H_S_EP_X, &
                 J_E, H_E_1, H_E_2, H_E_3, N_E, &
                 E_E, S_E_1, S_E_2, S_E_3, D_E, &
                 J_Eq_E, N_Eq_E, J_Rd_E, N_Rd_E, &
                 T_Nu_E, Eta_Nu_E, E_Ave_E, F_Ave_E, &
                 J_EB, H_EB_1, H_EB_2, H_EB_3, N_EB, &
                 E_EB, S_EB_1, S_EB_2, S_EB_3, D_EB, &
                 J_Eq_EB, N_Eq_EB, J_Rd_EB, N_Rd_EB, &
                 T_Nu_EB, Eta_Nu_EB, E_Ave_EB, F_Ave_EB, &
                 J_X, H_X_1, H_X_2, H_X_3, N_X, &
                 E_X, S_X_1, S_X_2, S_X_3, D_X, &
                 J_Eq_X, N_Eq_X, J_Rd_X, N_Rd_X, &
                 T_Nu_X, Eta_Nu_X, E_Ave_X, F_Ave_X, &
                 E_F, S_F_1, S_F_2, S_F_3, D_F, DB_F, &
                 N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F, &
                 M_F, T_F, P_F, SB_F, SS_F, X_AA_F, X_n_F, X_p_F, X_A_F, &
                 Z_F, A_F, Mu_e_F, Mu_n_p_F, Mu_p_F, Mu_n_F, G_F, &
                 Error, nIterations, Omega, Residual, &
                 C % ProperCell, &
                 ApplyImplicit_CS, &
                 EOS, &
                 E_E_0,  S_E_1_0,  S_E_2_0,  S_E_3_0,  D_E_0, &
                 E_EB_0, S_EB_1_0, S_EB_2_0, S_EB_3_0, D_EB_0, &
                 E_X_0,  S_X_1_0,  S_X_2_0,  S_X_3_0,  D_X_0, &
                 E_F_0,  S_F_1_0,  S_F_2_0,  S_F_3_0,  D_F_0, &
                 M_DD_11, M_DD_22, M_DD_33, &
                 M_UU_11, M_UU_22, M_UU_33, &
                 T_L_N, T_L_T, T_Ye, & 
                 AA, Tol, dT, Rho_DB, &
                 M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, E_Shift, &
                 ia_F_I, ia_F_O, ia_E, &
                 mRI, mII, iC, iSolve, &
                 KK_E_E,  KK_E_S_1,  KK_E_S_2,  KK_E_S_3,  KK_E_D, & 
                 KK_EB_E, KK_EB_S_1, KK_EB_S_2, KK_EB_S_3, KK_EB_D, & 
                 KK_X_E,  KK_X_S_1,  KK_X_S_2,  KK_X_S_3,  KK_X_D, & 
                 KK_F_E,  KK_F_S_1,  KK_F_S_2,  KK_F_S_3,  KK_F_D, &
                 Res_J_Eq_E,  Res_N_Eq_E, &
                 Res_J_Eq_EB, Res_N_Eq_EB, &
                 Res_J_Eq_X,  Res_N_Eq_X )
      end if
      
      end associate   !-- Rho_DB, ...

! call Show ( KK_E_V, '>>> KK_E_V' )
! call Show ( KK_EB_V, '>>> KK_EB_V' )
! call Show ( KK_X_V, '>>> KK_X_V' )
! call Show ( KK_F_V, '>>> KK_F_V' )

! iV  =  40
! associate &
!   (  Y_E    =>  S_R ( 1 ) % Solution, &
!      K_E_1  =>  S_R ( 1 ) % SlopeStageExplicit ( 1 ) % Element, &
!      K_E_2  =>  S_R ( 1 ) % SlopeStageExplicit ( 2 ) % Element, &
!     KK_E_2  =>  S_R ( 1 ) % SlopeStageImplicit ( 2 ) % Element, &
!     KK_E_3  =>  S_R ( 1 ) % SlopeStageImplicit ( 3 ) % Element )
! associate &
!   (  Y_E_V    =>   Y_E   % Storage ( iC ) % Value, &
!      K_E_1_V  =>   K_E_1 % Storage ( iC ) % Value, &
!      K_E_2_V  =>   K_E_2 % Storage ( iC ) % Value, &
!     KK_E_2_V  =>  KK_E_2 % Storage ( iC ) % Value, &
!     KK_E_3_V  =>  KK_E_3 % Storage ( iC ) % Value )
! associate &
!   (  R    =>  G_V ( :, G_N % CENTER_U_1 ), &
!      J_N  =>  Y_E_V ( :, iEnergy_R ), &
!      H_N  =>  Y_E_V ( :, iMomentum_R ( 1 ) ), &
!      K_E_1_E  =>   K_E_1_V ( :, iEnergy_R ), &
!      K_E_1_H  =>   K_E_1_V ( :, iMomentum_R ( 1 ) ), &
!      K_E_2_E  =>   K_E_2_V ( :, iEnergy_R ), &
!      K_E_2_H  =>   K_E_2_V ( :, iMomentum_R ( 1 ) ), &
!     KK_E_2_E  =>  KK_E_2_V ( :, iEnergy_R ), &
!     KK_E_2_H  =>  KK_E_2_V ( :, iMomentum_R ( 1 ) ), &
!     KK_E_3_E  =>  KK_E_3_V ( :, iEnergy_R ), &
!     KK_E_3_H  =>  KK_E_3_V ( :, iMomentum_R ( 1 ) ) )
! if ( iS == 2 ) then
!   call Show ( '>>> Stage' )
!   call Show ( iS, '>>> iS' )
!   call Show ( dT, '>>> dT' )
!   call Show ( '>>> Cell' )
!   call Show ( iV, '>>> iV' )
!   call Show ( R ( iV ), UNIT % KILOMETER, '>>> R' )
!   call Show ( '>>> Before implicit solve' )
!     call Show ( '>>> Energy' )
!     call Show ( J_N ( iV ), '>>> J_(1)' )
!     call Show ( dT * K_E_1_E ( iV ), '>>> dT * K_E_E_(1)' )
!     call Show ( E_E_0 ( iV ), '>>> J_(1+)' )
!     call Show ( J_N ( iV ) + dT * K_E_1_E ( iV ), '>>> J_(1+) check' )
!     call Show ( K_E_1_E ( iV ), '>>> K_E_E_(1)' )
!     call Show ( - 2. * H_N ( iV ) / R ( iV ), '>>> - 2H/R' )
!     call Show ( - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
!                   / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
!                 '>>> - dH/dR' )
!     call Show ( - 2. * H_N ( iV ) / R ( iV ) &
!                 - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
!                   / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
!                 '>>> - ( 2H/R + dH/dR )' )
!     call Show ( '>>> Momentum' )
!     call Show ( H_N ( iV ), '>>> H_(1)' )
!     call Show ( dT * K_E_1_H ( iV ), '>>> dT * K_E_H_(1)' )
!     call Show ( S_E_1_0 ( iV ), '>>> H_(1+)' )
!     call Show ( H_N ( iV ) + dT * K_E_1_H ( iV ), '>>> H_(1+) check' )
!     call Show ( K_E_1_H ( iV ), '>>> K_E_H_(1)' )
!     call Show ( - ( J_N ( iV + 1 ) - J_N ( iV - 1 ) ) &
!                   / ( 3. * ( R ( iV + 1 ) - R ( iV - 1 ) ) ), &
!                 '>>> -(1/3) dJ/dR' )
!   call Show ( '>>> After implicit solve' )
!     call Show ( '>>> Energy' )
!     call Show ( dT * KK_E_2_E ( iV ), '>>> dT * KK_E_E_(2)' )
!     call Show ( E_E_0 ( iV ) +  dT * KK_E_2_E ( iV ), '>>> J_(2)' )
!     call Show ( '>>> Momentum' )
!     call Show ( dT * KK_E_2_H ( iV ), '>>> dT * KK_E_H_(2)' )
!     call Show ( S_E_1_0 ( iV ) + dT * KK_E_2_H ( iV ), '>>> H_(2)' )
!     call Show ( K_E_1_H ( iV ) / Chi_H_E ( iV ), 'K_E_1_H / Chi_H_E' )
!     call Show ( Chi_H_E ( iV ), '>>> Chi_H_E' )
! else if ( iS == 3 ) then
!   call Show ( '>>> Stage' )
!   call Show ( iS, '>>> iS' )
!   call Show ( dT, '>>> dT' )
!   call Show ( '>>> Cell' )
!   call Show ( iV, '>>> iV' )
!   call Show ( R ( iV ), UNIT % KILOMETER, '>>> R' )
!   call Show ( '>>> Before implicit solve' )
!     call Show ( '>>> Energy' )
!     call Show ( J_N ( iV ), '>>> J_(1)' )
!     call Show ( 0.5 * dT * K_E_1_E ( iV ), '>>> (1/2) dT * K_E_E_(1)' )
!     call Show ( 0.5 * dT * K_E_2_E ( iV ), '>>> (1/2) dT * K_E_E_(2)' )
!     call Show ( 0.5 * dT * KK_E_2_E ( iV ), '>>> (1/2) dT * KK_E_E_(2)' )
!     call Show ( E_E_0 ( iV ), '>>> J_(2+)' )
!     call Show ( J_N ( iV ) + 0.5 * dT * K_E_1_E ( iV ) &
!                  + 0.5 * dT * K_E_2_E ( iV ) + 0.5 * dT * KK_E_2_E ( iV ), &
!                 '>>> J_(2+) check' )
!   !   call Show ( K_E_1_E ( iV ), '>>> K_E_E_(1)' )
!   !   call Show ( - 2. * H_N ( iV ) / R ( iV ), '>>> - 2H/R' )
!   !   call Show ( - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
!   !                 / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
!   !               '>>> - dH/dR' )
!   !   call Show ( - 2. * H_N ( iV ) / R ( iV ) &
!   !               - ( H_N ( iV + 1 ) - H_N ( iV - 1 ) ) &
!   !                 / ( R ( iV + 1 ) - R ( iV - 1 ) ), &
!   !               '>>> - ( 2H/R + dH/dR )' )
!     call Show ( '>>> Momentum' )
!     call Show ( H_N ( iV ), '>>> H_(1)' )
!     call Show ( 0.5 * dT * K_E_1_H ( iV ), '>>> (1/2) dT * K_E_H_(1)' )
!     call Show ( 0.5 * dT * K_E_2_H ( iV ), '>>> (1/2) dT * K_E_H_(2)' )
!     call Show ( 0.5 * dT * KK_E_2_H ( iV ), '>>> (1/2) dT * KK_E_H_(2)' )
!     call Show ( S_E_1_0 ( iV ), '>>> H_(2+)' )
!     call Show ( H_N ( iV ) + 0.5 * dT * K_E_1_H ( iV ) &
!                  + 0.5 * dT * K_E_2_H ( iV ) + 0.5 * dT * KK_E_2_H ( iV ), &
!                 '>>> H_(2+) check' )
!   !   call Show ( H_N ( iV ) + dT * K_E_1_H ( iV ), '>>> H_(1+) check' )
!   !   call Show ( K_E_1_H ( iV ), '>>> K_E_H_(1)' )
!   !   call Show ( - ( J_N ( iV + 1 ) - J_N ( iV - 1 ) ) &
!   !                 / ( 3. * ( R ( iV + 1 ) - R ( iV - 1 ) ) ), &
!   !               '>>> -(1/3) dJ/dR' )
!   call Show ( '>>> After implicit solve' )
!     call Show ( '>>> Energy' )
!     call Show ( 0.5 * dT * KK_E_3_E ( iV ), '>>> (1/2) dT * KK_E_E_(3)' )
!     call Show ( E_E_0 ( iV ) + 0.5 * dT * KK_E_3_E ( iV ), '>>> J_(3)' )
!     call Show ( '>>> Momentum' )
!     call Show ( 0.5 * dT * KK_E_3_H ( iV ), '>>> (1/2) dT * KK_E_H_(3)' )
!     call Show ( S_E_1_0 ( iV ) + 0.5 * dT * KK_E_3_H ( iV ), '>>> H_(3)' )
!   !   call Show ( K_E_1_H ( iV ) / Chi_H_E ( iV ), 'K_E_1_H / Chi_H_E' )
!   !   call Show ( Chi_H_E ( iV ), '>>> Chi_H_E' )
! end if
! end associate !-- R, etc.
! end associate !-- Y_E_V, etc.
! end associate !-- Y_E, etc.

      class default
        call Show ( 'Chart type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_NM_G_1D__Form', 'module', CONSOLE % ERROR )
        call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- C
    end do !-- iC

    end associate !-- S_R, etc.

  end subroutine SolveUpdateImplicit


  subroutine SetBalancedIndices &
               ( R, F, iMomentum_R, iMomentum_F, iEnergy_R, iEnergy_F, &
                 iNumber_R, iNumber_F )

    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R
    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    integer ( KDI ), dimension ( 3 ), intent ( out ) :: &
      iMomentum_R, &
      iMomentum_F
    integer ( KDI ), intent ( out ) :: &
      iEnergy_R, &
      iEnergy_F, &
      iNumber_R, &
      iNumber_F

    call Search &
           ( R % iaBalanced, R % ENERGY_DENSITY_B, iEnergy_R )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_1, iMomentum_R ( 1 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_2, iMomentum_R ( 2 ) )
    call Search &
           ( R % iaBalanced, R % MOMENTUM_DENSITY_B_D_3, iMomentum_R ( 3 ) )
    call Search &
           ( R % iaBalanced, R % NUMBER_DENSITY_B, iNumber_R )

    call Search &
           ( F % iaBalanced, F % ENERGY_DENSITY_B, iEnergy_F )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_1, iMomentum_F ( 1 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_2, iMomentum_F ( 2 ) )
    call Search &
           ( F % iaBalanced, F % MOMENTUM_DENSITY_D_3, iMomentum_F ( 3 ) )
    call Search &
           ( F % iaBalanced, F % ELECTRON_DENSITY_B, iNumber_F )

  end subroutine SetBalancedIndices


  subroutine SetStoragePointers_F ( F, Y_I, KK, iC, F_V, Y_I_V, KK_V )

    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      F, Y_I, KK
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      F_V, Y_I_V, KK_V

      F_V  =>    F % Storage ( iC ) % Value
    Y_I_V  =>  Y_I % Storage ( iC ) % Value
     KK_V  =>   KK % Storage ( iC ) % Value

  end subroutine SetStoragePointers_F


  subroutine SetStoragePointers_R &
               ( I, R, Y_I, KK, iC, I_V, R_V, Y_I_V, KK_V )

    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      I, R, Y_I, KK
    integer ( KDI ), intent ( in ) :: &
      iC
    real ( KDR ), dimension ( :, : ), pointer, intent ( out ) :: &
      I_V, R_V, Y_I_V, KK_V

      I_V  =>    I % Storage ( iC ) % Value
      R_V  =>    R % Storage ( iC ) % Value
    Y_I_V  =>  Y_I % Storage ( iC ) % Value
     KK_V  =>   KK % Storage ( iC ) % Value

  end subroutine SetStoragePointers_R


  subroutine SetFieldPointers_FS_B &  !-- FieldSet_Balanced
               ( FS_V, iMomentum, iEnergy, iNumber, &
                 FS_S_1, FS_S_2, FS_S_3, FS_E, FS_D )
    
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      FS_V
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      iMomentum
    integer ( KDI ), intent ( in ) :: &
      iEnergy, &
      iNumber
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      FS_S_1, FS_S_2, FS_S_3, &
      FS_E, &
      FS_D

    FS_S_1  =>  FS_V ( :, iMomentum ( 1 ) ) 
    FS_S_2  =>  FS_V ( :, iMomentum ( 2 ) ) 
    FS_S_3  =>  FS_V ( :, iMomentum ( 3 ) ) 
    
    FS_E    =>  FS_V ( :, iEnergy )
    FS_D    =>  FS_V ( :, iNumber )

  end subroutine SetFieldPointers_FS_B


  subroutine SetFieldPointers_F &
               ( F, F_V, E, S_1, S_2, S_3, D, DB, &
                 M, N, V_1, V_2, V_3, EC, YE, T, P, SB, SS, &
                 X_AA, X_p, X_n, X_A, Z, A, Mu_e, Mu_n_p, Mu_p, Mu_n, G )

    class ( Fluid_P_HN_Form ), intent ( in ) :: &
      F
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      F_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      E, S_1, S_2, S_3, D, DB, &
      M, N, V_1, V_2, V_3, EC, YE, &
      T, P, SB, SS, X_AA, X_p, X_n, X_A, &
      Z, A, Mu_e, Mu_n_p, Mu_p, Mu_n, G

      E     =>  F_V ( :, F % ENERGY_DENSITY_B )
      S_1   =>  F_V ( :, F % MOMENTUM_DENSITY_D_1 )
      S_2   =>  F_V ( :, F % MOMENTUM_DENSITY_D_2 )
      S_3   =>  F_V ( :, F % MOMENTUM_DENSITY_D_3 )
      D     =>  F_V ( :, F % ELECTRON_DENSITY_B )
      DB    =>  F_V ( :, F % BARYON_DENSITY_B )
      
      M       =>  F_V ( :, F % BARYON_MASS )
      N       =>  F_V ( :, F % BARYON_DENSITY_C )
      V_1     =>  F_V ( :, F % VELOCITY_U_1 )
      V_2     =>  F_V ( :, F % VELOCITY_U_1 )
      V_3     =>  F_V ( :, F % VELOCITY_U_1 )
      EC      =>  F_V ( :, F % ENERGY_DENSITY_C )
      YE      =>  F_V ( :, F % ELECTRON_FRACTION )
      T       =>  F_V ( :, F % TEMPERATURE )
      P       =>  F_V ( :, F % PRESSURE )
      SB      =>  F_V ( :, F % ENTROPY_PER_BARYON )
      SS      =>  F_V ( :, F % SOUND_SPEED )
      X_AA    =>  F_V ( :, F % MASS_FRACTION_ALPHA )
      X_p     =>  F_V ( :, F % MASS_FRACTION_PROTON )
      X_n     =>  F_V ( :, F % MASS_FRACTION_NEUTRON )
      X_A     =>  F_V ( :, F % MASS_FRACTION_HEAVY )
      Z       =>  F_V ( :, F % ATOMIC_NUMBER_HEAVY )
      A       =>  F_V ( :, F % MASS_NUMBER_HEAVY )
      Mu_e    =>  F_V ( :, F % CHEMICAL_POTENTIAL_E )
      Mu_n_p  =>  F_V ( :, F % CHEMICAL_POTENTIAL_N_P )
      Mu_p    =>  F_V ( :, F % CHEMICAL_POTENTIAL_P )
      Mu_n    =>  F_V ( :, F % CHEMICAL_POTENTIAL_N )
      G       =>  F_V ( :, F % ADIABATIC_INDEX )

  end subroutine SetFieldPointers_F


  subroutine SetFieldPointers_R &
               ( R, R_V, J, H_1, H_2, H_3, N, E, S_1, S_2, S_3, D, &
                 J_Eq, N_Eq, J_Rd, N_Rd, T_Nu, Eta_Nu, E_Ave, F_Ave )

    class ( NeutrinoMoments_G_Form ), intent ( in ) :: &
      R
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      R_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
      J, H_1, H_2, H_3, N, &
      E, S_1, S_2, S_3, D, &
      J_Eq, N_Eq, &
      J_Rd, N_Rd, &
      T_Nu, Eta_Nu, &
      E_Ave, F_Ave

      J       =>  R_V ( :, R % ENERGY_DENSITY_C )
      H_1     =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_1 )
      H_2     =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_2 )
      H_3     =>  R_V ( :, R % MOMENTUM_DENSITY_C_U_3 )
      N       =>  R_V ( :, R % NUMBER_DENSITY_C )
      E       =>  R_V ( :, R % ENERGY_DENSITY_B )
      S_1     =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_1 )
      S_2     =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_2 )
      S_3     =>  R_V ( :, R % MOMENTUM_DENSITY_B_D_3 )
      D       =>  R_V ( :, R % NUMBER_DENSITY_B )
      J_Eq    =>  R_V ( :, R % ENERGY_DENSITY_C_EQ )
      N_Eq    =>  R_V ( :, R % NUMBER_DENSITY_C_EQ )
      J_Rd    =>  R_V ( :, R % ENERGY_DENSITY_C_RD )
      N_Rd    =>  R_V ( :, R % NUMBER_DENSITY_C_RD )
      T_Nu    =>  R_V ( :, R % TEMPERATURE_GREY )
      Eta_Nu  =>  R_V ( :, R % DEGENERACY_GREY )
      E_Ave   =>  R_V ( :, R % ENERGY_AVERAGE )
      F_Ave   =>  R_V ( :, R % OCCUPANCY_AVERAGE )

  end subroutine SetFieldPointers_R


  subroutine SetFieldPointers_I &
               ( I, I_V, &
                 Xi_J, Xi_H, Xi_N, &
                 Xi_J_EA_N, Xi_J_EA_A, Xi_J_P_EP, Xi_J_S_EP, &
                 Chi_J, Chi_H, Chi_N, &
                 Chi_J_EA_N, Chi_J_EA_A, Chi_J_P_EP, Chi_J_S_EP, &
                 Chi_H_S_N, Chi_H_S_A, Chi_H_S_EP )

    class ( Interactions_NM_G_Form ), intent ( in ) :: &
      I
    real ( KDR ), dimension ( :, : ), intent ( in ), target :: &
      I_V
    real ( KDR ), dimension ( : ), intent ( out ), pointer :: &
       Xi_J,  Xi_H,  Xi_N, &
       Xi_J_EA_N, Xi_J_EA_A, Xi_J_P_EP, Xi_J_S_EP, &
      Chi_J, Chi_H, Chi_N, &
      Chi_J_EA_N, Chi_J_EA_A, Chi_J_P_EP, Chi_J_S_EP, &
      Chi_H_S_N, Chi_H_S_A, Chi_H_S_EP

     Xi_J         =>  I_V ( :, I % EMISSIVITY_J )
     Xi_H         =>  I_V ( :, I % EMISSIVITY_H )
     Xi_N         =>  I_V ( :, I % EMISSIVITY_N )
     Xi_J_EA_N    =>  I_V ( :, I % EMISSIVITY_J_EA_N )
     Xi_J_EA_A    =>  I_V ( :, I % EMISSIVITY_J_EA_A )
     Xi_J_P_EP    =>  I_V ( :, I % EMISSIVITY_J_P_EP )
     Xi_J_S_EP    =>  I_V ( :, I % EMISSIVITY_J_S_EP )
     
    Chi_J         =>  I_V ( :, I % OPACITY_J )
    Chi_H         =>  I_V ( :, I % OPACITY_H )
    Chi_N         =>  I_V ( :, I % OPACITY_N )
    Chi_J_EA_N    =>  I_V ( :, I % OPACITY_J_EA_N )
    Chi_J_EA_A    =>  I_V ( :, I % OPACITY_J_EA_A )
    Chi_J_P_EP    =>  I_V ( :, I % OPACITY_J_P_EP )
    Chi_J_S_EP    =>  I_V ( :, I % OPACITY_J_S_EP )
    Chi_H_S_N     =>  I_V ( :, I % OPACITY_H_S_N )
    Chi_H_S_A     =>  I_V ( :, I % OPACITY_H_S_A )
    Chi_H_S_EP    =>  I_V ( :, I % OPACITY_H_S_EP )
                     
  end subroutine SetFieldPointers_I


end module Step_RK_NM_G_1D_C__Form
