#include "Preprocessor"

submodule ( Interactions_NM_G__Form ) Interactions_NM_G__Kernel

  use Basics 
  
  implicit none

  real ( KDR ), parameter :: &
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
  real ( KDR ), parameter :: &
     Pi      =  ( CONSTANT % PI ), &
     Pi_2    =  ( CONSTANT % PI ) ** 2, &
     Pi_3    =  ( CONSTANT % PI ) ** 3, &
     Pi_5    =  ( CONSTANT % PI ) ** 5, &
      amu    =  ( CONSTANT % ATOMIC_MASS_UNIT ), &
        Q    =  ( CONSTANT % NEUTRON_MASS )  -  ( CONSTANT % PROTON_MASS ), &
    G_F_2    =  ( CONSTANT % FERMI_COUPLING ) ** 2, &
    g_A_2    =  ( CONSTANT % NEUTRON_AXIAL_COUPLING ) ** 2, &
    S_2_T_W  =  CONSTANT % SIN_2_WEINBERG

contains


  module procedure ComputeInteractions_E_S_Kernel

    !-- ComputeInteractions_Electron_Single_Kernel

    !--   Spectral & Equilibrium
    call Compute_SP_S_Kernel &
           ( T_Nu, Eta_Nu, E_Ave, F_Ave, &
             J, N, nSpecies = 1, iV = iV )
    call Compute_Eq_E_S_Kernel &
           ( J_Eq, N_Eq, J_Rd, N_Rd, J, N, &
             T_F, Mu_e_F, Mu_n_p_F, &
             Sign = +1.0_KDR, iV = iV )
            
    !-- Emission / Absorption
    call Compute_EA_E_S_Kernel &
           ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
             J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
             M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, Mu_e_F, Mu_n_p_F, &
             Rho_DB, iV )

    !-- Pair emission
    call Compute_P_S_Kernel &
           ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_P_EP, Chi_J_P_EP, &
             J_Eq, N_Eq, T_nu, Eta_nu, T_nuB, Eta_nuB, &
             M_F, N_F, T_F, Mu_e_F, &
             Sign = +1, nSpecies = 1, Rho_DB = Rho_DB, iV = iV )

    !-- Scattering on nucleons and nuclei
    call Compute_S_B_S_Kernel &
           ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
             M_F, N_F, X_p_F, X_n_F, X_A_F, Z_F, A_F, iV )

    !-- Scattering on electrons and positrons
    call Compute_S_EP_E_EB_S_Kernel &
           ( Xi_J, Chi_J, Chi_H, &
             Xi_J_S_EP, Chi_J_S_EP, Chi_H_S_EP, &
             J_Eq, T_nu, Eta_nu, M_F, N_F, T_F, Mu_e_F, &
             Sign = +1, Rho_DB = Rho_DB, iV = iV )

  end procedure ComputeInteractions_E_S_Kernel


  module procedure ComputeInteractions_EB_S_Kernel

    !-- ComputeInteractions_ElectronBar_Single_Kernel

    !--   Spectral & Equilibrium
    call Compute_SP_S_Kernel &
           ( T_Nu, Eta_Nu, E_Ave, F_Ave, &
             J, N, nSpecies = 1, iV = iV )
    call Compute_Eq_E_S_Kernel &
           ( J_Eq, N_Eq, J_Rd, N_Rd, J, N, &
             T_F, Mu_e_F, Mu_n_p_F, &
             Sign = -1.0_KDR, iV = iV )
            
    !-- Emission / Absorption
    call Compute_EA_EB_S_Kernel &
           ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
             J_Eq, N_Eq, T_nu, Eta_nu, E_Ave, F_Ave, &
             M_F, N_F, T_F, X_n_F, X_p_F, Mu_e_F, Mu_n_p_F, &
             Rho_DB, iV )

    !-- Pair emission
    call Compute_P_S_Kernel &
           ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_P_EP, Chi_J_P_EP, &
             J_Eq, N_Eq, T_nu, Eta_nu, T_nuB, Eta_nuB, &
             M_F, N_F, T_F, Mu_e_F, &
             Sign = +1, nSpecies = 1, Rho_DB = Rho_DB, iV = iV )

    !-- Scattering on nucleons and nuclei
    call Compute_S_B_S_Kernel &
           ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
             M_F, N_F, X_p_F, X_n_F, X_A_F, Z_F, A_F, iV )

    !-- Scattering on electrons and positrons
    call Compute_S_EP_E_EB_S_Kernel &
           ( Xi_J, Chi_J, Chi_H, &
             Xi_J_S_EP, Chi_J_S_EP, Chi_H_S_EP, &
             J_Eq, T_nu, Eta_nu, M_F, N_F, T_F, Mu_e_F, &
             Sign = -1, Rho_DB = Rho_DB, iV = iV )

  end procedure ComputeInteractions_EB_S_Kernel


  module procedure ComputeInteractions_HL_S_Kernel

    !-- ComputeInteractions_HeavyLepton_Single_Kernel

    !--   Spectral & Equilibrium
    call Compute_SP_S_Kernel &
           ( T_Nu, Eta_Nu, E_Ave, F_Ave, &
             J, N, nSpecies = 4, iV = iV )
    call Compute_Eq_HL_S_Kernel &
           ( J_Eq, N_Eq, J_RD, N_RD, J, N, T_F, &
               nSpecies = 4, iV = iV )

    !-- Emission / Absorption
    call Compute_EA_HL_S_Kernel &
           ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_EA_N, Xi_J_EA_A, Chi_J_EA_N, Chi_J_EA_A, &
             Chi_H_S_N, Chi_H_S_A, &
             iV = iV )

    !-- Pair emission
    call Compute_P_S_Kernel &
           ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, &
             Xi_J_P_EP, Chi_J_P_EP, &
             J_Eq, N_Eq, T_nu, Eta_nu, T_nuB, Eta_nuB, &
             M_F, N_F, T_F, Mu_e_F, &
             Sign = -1, nSpecies = 4, Rho_DB = Rho_DB, iV = iV )

    !-- Scattering on nucleons and nuclei
    call Compute_S_B_S_Kernel &
           ( Chi_H, Chi_H_S_N, Chi_H_S_A, T_nu, Eta_nu, &
             M_F, N_F, X_p_F, X_n_F, X_A_F, Z_F, A_F, iV )

    !-- Scattering on electrons and positrons
    call Compute_S_EP_HL_S_Kernel &
           ( Xi_J, Chi_J, Chi_H, &
             Xi_J_S_EP, Chi_J_S_EP, Chi_H_S_EP, &
             J_Eq, T_nu, Eta_nu, M_F, N_F, T_F, Mu_e_F, &
             nSpecies = 4, Rho_DB = Rho_DB, iV = iV )

  end procedure ComputeInteractions_HL_S_Kernel


  module procedure Compute_EA_E_A_Kernel

    !-- Compute_EmissionAbsorption_Electron_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor_p, Factor_n, Factor_A, Dlta, &
      N_p, N_n, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, F_e, &
      Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp, &
      Fermi_2_nu,  Fermi_3_nu,  Fermi_4_nu,  Fermi_5_nu, &
       Xi_J_p,  Xi_H_p,  Xi_N_p, &
      Chi_J_n, Chi_H_n, Chi_N_n, &
       Xi_J_A,  Xi_H_A,  Xi_N_A, &
      Chi_J_A, Chi_H_A, Chi_N_A
    logical ( KDL ) :: &
      UseDevice      

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    Factor_p  =  G_F_2 / ( 2 * Pi_3 )  *  ( 1  +  3 * g_A_2 )
    Factor_n  =  G_F_2 / Pi            *  ( 1  +  3 * g_A_2 )
    Factor_A  =  G_F_2 / ( 2 * Pi_3 )  *  ( 2.0_KDR / 7.0_KDR )  *  g_A_2
      Dlta    =  3.0_KDR !-- MeV

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor_p, Factor_n, Factor_A, Dlta ) &
      !$OMP private &
      !$OMP   ( N_p, N_n, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, F_e, &
      !$OMP     Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      !$OMP     Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp, &
      !$OMP     Fermi_2_nu,  Fermi_3_nu,  Fermi_4_nu,  Fermi_5_nu, &
      !$OMP     Xi_J_p,  Xi_H_p,  Xi_N_p, &
      !$OMP     Chi_J_n, Chi_H_n, Chi_N_n, &
      !$OMP     Xi_J_A,  Xi_H_A,  Xi_N_A, &
      !$OMP     Chi_J_A, Chi_H_A, Chi_N_A )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        !-- e- + p  ->  n + nu_e 

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
        Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
        Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
        Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

        Xi_J_p  =  Factor_p  *  N_p  *  T ( iV ) ** 4  &
                   *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                        +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                        +  Q ** 2            *  Fermi_3_e_Q )!  &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        Xi_H_p  =  0.0_KDR

        Xi_N_p  =  Factor_p  *  N_p  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                        +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                        +  Q ** 2            *  Fermi_2_e_Q )!  &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then

          !-- nu_e + n  ->  p + e-, detailed balance

          Chi_J_n  =  Xi_J_p / J_Eq ( iV )

          Chi_H_n  =  Chi_J_n

          Chi_N_n  =  Xi_N_p / N_Eq ( iV )

        else

          !-- nu_e + n  ->  p + e-, direct

          N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

          Fermi_2_nu  =  Fermi_2 ( Eta_nu ( iV ) )
          Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
          Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )
          Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

          F_e  =  1.  &
                  /  ( 1. +  &
                       exp ( ( E_Ave ( iV )  -  Mu_e ( iV ) ) / T ( iV ) ) )

          Chi_J_n  =  Factor_n  *  N_n  /  max ( Fermi_3_nu, SqrtTiny )  &
                     *  (    T_nu ( iV ) ** 2     *  Fermi_5_nu  &
                          +  2 * Q * T_nu ( iV )  *  Fermi_4_nu  &
                          +  Q ** 2               *  Fermi_3_nu )!  &
  !                   *  ( 1.0_KDR  -  F_e )

          Chi_H_n  =  Chi_J_n

          Chi_N_n  =  Factor_n  *  N_n  /  max ( Fermi_2_nu, SqrtTiny )  &
                     *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                          +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                          +  Q ** 2               *  Fermi_2_nu )!  &
  !                   *  ( 1.0_KDR  -  F_e )

        end if

        !-- e- + A  ->  A' + nu_e 

        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )
        
        if ( Z ( iV )  <=  20.0_KDR ) then
          N_p_Z  =  0.0_KDR
        else if ( Z ( iV )  <=  28.0_KDR ) then
          N_p_Z  =  Z ( iV )  -  20.0_KDR
        else
          N_p_Z  =  8.0_KDR
        end if

        if ( A ( iV )  -  Z ( iV )  <=  34.0_KDR ) then
          N_h_N  =  6.0_KDR
        else if ( A ( iV )  -  Z ( iV )  <=  40.0_KDR ) then
          N_h_N  =  40.0_KDR  -  ( A ( iV )  -  Z ( iV ) )
        else
          N_h_N  =  0.0_KDR
        end if

        Qp  =  Mu_n_p ( iV )  +  Dlta

        Eta_e_Qp  =  ( Mu_e ( iV )  -  Qp )  /  T ( iV )

        Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
        Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
        Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
        Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

        Xi_J_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 4  &
                   *  (    T ( iV ) ** 2      *  Fermi_5_e_Qp  &
                        +  2 * Qp * T ( iV )  *  Fermi_4_e_Qp  &
                        +  Qp ** 2            *  Fermi_3_e_Qp )

        Xi_H_A  =  0.0_KDR

        Xi_N_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 2      *  Fermi_4_e_Qp  &
                        +  2 * Qp * T ( iV )  *  Fermi_3_e_Qp  &
                        +  Qp ** 2            *  Fermi_2_e_Qp )

        !-- nu_e + A'  ->  A + e-, detailed balance

        Chi_J_A  =  Xi_J_A / J_Eq ( iV )

        Chi_H_A  =  Chi_J_A

        Chi_N_A  =  Xi_N_A / N_Eq ( iV )

        !-- Total

         Xi_J_EA_N ( iV )  =   Xi_J_p
        Chi_J_EA_N ( iV )  =  Chi_J_n

         Xi_J_EA_A ( iV )  =   Xi_J_A
        Chi_J_EA_A ( iV )  =  Chi_J_A

        Xi_J ( iV )  =  Xi_J_p  +  Xi_J_A
        Xi_H ( iV )  =  Xi_H_p  +  Xi_H_A
        Xi_N ( iV )  =  Xi_N_p  +  Xi_N_A
        
        Chi_J ( iV )  =  Chi_J_n  +  Chi_J_A
        Chi_H ( iV )  =  Chi_H_n  +  Chi_H_A
        Chi_N ( iV )  =  Chi_N_n  +  Chi_N_A

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor_p, Factor_n, Factor_A, Dlta ) &
      !$OMP private &
      !$OMP   ( N_p, N_n, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, F_e, &
      !$OMP     Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      !$OMP     Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp, &
      !$OMP     Fermi_2_nu,  Fermi_3_nu,  Fermi_4_nu,  Fermi_5_nu, &
      !$OMP     Xi_J_p,  Xi_H_p,  Xi_N_p, &
      !$OMP     Chi_J_n, Chi_H_n, Chi_N_n, &
      !$OMP     Xi_J_A,  Xi_H_A,  Xi_N_A, &
      !$OMP     Chi_J_A, Chi_H_A, Chi_N_A )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        !-- e- + p  ->  n + nu_e 

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
        Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
        Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
        Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

        Xi_J_p  =  Factor_p  *  N_p  *  T ( iV ) ** 4  &
                   *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                        +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                        +  Q ** 2            *  Fermi_3_e_Q )!  &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        Xi_H_p  =  0.0_KDR

        Xi_N_p  =  Factor_p  *  N_p  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                        +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                        +  Q ** 2            *  Fermi_2_e_Q )!  &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then

          !-- nu_e + n  ->  p + e-, detailed balance

          Chi_J_n  =  Xi_J_p / J_Eq ( iV )

          Chi_H_n  =  Chi_J_n

          Chi_N_n  =  Xi_N_p / N_Eq ( iV )

        else

          !-- nu_e + n  ->  p + e-, direct

          N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

          Fermi_2_nu  =  Fermi_2 ( Eta_nu ( iV ) )
          Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
          Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )
          Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

          F_e  =  1.  &
                  /  ( 1. +  &
                       exp ( ( E_Ave ( iV )  -  Mu_e ( iV ) ) / T ( iV ) ) )

          Chi_J_n  =  Factor_n  *  N_n  /  Fermi_3_nu  &
                     *  (    T_nu ( iV ) ** 2     *  Fermi_5_nu  &
                          +  2 * Q * T_nu ( iV )  *  Fermi_4_nu  &
                          +  Q ** 2               *  Fermi_3_nu )!  &
  !                   *  ( 1.0_KDR  -  F_e )

          Chi_H_n  =  Chi_J_n

          Chi_N_n  =  Factor_n  *  N_n  /  Fermi_2_nu  &
                     *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                          +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                          +  Q ** 2               *  Fermi_2_nu )!  &
  !                   *  ( 1.0_KDR  -  F_e )

        end if

        !-- e- + A  ->  A' + nu_e 

        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )
        
        if ( Z ( iV )  <=  20.0_KDR ) then
          N_p_Z  =  0.0_KDR
        else if ( Z ( iV )  <=  28.0_KDR ) then
          N_p_Z  =  Z ( iV )  -  20.0_KDR
        else
          N_p_Z  =  8.0_KDR
        end if

        if ( A ( iV )  -  Z ( iV )  <=  34.0_KDR ) then
          N_h_N  =  6.0_KDR
        else if ( A ( iV )  -  Z ( iV )  <=  40.0_KDR ) then
          N_h_N  =  40.0_KDR  -  ( A ( iV )  -  Z ( iV ) )
        else
          N_h_N  =  0.0_KDR
        end if

        Qp  =  Mu_n_p ( iV )  +  Dlta

        Eta_e_Qp  =  ( Mu_e ( iV )  -  Qp )  /  T ( iV )

        Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
        Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
        Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
        Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

        Xi_J_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 4  &
                   *  (    T ( iV ) ** 2      *  Fermi_5_e_Qp  &
                        +  2 * Qp * T ( iV )  *  Fermi_4_e_Qp  &
                        +  Qp ** 2            *  Fermi_3_e_Qp )

        Xi_H_A  =  0.0_KDR

        Xi_N_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 2      *  Fermi_4_e_Qp  &
                        +  2 * Qp * T ( iV )  *  Fermi_3_e_Qp  &
                        +  Qp ** 2            *  Fermi_2_e_Qp )

        !-- nu_e + A'  ->  A + e-, detailed balance

        Chi_J_A  =  Xi_J_A / J_Eq ( iV )

        Chi_H_A  =  Chi_J_A

        Chi_N_A  =  Xi_N_A / N_Eq ( iV )

        !-- Total

         Xi_J_EA_N ( iV )  =   Xi_J_p
        Chi_J_EA_N ( iV )  =  Chi_J_n

         Xi_J_EA_A ( iV )  =   Xi_J_A
        Chi_J_EA_A ( iV )  =  Chi_J_A

        Xi_J ( iV )  =  Xi_J_p  +  Xi_J_A
        Xi_H ( iV )  =  Xi_H_p  +  Xi_H_A
        Xi_N ( iV )  =  Xi_N_p  +  Xi_N_A
        
        Chi_J ( iV )  =  Chi_J_n  +  Chi_J_A
        Chi_H ( iV )  =  Chi_H_n  +  Chi_H_A
        Chi_N ( iV )  =  Chi_N_n  +  Chi_N_A

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_E_A_Kernel


  module procedure Compute_EA_E_S_Kernel

    !-- Compute_EmissionAbsorption_Electron_Single_Kernel

    real ( KDR ) :: &
      Factor_p, Factor_n, Factor_A, Dlta, &
      N_p, N_n, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, F_e, &
      Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp, &
      Fermi_2_nu,  Fermi_3_nu,  Fermi_4_nu,  Fermi_5_nu, &
       Xi_J_p,  Xi_H_p,  Xi_N_p, &
      Chi_J_n, Chi_H_n, Chi_N_n, &
       Xi_J_A,  Xi_H_A,  Xi_N_A, &
      Chi_J_A, Chi_H_A, Chi_N_A

    !$OMP_DECLARE_TARGET
    
    Factor_p  =  G_F_2 / ( 2 * Pi_3 )  *  ( 1  +  3 * g_A_2 )
    Factor_n  =  G_F_2 / Pi            *  ( 1  +  3 * g_A_2 )
    Factor_A  =  G_F_2 / ( 2 * Pi_3 )  *  ( 2.0_KDR / 7.0_KDR )  *  g_A_2
      Dlta    =  3.0_KDR !-- MeV

    if ( T ( iV ) == 0.0_KDR ) &
      return

    !-- e- + p  ->  n + nu_e 

    N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

    Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

    Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
    Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
    Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
    Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

    Xi_J_p  =  Factor_p  *  N_p  *  T ( iV ) ** 4  &
               *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                    +  Q ** 2            *  Fermi_3_e_Q )!  &
!               *  ( 1.0_KDR  -  F_Ave ( iV ) )

    Xi_H_p  =  0.0_KDR

    Xi_N_p  =  Factor_p  *  N_p  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                    +  Q ** 2            *  Fermi_2_e_Q )!  &
!               *  ( 1.0_KDR  -  F_Ave ( iV ) )

    if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then

      !-- nu_e + n  ->  p + e-, detailed balance

      Chi_J_n  =  Xi_J_p / J_Eq ( iV )

      Chi_H_n  =  Chi_J_n

      Chi_N_n  =  Xi_N_p / N_Eq ( iV )

    else

      !-- nu_e + n  ->  p + e-, direct

      N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

      Fermi_2_nu  =  Fermi_2 ( Eta_nu ( iV ) )
      Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
      Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )
      Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

      F_e  =  1.  &
              /  ( 1. +  &
                   exp ( ( E_Ave ( iV )  -  Mu_e ( iV ) ) / T ( iV ) ) )

      Chi_J_n  =  Factor_n  *  N_n  /  max ( Fermi_3_nu, SqrtTiny )  &
                 *  (    T_nu ( iV ) ** 2     *  Fermi_5_nu  &
                      +  2 * Q * T_nu ( iV )  *  Fermi_4_nu  &
                      +  Q ** 2               *  Fermi_3_nu )!  &
!                  *  ( 1.0_KDR  -  F_e )

      Chi_H_n  =  Chi_J_n

      Chi_N_n  =  Factor_n  *  N_n  /  max ( Fermi_2_nu, SqrtTiny )  &
                 *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                      +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                      +  Q ** 2               *  Fermi_2_nu )!  &
!                  *  ( 1.0_KDR  -  F_e )

    end if

    !-- e- + A  ->  A' + nu_e 

    N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
            /  ( max ( A ( iV ), SqrtTiny ) * amu )
    
    if ( Z ( iV )  <=  20.0_KDR ) then
      N_p_Z  =  0.0_KDR
    else if ( Z ( iV )  <=  28.0_KDR ) then
      N_p_Z  =  Z ( iV )  -  20.0_KDR
    else
      N_p_Z  =  8.0_KDR
    end if

    if ( A ( iV )  -  Z ( iV )  <=  34.0_KDR ) then
      N_h_N  =  6.0_KDR
    else if ( A ( iV )  -  Z ( iV )  <=  40.0_KDR ) then
      N_h_N  =  40.0_KDR  -  ( A ( iV )  -  Z ( iV ) )
    else
      N_h_N  =  0.0_KDR
    end if

    Qp  =  Mu_n_p ( iV )  +  Dlta

    Eta_e_Qp  =  ( Mu_e ( iV )  -  Qp )  /  T ( iV )

    Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
    Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
    Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
    Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

    Xi_J_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 4  &
               *  (    T ( iV ) ** 2      *  Fermi_5_e_Qp  &
                    +  2 * Qp * T ( iV )  *  Fermi_4_e_Qp  &
                    +  Qp ** 2            *  Fermi_3_e_Qp )

    Xi_H_A  =  0.0_KDR

    Xi_N_A  =  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2      *  Fermi_4_e_Qp  &
                    +  2 * Qp * T ( iV )  *  Fermi_3_e_Qp  &
                    +  Qp ** 2            *  Fermi_2_e_Qp )

    !-- nu_e + A'  ->  A + e-, detailed balance

    Chi_J_A  =  Xi_J_A / J_Eq ( iV )

    Chi_H_A  =  Chi_J_A

    Chi_N_A  =  Xi_N_A / N_Eq ( iV )

    !-- Total

     Xi_J_EA_N ( iV )  =   Xi_J_p
    Chi_J_EA_N ( iV )  =  Chi_J_n

     Xi_J_EA_A ( iV )  =   Xi_J_A
    Chi_J_EA_A ( iV )  =  Chi_J_A

    Xi_J ( iV )  =  Xi_J_p  +  Xi_J_A
    Xi_H ( iV )  =  Xi_H_p  +  Xi_H_A
    Xi_N ( iV )  =  Xi_N_p  +  Xi_N_A
    
    Chi_J ( iV )  =  Chi_J_n  +  Chi_J_A
    Chi_H ( iV )  =  Chi_H_n  +  Chi_H_A
    Chi_N ( iV )  =  Chi_N_n  +  Chi_N_A

  end procedure Compute_EA_E_S_Kernel


  module procedure Compute_EA_EB_A_Kernel

    !-- Compute_EmissionAbsorption_ElectronBar_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor_n, Factor_p, &
      N_n, N_p, Eta_e, &
      Fermi_2_e,  Fermi_3_e,  Fermi_4_e,  Fermi_5_e, &
      Fermi_2_nu, Fermi_3_nu, Fermi_4_nu, Fermi_5_nu, &
       Xi_J_n,  Xi_H_n,  Xi_N_n, &
      Chi_J_p, Chi_H_p, Chi_N_p
    logical ( KDL ) :: &
      UseDevice      

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    Factor_n  =  G_F_2 / ( 2 * Pi_3 )  *  ( 1  +  3 * g_A_2 )
    Factor_p  =  G_F_2 / Pi            *  ( 1  +  3 * g_A_2 )

    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor_n ) &
      !$OMP private ( N_n, Eta_e, &
      !$OMP           Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e, &
      !$OMP            Xi_J_n,  Xi_H_n,  Xi_N_n, &
      !$OMP           Chi_J_p, Chi_H_p, Chi_N_p )
      
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        !-- e+ + n  ->  p + nu_e_bar 

        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        Fermi_2_e  =  Fermi_2 ( -Eta_e )
        Fermi_3_e  =  Fermi_3 ( -Eta_e )
        Fermi_4_e  =  Fermi_4 ( -Eta_e )
        Fermi_5_e  =  Fermi_5 ( -Eta_e )

        Xi_J_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 3             *  Fermi_5_e  &
                        +  3 *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                        +  3 *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                        +  Q ** 3                    *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        Xi_H_n  =  0.0_KDR

        Xi_N_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then

          !-- nu_e_bar + p  ->  n + e+, detailed balance

          Chi_J_p  =  Xi_J_n / J_Eq ( iV )

          Chi_H_p  =  Chi_J_p

          Chi_N_p  =  Xi_N_n / N_Eq ( iV )

        else

          !-- nu_e_bar + p  ->  n + e+, direct

          N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

          Fermi_2_nu  =  Fermi_2 ( Eta_nu ( iV ) )
          Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
          Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )
          Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

          Chi_J_p  =  Factor_p  *  N_p  /  max ( Fermi_3_nu, SqrtTiny )  &
                      *  (    T_nu ( iV ) ** 2        *  Fermi_5_nu  &
                           +  3 *  Q  *  T_nu ( iV )  *  Fermi_4_nu  &
                           +  3 *  Q ** 2             *  Fermi_3_nu  &
                           +  Q ** 3  /  T_nu ( iV )  *  Fermi_2_nu )
                  
          Chi_H_p  =  Chi_J_p

          Chi_N_p  =  Factor_p  *  N_p  /  max ( Fermi_2_nu, SqrtTiny )  &
                      *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                           +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                           +  Q ** 2               *  Fermi_2_nu )

        end if

        !-- Total

         Xi_J_EA_N ( iV )  =   Xi_J_n
        Chi_J_EA_N ( iV )  =  Chi_J_p

         Xi_J_EA_A ( iV )  =  0.0_KDR
        Chi_J_EA_A ( iV )  =  0.0_KDR

        Xi_J ( iV )  =  Xi_J_n
        Xi_H ( iV )  =  Xi_H_n
        Xi_N ( iV )  =  Xi_N_n
        
        Chi_J ( iV )  =  Chi_J_p
        Chi_H ( iV )  =  Chi_H_p
        Chi_N ( iV )  =  Chi_N_p

      end do
      
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
      
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor_n ) &
      !$OMP private ( N_n, Eta_e, &
      !$OMP           Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e, &
      !$OMP            Xi_J_n,  Xi_H_n,  Xi_N_n, &
      !$OMP           Chi_J_p, Chi_H_p, Chi_N_p )
      
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        !-- e+ + n  ->  p + nu_e_bar 

        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        Fermi_2_e  =  Fermi_2 ( -Eta_e )
        Fermi_3_e  =  Fermi_3 ( -Eta_e )
        Fermi_4_e  =  Fermi_4 ( -Eta_e )
        Fermi_5_e  =  Fermi_5 ( -Eta_e )

        Xi_J_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
                   *  (    T ( iV ) ** 3             *  Fermi_5_e  &
                        +  3 *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                        +  3 *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                        +  Q ** 3                    *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        Xi_H_n  =  0.0_KDR

        Xi_N_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

        !-- nu_e_bar + p  ->  n + e+, detailed balance

        Chi_J_p  =  Xi_J_n / J_Eq ( iV )

        Chi_H_p  =  Chi_J_p

        Chi_N_p  =  Xi_N_n / N_Eq ( iV )

        !-- Total

         Xi_J_EA_N ( iV )  =   Xi_J_n
        Chi_J_EA_N ( iV )  =  Chi_J_p

         Xi_J_EA_A ( iV )  =  0.0_KDR
        Chi_J_EA_A ( iV )  =  0.0_KDR

        Xi_J ( iV )  =  Xi_J_n
        Xi_H ( iV )  =  Xi_H_n
        Xi_N ( iV )  =  Xi_N_n
        
        Chi_J ( iV )  =  Chi_J_p
        Chi_H ( iV )  =  Chi_H_p
        Chi_N ( iV )  =  Chi_N_p

      end do
      
      !$OMP end parallel do
      
    end if

  end procedure Compute_EA_EB_A_Kernel


  module procedure Compute_EA_EB_S_Kernel

    !-- Compute_EmissionAbsorption_ElectronBar_Single_Kernel

    real ( KDR ) :: &
      Factor_n, Factor_p, &
      N_n, N_p, Eta_e, &
      Fermi_2_e,  Fermi_3_e,  Fermi_4_e,  Fermi_5_e, &
      Fermi_2_nu, Fermi_3_nu, Fermi_4_nu, Fermi_5_nu, &
       Xi_J_n,  Xi_H_n,  Xi_N_n, &
      Chi_J_p, Chi_H_p, Chi_N_p
      
    !$OMP_DECLARE_TARGET

    Factor_n  =  G_F_2 / ( 2 * Pi_3 )  *  ( 1  +  3 * g_A_2 )
    Factor_p  =  G_F_2 / Pi            *  ( 1  +  3 * g_A_2 )

    !-- e+ + n  ->  p + nu_e_bar 

    N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu

    Eta_e  =  Mu_e ( iV )  /  T ( iV )

    Fermi_2_e  =  Fermi_2 ( -Eta_e )
    Fermi_3_e  =  Fermi_3 ( -Eta_e )
    Fermi_4_e  =  Fermi_4 ( -Eta_e )
    Fermi_5_e  =  Fermi_5 ( -Eta_e )

    Xi_J_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 3             *  Fermi_5_e  &
                    +  3 *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                    +  3 *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                    +  Q ** 3                    *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

    Xi_H_n  =  0.0_KDR

    Xi_N_n  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
           *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                +  2 * Q * T ( iV )  *  Fermi_3_e  &
                +  Q ** 2            *  Fermi_2_e )! &
!                   *  ( 1.0_KDR  -  F_Ave ( iV ) )

    if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then

      !-- nu_e_bar + p  ->  n + e+, detailed balance

      Chi_J_p  =  Xi_J_n / max ( J_Eq ( iV ), SqrtTiny )

      Chi_H_p  =  Chi_J_p

      Chi_N_p  =  Xi_N_n / max ( N_Eq ( iV ), SqrtTiny )

    else

      !-- nu_e_bar + p  ->  n + e+, direct

      N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

      Fermi_2_nu  =  Fermi_2 ( Eta_nu ( iV ) )
      Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
      Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )
      Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

      Chi_J_p  =  Factor_p  *  N_p  /  max ( Fermi_3_nu, SqrtTiny )  &
                  *  (    T_nu ( iV ) ** 2        *  Fermi_5_nu  &
                       +  3 *  Q  *  T_nu ( iV )  *  Fermi_4_nu  &
                       +  3 *  Q ** 2             *  Fermi_3_nu  &
                       +  Q ** 3  /  T_nu ( iV )  *  Fermi_2_nu )
              
      Chi_H_p  =  Chi_J_p

      Chi_N_p  =  Factor_p  *  N_p  /  max ( Fermi_2_nu, SqrtTiny )  &
                  *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                       +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                       +  Q ** 2               *  Fermi_2_nu )

    end if

    !-- Total

     Xi_J_EA_N ( iV )  =   Xi_J_n
    Chi_J_EA_N ( iV )  =  Chi_J_p

     Xi_J_EA_A ( iV )  =  0.0_KDR
    Chi_J_EA_A ( iV )  =  0.0_KDR

    Xi_J ( iV )  =  Xi_J_n
    Xi_H ( iV )  =  Xi_H_n
    Xi_N ( iV )  =  Xi_N_n
        
    Chi_J ( iV )  =  Chi_J_p
    Chi_H ( iV )  =  Chi_H_p
    Chi_N ( iV )  =  Chi_N_p

  end procedure Compute_EA_EB_S_Kernel


  module procedure Compute_EA_HL_A_Kernel

    !-- Compute_EmissionAbsorption_HeavyLepton_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Xi_J )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

         Xi_J_EA_N ( iV )  =  0.0_KDR
        Chi_J_EA_N ( iV )  =  0.0_KDR

         Xi_J_EA_A ( iV )  =  0.0_KDR
        Chi_J_EA_A ( iV )  =  0.0_KDR

        Chi_H_S_N ( iV )  =  0.0_KDR
        Chi_H_S_A ( iV )  =  0.0_KDR

        Xi_J ( iV )  =  0.0_KDR
        Xi_H ( iV )  =  0.0_KDR
        Xi_N ( iV )  =  0.0_KDR
        
        Chi_J ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  0.0_KDR
        Chi_N ( iV )  =  0.0_KDR

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

         Xi_J_EA_N ( iV )  =  0.0_KDR
        Chi_J_EA_N ( iV )  =  0.0_KDR

         Xi_J_EA_A ( iV )  =  0.0_KDR
        Chi_J_EA_A ( iV )  =  0.0_KDR

        Chi_H_S_N ( iV )  =  0.0_KDR
        Chi_H_S_A ( iV )  =  0.0_KDR

        Xi_J ( iV )  =  0.0_KDR
        Xi_H ( iV )  =  0.0_KDR
        Xi_N ( iV )  =  0.0_KDR
        
        Chi_J ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  0.0_KDR
        Chi_N ( iV )  =  0.0_KDR

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_HL_A_Kernel


  module procedure Compute_EA_HL_S_Kernel

    !-- Compute_EmissionAbsorption_HeavyLepton_Single_Kernel
    
    !$OMP_DECLARE_TARGET

     Xi_J_EA_N ( iV )  =  0.0_KDR
    Chi_J_EA_N ( iV )  =  0.0_KDR

     Xi_J_EA_A ( iV )  =  0.0_KDR
    Chi_J_EA_A ( iV )  =  0.0_KDR

    Chi_H_S_N ( iV )  =  0.0_KDR
    Chi_H_S_A ( iV )  =  0.0_KDR

    Xi_J ( iV )  =  0.0_KDR
    Xi_H ( iV )  =  0.0_KDR
    Xi_N ( iV )  =  0.0_KDR
    
    Chi_J ( iV )  =  0.0_KDR
    Chi_H ( iV )  =  0.0_KDR
    Chi_N ( iV )  =  0.0_KDR

  end procedure Compute_EA_HL_S_Kernel


  module procedure Compute_P_A_Kernel

    !-- Compute_Pair_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor, &
      Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu, Fermi_3_nuB, Fermi_4_nuB, &
       Xi_J_P,           Xi_N_P, &
      Chi_J_P, Chi_H_P, Chi_N_P
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Xi_J )

    Factor  =  nSpecies * G_F_2  /  ( 9.  *  Pi_5 )  &
               *  ( 1.  +  Sign * 4. * S_2_T_W  +  8. * S_2_T_W ** 2 )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor ) &
      !$OMP private ( Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, &
      !$OMP           Xi_J_P,           Xi_N_P, &
      !$OMP            Chi_J_P, Chi_H_P, Chi_N_P )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
        Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

        Xi_J_P  =  Factor  *  T ( iV ) ** 9  &
                   *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                              +  Fermi_3_eP * Fermi_4_eM )

        Xi_N_P  =  Factor  *  T ( iV ) ** 8  &
                   *  Fermi_3_eM * Fermi_3_eP

        if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

          !-- detailed balance

          Chi_J_P  =  Xi_J_P  /  max ( J_Eq ( iV ), SqrtTiny )
          Chi_H_P  =  Chi_J_P
          Chi_N_P  =  Xi_N_P  /  max ( N_Eq ( iV ), SqrtTiny )

        else

          !-- direct

          Chi_J_P  =  0.0_KDR
          Chi_H_P  =  0.0_KDR
          Chi_N_P  =  0.0_KDR

          ! Fermi_3_nu   =  Fermi_3 ( Eta_nu  ( iV ) )
          ! Fermi_4_nu   =  Fermi_4 ( Eta_nu  ( iV ) )
          ! Fermi_3_nuB  =  Fermi_3 ( Eta_nuB ( iV ) )
          ! Fermi_4_nuB  =  Fermi_4 ( Eta_nuB ( iV ) )

          ! Chi_J_P  =  Factor  *  T_nu ( iV ) ** 4  *  T_nuB ( iV ) ** 4  &
          !             *  0.5 * (    T_nu  ( iV ) * Fermi_3_nu  * Fermi_4_nu  &
          !                        +  T_nuB ( iV ) * Fermi_3_nuB * Fermi_4_nuB )

          ! Chi_H_P  =  Chi_J_P

          ! Chi_N_P  =  Factor  *  T_nu ( iV ) ** 4  *  T_nuB ( iV ) ** 4  &
          !             *  Fermi_3_nu * Fermi_3_nuB

        end if
          
        !-- Total

         Xi_J_P_EP ( iV )  =   Xi_J_P
        Chi_J_P_EP ( iV )  =  Chi_J_P

        Xi_J ( iV )  =  Xi_J ( iV )  +  Xi_J_P
        Xi_N ( iV )  =  Xi_N ( iV )  +  Xi_N_P

        Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_P
        Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_P
        Chi_N ( iV )  =  Chi_N ( iV )  +  Chi_N_P

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Rho_DB ) &
      !$OMP firstprivate ( Factor ) &
      !$OMP private ( Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, &
      !$OMP           Xi_J_P,           Xi_N_P, &
      !$OMP            Chi_J_P, Chi_H_P, Chi_N_P )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
        Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

        Xi_J_P  =  Factor  *  T ( iV ) ** 9  &
                   *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                              +  Fermi_3_eP * Fermi_4_eM )

        Xi_N_P  =  Factor  *  T ( iV ) ** 8  &
                   *  Fermi_3_eM * Fermi_3_eP

        if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 
          Chi_J_P  =  Xi_J_P  /  J_Eq ( iV )
          Chi_H_P  =  Xi_J_P  /  J_Eq ( iV )
          Chi_N_P  =  Xi_N_P  /  N_Eq ( iV )
        else
          Chi_J_P  =  0.0_KDR
          Chi_H_P  =  0.0_KDR
          Chi_N_P  =  0.0_KDR
        end if
          
        !-- Total

         Xi_J_P_EP ( iV )  =   Xi_J_P
        Chi_J_P_EP ( iV )  =  Chi_J_P

        Xi_J ( iV )  =  Xi_J ( iV )  +  Xi_J_P
        Xi_N ( iV )  =  Xi_N ( iV )  +  Xi_N_P

        Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_P
        Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_P
        Chi_N ( iV )  =  Chi_N ( iV )  +  Chi_N_P

      end do !-- iV
      !$OMP end parallel do
    end if
   
  end procedure Compute_P_A_Kernel


  module procedure Compute_P_S_Kernel

    !-- Compute_Pair_Single_Kernel

    real ( KDR ) :: &
      Factor, &
      Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu, Fermi_3_nuB, Fermi_4_nuB, &
       Xi_J_P,           Xi_N_P, &
      Chi_J_P, Chi_H_P, Chi_N_P
      
    !$OMP_DECLARE_TARGET
          
    Factor  =  nSpecies * G_F_2  /  ( 9.  *  Pi_5 )  &
               *  ( 1.  +  Sign * 4. * S_2_T_W  +  8. * S_2_T_W ** 2 )

    if ( T ( iV ) == 0.0_KDR ) &
      return

    Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
    Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

    Xi_J_P  =  Factor  *  T ( iV ) ** 9  &
               *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                          +  Fermi_3_eP * Fermi_4_eM )

    Xi_N_P  =  Factor  *  T ( iV ) ** 8  &
               *  Fermi_3_eM * Fermi_3_eP

    if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

      !-- detailed balance
      
      Chi_J_P  =  Xi_J_P  /  max ( J_Eq ( iV ), SqrtTiny )
      Chi_H_P  =  Chi_J_P
      Chi_N_P  =  Xi_N_P  /  max ( N_Eq ( iV ), SqrtTiny )

    else

      !-- direct

      Chi_J_P  =  0.0_KDR
      Chi_H_P  =  0.0_KDR
      Chi_N_P  =  0.0_KDR

      ! Fermi_3_nu   =  Fermi_3 ( Eta_nu  ( iV ) )
      ! Fermi_4_nu   =  Fermi_4 ( Eta_nu  ( iV ) )
      ! Fermi_3_nuB  =  Fermi_3 ( Eta_nuB ( iV ) )
      ! Fermi_4_nuB  =  Fermi_4 ( Eta_nuB ( iV ) )

      ! Chi_J_P  =  Factor  *  T_nu ( iV ) ** 4  *  T_nuB ( iV ) ** 4  &
      !             *  0.5 * (    T_nu  ( iV ) * Fermi_3_nu  * Fermi_4_nu  &
      !                        +  T_nuB ( iV ) * Fermi_3_nuB * Fermi_4_nuB )

      ! Chi_H_P  =  Chi_J_P

      ! Chi_N_P  =  Factor  *  T_nu ( iV ) ** 4  *  T_nuB ( iV ) ** 4  &
      !             *  Fermi_3_nu * Fermi_3_nuB

    end if
          
    !-- Total

     Xi_J_P_EP ( iV )  =   Xi_J_P
    Chi_J_P_EP ( iV )  =  Chi_J_P

    Xi_J ( iV )  =  Xi_J ( iV )  +  Xi_J_P
    Xi_N ( iV )  =  Xi_N ( iV )  +  Xi_N_P

    Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_P
    Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_P
    Chi_N ( iV )  =  Chi_N ( iV )  +  Chi_N_P

  end procedure Compute_P_S_Kernel


  module procedure Compute_S_B_A_Kernel

    !-- Compute_Scattering_Baryons_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor_p, Factor_n, Factor_A, &
      N_p, N_n, N_A, &
      Fermi_3_nu, Fermi_5_nu
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Chi_H )

    Factor_p  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                 *  ( ( 1. / 2.  -  2. * S_2_T_W ) ** 2  +  5. / 4. * g_A_2 )

    Factor_n  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                 *  ( 1. / 4.  +  5. / 4. * g_A_2 )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP firstprivate ( Factor_p, Factor_n ) &
      !$OMP private ( Factor_A, &
      !$OMP           N_p, N_n, N_A, &
      !$OMP           Fermi_3_nu, Fermi_5_nu )
      do iV = 1, nV

        Factor_A  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                     *  ( -  A ( iV )  /  2. &
                          +  Z ( iV )  *  ( 1.  -  2. * S_2_T_W ) ) ** 2  

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )

        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

        Chi_H_S_N ( iV )  &
          =  ( Factor_p * N_p  +  Factor_n * N_n )  &
             *  T_nu ( iV ) ** 2  *  Fermi_5_nu / max ( Fermi_3_nu, SqrtTiny )

        Chi_H_S_A ( iV )  &
          =  Factor_A * N_A  &
             *  T_nu ( iV ) ** 2  *  Fermi_5_nu / max ( Fermi_3_nu, SqrtTiny ) 

        Chi_H ( iV )  &
          =  Chi_H ( iV )  +  Chi_H_S_N ( iV )  +  Chi_H_S_A ( iV )

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP firstprivate ( Factor_p, Factor_n ) &
      !$OMP private ( Factor_A, &
      !$OMP           N_p, N_n, N_A, &
      !$OMP           Fermi_3_nu, Fermi_5_nu )
      do iV = 1, nV

        Factor_A  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                     *  (    A ( iV )  *  ( 1. / 2.  -  2. * S_2_T_W ) &
                          -  Z ( iV )  *  ( 1.  -  2. * S_2_T_W ) ) ** 2  

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )

        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

        !-- Elastic scattering on nucleons

        Chi_H_S_N ( iV )  &
          =  ( Factor_p * N_p  +  Factor_n * N_n )  &
             *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

        Chi_H_S_A ( iV )  &
          =  Factor_A * N_A  &
             *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

        Chi_H ( iV )  &
          =  Chi_H ( iV )  +  Chi_H_S_N ( iV )  +  Chi_H_S_A ( iV )

      end do
      !$OMP end parallel do
      
    end if

  end procedure Compute_S_B_A_Kernel


  module procedure Compute_S_B_S_Kernel

    !-- Compute_Scattering_Baryons_Single_Kernel

    real ( KDR ) :: &
      Factor_p, Factor_n, Factor_A, &
      N_p, N_n, N_A, &
      Fermi_3_nu, Fermi_5_nu
      
    !$OMP_DECLARE_TARGET

    Factor_p  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                 *  ( ( 1. / 2.  -  2. * S_2_T_W ) ** 2  +  5. / 4. * g_A_2 )

    Factor_n  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                 *  ( 1. / 4.  +  5. / 4. * g_A_2 )

    Factor_A  =  2.  *  G_F_2  /  ( 3. * Pi )  &
                 *  ( -  A ( iV )  /  2. &
                      +  Z ( iV )  *  ( 1.  -  2. * S_2_T_W ) ) ** 2  

    N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
    N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
    N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
            /  ( max ( A ( iV ), SqrtTiny ) * amu )

    Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
    Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

    Chi_H_S_N ( iV )  &
      =  ( Factor_p * N_p  +  Factor_n * N_n )  &
         *  T_nu ( iV ) ** 2  *  Fermi_5_nu / max ( Fermi_3_nu, SqrtTiny ) 

    Chi_H_S_A ( iV )  &
      =  Factor_A * N_A  &
         *  T_nu ( iV ) ** 2  *  Fermi_5_nu / max ( Fermi_3_nu, SqrtTiny ) 

    Chi_H ( iV )  &
      =  Chi_H ( iV )  +  Chi_H_S_N ( iV )  +  Chi_H_S_A ( iV )

  end procedure Compute_S_B_S_Kernel


  module procedure Compute_S_EP_E_EB_A_Kernel

    !-- Compute_Scattering_ElectronPositron_Electron_ElectronBar_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor_Xi, Factor_Chi_J, Factor_Chi_H, c1, c2, &
      Fermi_3_eM, Fermi_4_eM, &
      Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Chi_H )

    Factor_Xi     =  G_F_2  /  ( 6. * Pi_5 )
    Factor_Chi_J  =  G_F_2  /  ( 3. * Pi_3 )
    Factor_Chi_H  =  G_F_2  /  ( 6. * Pi_3 )

    c1  =  ( 1  +  2. * S_2_T_W ) ** 2
    c2  =  4. * S_2_T_W ** 2

    if ( UseDevice ) then
    else
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        Fermi_3_eM  =  Fermi_3 ( + Sign * Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eM  =  Fermi_4 ( + Sign * Mu_e ( iV ) / T ( iV ) )
        Fermi_3_eP  =  Fermi_3 ( - Sign * Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eP  =  Fermi_4 ( - Sign * Mu_e ( iV ) / T ( iV ) )

        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )

        Xi_J_S_EP ( iV )  &
          =  Factor_Xi  *  T ( iV ) ** 5  *  T_nu ( iV ) ** 4  *  Fermi_3_nu  &
             *  (     ( c1 + c2 / 6. ) * Fermi_4_eM  &
                   +  ( c2 + c1 / 6. ) * Fermi_4_eP  )

        ! if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

        !   !-- detailed balance

        !   Chi_J_S_EP ( iV )  =  Xi_J_S_EP ( iV )  /  J_Eq ( iV )

        ! else

          !-- direct 

          Chi_J_S_EP ( iV )  &  
            =  Factor_Chi_J  *  T ( iV ) ** 4  *  T_nu ( iV )  &
               *  Fermi_4_nu / max ( Fermi_3_nu, SqrtTiny )  &
               *  (     ( c1 + c2 / 6. ) * Fermi_3_eM  &
                     +  ( c2 + c1 / 6. ) * Fermi_3_eP  )

        ! end if 

        Chi_H_S_EP ( iV )  &  
          =  Chi_J_S_EP  ( iV )  &
             +  Factor_Chi_H  *  T ( iV ) ** 5  &
                *  (     ( c1 + c2 / 6. ) * Fermi_4_eM  &
                      +  ( c2 + c1 / 6. ) * Fermi_4_eP  )

        if ( M ( iV )  *  N ( iV )  <  Rho_DB ) then 
           Xi_J ( iV )  =   Xi_J ( iV )  +   Xi_J_S_EP ( iV )
          Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_S_EP ( iV )
        end if

        Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_S_EP ( iV )


      end do
    end if

  end procedure Compute_S_EP_E_EB_A_Kernel


  module procedure Compute_S_EP_E_EB_S_Kernel

    !-- Compute_Scattering_ElectronPositron_Electron_ElectronBar_Single_Kernel

    real ( KDR ) :: &
      Factor_Xi, Factor_Chi_J, Factor_Chi_H, c1, c2, &
      Fermi_3_eM, Fermi_4_eM, &
      Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu

    Factor_Xi     =  G_F_2  /  ( 6. * Pi_5 )
    Factor_Chi_J  =  G_F_2  /  ( 3. * Pi_3 )
    Factor_Chi_H  =  G_F_2  /  ( 6. * Pi_3 )

    c1  =  ( 1  +  2. * S_2_T_W ) ** 2
    c2  =  4. * S_2_T_W ** 2

    if ( T ( iV ) == 0.0_KDR ) &
      return

    Fermi_3_eM  =  Fermi_3 ( + Sign * Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eM  =  Fermi_4 ( + Sign * Mu_e ( iV ) / T ( iV ) )
    Fermi_3_eP  =  Fermi_3 ( - Sign * Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eP  =  Fermi_4 ( - Sign * Mu_e ( iV ) / T ( iV ) )

    Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
    Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )

    Xi_J_S_EP ( iV )  &
      =  Factor_Xi  *  T ( iV ) ** 5  *  T_nu ( iV ) ** 4  *  Fermi_3_nu  &
         *  (     ( c1 + c2 / 6. ) * Fermi_4_eM  &
               +  ( c2 + c1 / 6. ) * Fermi_4_eP  )

    ! if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

    !   !-- detailed balance

    !   Chi_J_S_EP ( iV )  =  Xi_J_S_EP ( iV )  /  J_Eq ( iV )

    ! else

      !-- direct 

      Chi_J_S_EP ( iV )  &  
        =  Factor_Chi_J  *  T ( iV ) ** 4  *  T_nu ( iV )  &
           *  Fermi_4_nu / max ( Fermi_3_nu, SqrtTiny )  &
           *  (     ( c1 + c2 / 6. ) * Fermi_3_eM  &
                 +  ( c2 + c1 / 6. ) * Fermi_3_eP  )

    ! end if 

    Chi_H_S_EP ( iV )  &  
      =  Chi_J_S_EP  ( iV )  &
         +  Factor_Chi_H  *  T ( iV ) ** 5  &
            *  (     ( c1 + c2 / 6. ) * Fermi_4_eM  &
                  +  ( c2 + c1 / 6. ) * Fermi_4_eP  )

    if ( M ( iV )  *  N ( iV )  <  Rho_DB ) then 
       Xi_J ( iV )  =   Xi_J ( iV )  +   Xi_J_S_EP ( iV )
      Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_S_EP ( iV )
    end if

    Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_S_EP ( iV )

  end procedure Compute_S_EP_E_EB_S_Kernel


  module procedure Compute_S_EP_HL_A_Kernel

    !-- Compute_Scattering_ElectronPositron_HeavyLepton_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Factor_Xi, Factor_Chi_J, Factor_Chi_H, c1, c2, cA, &
      Fermi_3_eM, Fermi_4_eM, &
      Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Chi_H )

    Factor_Xi     =  7. * G_F_2  /  ( 36. * Pi_5 )
    Factor_Chi_J  =  7. * G_F_2  /  ( 18. * Pi_3 )
    Factor_Chi_H  =  7. * G_F_2  /  ( 36. * Pi_3 )

    c1  =  ( 1  -  2. * S_2_T_W ) ** 2
    c2  =  4. * S_2_T_W ** 2
    cA  =  0.5_KDR * ( c1 + c2 )

    if ( UseDevice ) then
    else
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        Fermi_3_eM  =  Fermi_3 (   Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eM  =  Fermi_4 (   Mu_e ( iV ) / T ( iV ) )
        Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
        Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )

        Xi_J_S_EP ( iV )  &
          =  nSpecies  *  Factor_Xi  *  T ( iV ) ** 5  *  T_nu ( iV ) ** 4  &
             *  Fermi_3_nu  *  cA  *  ( Fermi_4_eM  +  Fermi_4_eP  )

        ! if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

        !   !-- detailed balance

        !   Chi_J_S_EP ( iV )  =  Xi_J_S_EP ( iV )  /  J_Eq ( iV )

        ! else

          !-- direct 

          Chi_J_S_EP ( iV )  &  
            =  Factor_Chi_J  *  T ( iV ) ** 4  *  T_nu ( iV )  &
               *  Fermi_4_nu / max ( Fermi_3_nu, SqrtTiny )  &
               *  cA  *  ( Fermi_3_eM  +  Fermi_3_eP  )

        ! end if 

        Chi_H_S_EP ( iV )  &  
          =  Chi_J_S_EP  ( iV )  &
             +  Factor_Chi_H  *  T ( iV ) ** 5  &
                *  cA  *  ( Fermi_4_eM  +  Fermi_4_eP )

        if ( M ( iV )  *  N ( iV )  <  Rho_DB ) then 
           Xi_J ( iV )  =   Xi_J ( iV )  +   Xi_J_S_EP ( iV )
          Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_S_EP ( iV )
        end if

        Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_S_EP ( iV )

      end do
    end if

  end procedure Compute_S_EP_HL_A_Kernel


  module procedure Compute_S_EP_HL_S_Kernel

    !-- Compute_Scattering_ElectronPositron_HeavyLepton_Single_Kernel

    real ( KDR ) :: &
      Factor_Xi, Factor_Chi_J, Factor_Chi_H, c1, c2, cA, &
      Fermi_3_eM, Fermi_4_eM, &
      Fermi_3_eP, Fermi_4_eP, &
      Fermi_3_nu, Fermi_4_nu

    Factor_Xi     =  7. * G_F_2  /  ( 36. * Pi_5 )
    Factor_Chi_J  =  7. * G_F_2  /  ( 18. * Pi_3 )
    Factor_Chi_H  =  7. * G_F_2  /  ( 36. * Pi_3 )

    c1  =  ( 1  -  2. * S_2_T_W ) ** 2
    c2  =  4. * S_2_T_W ** 2
    cA  =  0.5_KDR * ( c1 + c2 )

    if ( T ( iV ) == 0.0_KDR ) &
      return

    Fermi_3_eM  =  Fermi_3 (   Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eM  =  Fermi_4 (   Mu_e ( iV ) / T ( iV ) )
    Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
    Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

    Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
    Fermi_4_nu  =  Fermi_4 ( Eta_nu ( iV ) )

    Xi_J_S_EP ( iV )  &
      =  nSpecies  *  Factor_Xi  *  T ( iV ) ** 5  *  T_nu ( iV ) ** 4  &
         *  Fermi_3_nu  *  cA  *  ( Fermi_4_eM  +  Fermi_4_eP  )

    ! if ( M ( iV )  *  N ( iV )  >  Rho_DB ) then 

    !   !-- detailed balance

    !   Chi_J_S_EP ( iV )  =  Xi_J_S_EP ( iV )  /  J_Eq ( iV )

    ! else

      !-- direct 

      Chi_J_S_EP ( iV )  &  
        =  Factor_Chi_J  *  T ( iV ) ** 4  *  T_nu ( iV )  &
           *  Fermi_4_nu / max ( Fermi_3_nu, SqrtTiny )  &
           *  cA  *  ( Fermi_3_eM  +  Fermi_3_eP  )

    ! end if 

    Chi_H_S_EP ( iV )  &  
      =  Chi_J_S_EP  ( iV )  &
         +  Factor_Chi_H  *  T ( iV ) ** 5  &
            *  cA  *  ( Fermi_4_eM  +  Fermi_4_eP )

    if ( M ( iV )  *  N ( iV )  <  Rho_DB ) then 
       Xi_J ( iV )  =   Xi_J ( iV )  +   Xi_J_S_EP ( iV )
      Chi_J ( iV )  =  Chi_J ( iV )  +  Chi_J_S_EP ( iV )
    end if

    Chi_H ( iV )  =  Chi_H ( iV )  +  Chi_H_S_EP ( iV )

  end procedure Compute_S_EP_HL_S_Kernel


  function Fermi_2 ( Eta ) result ( F_2 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_2
      
    !$OMP_DECLARE_TARGET

    if ( Eta  >  0.0_KDR ) then
      F_2  =  Eta**3 / 3.  +  4. * Eta  +  2. * exp ( -Eta )
    else
      F_2  =  2. * exp ( Eta )
    end if

  end function Fermi_2


  function Fermi_3 ( Eta ) result ( F_3 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_3
    
    !$OMP_DECLARE_TARGET
    
    if ( Eta  >  0.0_KDR ) then
      F_3  =  Eta**4 / 4.  +  Pi_2 * Eta**2 / 2.  +  12.  -  6. * exp ( -Eta )
    else
      F_3  =  6. * exp ( Eta )
    end if

  end function Fermi_3


  function Fermi_4 ( Eta ) result ( F_4 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_4

    !$OMP_DECLARE_TARGET
    
    if ( Eta  >  0.0_KDR ) then
      F_4  =  Eta**5 / 5.  +  2. * Pi_2 * Eta**3 / 3.  +  48. * Eta  &
              +  24. * exp ( -Eta )
    else
      F_4  =  24. * exp ( Eta )
    end if

  end function Fermi_4


  function Fermi_5 ( Eta ) result ( F_5 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_5
      
    !$OMP_DECLARE_TARGET

    if ( Eta  >  0.0_KDR ) then
      F_5  =  Eta**6 / 6.  +  5. * Pi_2 * Eta**4 / 6.  +  110. * Eta**2  &
              +  240.  -  120. * exp ( -Eta )
    else
      F_5  =  120. * exp ( Eta )
    end if

  end function Fermi_5


end submodule Interactions_NM_G__Kernel
