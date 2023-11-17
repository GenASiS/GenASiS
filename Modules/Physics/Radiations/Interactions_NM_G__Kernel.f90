#include "Preprocessor"

submodule ( Interactions_NM_G__Form ) Interactions_NM_G__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_EA_E_Kernel

    !-- Compute_EmissionAbsorption_Electron_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_p, Q, Factor_A, Dlta, &
      N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi, &
      Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp, &
      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    Pi      =  CONSTANT % PI
    TwoPi   =  2.0_KDR * CONSTANT % PI
    FourPi  =  4.0_KDR * CONSTANT % PI
    G_F     =  CONSTANT % FERMI_COUPLING
    g_A     =  CONSTANT % NEUTRON_AXIAL_COUPLING
    m_n     =  CONSTANT % NEUTRON_MASS
    m_p     =  CONSTANT % PROTON_MASS
    amu     =  CONSTANT % ATOMIC_MASS_UNIT

    Factor_p  =  FourPi  /  TwoPi ** 3  &
                 *  G_F ** 2  /  Pi  *  ( 1  +  3 * g_A ** 2 )
         Q    =  m_n - m_p

    Factor_A  =  FourPi  /  TwoPi ** 3  &
                 *  ( 2.0_KDR / 7.0_KDR )  *  G_F ** 2  /  Pi  *  g_A ** 2
      Dlta    =  3.0_KDR  *  UNIT % MEGA_ELECTRON_VOLT

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu ) &
      !$OMP shared ( Factor_p, Q, Factor_A, Dlta ) &
      !$OMP private ( N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi ) &
      !$OMP private ( Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q  ) &
      !$OMP private ( Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp ) &
      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  /  ( A ( iV ) * amu )
        
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

        call DFERMI ( 2.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_2_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_3_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_4_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_5_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        !--- Energy

        !-- e- + p  ->  n + nu_e 
        Xi  =  Factor_p  *  N_p  *  T ( iV ) ** 4  &
               *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                    +  Q ** 2            *  Fermi_3_e_Q )

        !-- e- + A  ->  Ap + nu_e 
        Xi  =  Xi  &
               +  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 4  &
                  *  (    T ( iV ) ** 2      *  Fermi_5_e_Qp  &
                       +  2 * Qp * T ( iV )  *  Fermi_4_e_Qp  &
                       +  Qp ** 2            *  Fermi_3_e_Qp )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_J ( iV )  =  Xi
        Chi_J ( iV )  =  Xi / J_eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_eq ( iV )

        !-- Number

        !-- e- + p  ->  n + nu_e 
        Xi  =  Factor_p  *  N_p  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                    +  Q ** 2            *  Fermi_2_e_Q )

        !-- e- + A  ->  Ap + nu_e 
        Xi  =  Xi  &
               +  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 3  &
                  *  (    T ( iV ) ** 2      *  Fermi_4_e_Qp  &
                       +  2 * Qp * T ( iV )  *  Fermi_3_e_Qp  &
                       +  Qp ** 2            *  Fermi_2_e_Qp )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_eq ( iV )

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu ) &
      !$OMP shared ( Factor_p, Q, Factor_A, Dlta ) &
      !$OMP private ( N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi ) &
      !$OMP private ( Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q  ) &
      !$OMP private ( Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp ) &
      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  /  ( A ( iV ) * amu )
        
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

        call DFERMI ( 2.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_2_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_3_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_4_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_5_e_Qp, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        !--- Energy

        !-- e- + p  ->  n + nu_e 
        Xi  =  Factor_p  *  N_p  *  T ( iV ) ** 4  &
               *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                    +  Q ** 2            *  Fermi_3_e_Q )

        !-- e- + A  ->  Ap + nu_e 
        Xi  =  Xi  &
               +  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 4  &
                  *  (    T ( iV ) ** 2      *  Fermi_5_e_Qp  &
                       +  2 * Qp * T ( iV )  *  Fermi_4_e_Qp  &
                       +  Qp ** 2            *  Fermi_3_e_Qp )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_J ( iV )  =  Xi
        Chi_J ( iV )  =  Xi / J_eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_eq ( iV )

        !-- Number

        !-- e- + p  ->  n + nu_e 
        Xi  =  Factor_p  *  N_p  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                    +  Q ** 2            *  Fermi_2_e_Q )

        !-- e- + A  ->  Ap + nu_e 
        Xi  =  Xi  &
               +  Factor_A  *  N_A  *  N_p_Z  *  N_h_N  *  T ( iV ) ** 3  &
                  *  (    T ( iV ) ** 2      *  Fermi_4_e_Qp  &
                       +  2 * Qp * T ( iV )  *  Fermi_3_e_Qp  &
                       +  Qp ** 2            *  Fermi_2_e_Qp )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_eq ( iV )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_E_Kernel


  module procedure Compute_EA_E_Bar_Kernel

    !-- Compute_EmissionAbsorption_Electron_Bar_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_n, Q, &
      N_n, Eta_e, Xi, &
      Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e, &
      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    Pi      =  CONSTANT % PI
    TwoPi   =  2.0_KDR * CONSTANT % PI
    FourPi  =  4.0_KDR * CONSTANT % PI
    G_F     =  CONSTANT % FERMI_COUPLING
    g_A     =  CONSTANT % NEUTRON_AXIAL_COUPLING
    m_n     =  CONSTANT % NEUTRON_MASS
    m_p     =  CONSTANT % PROTON_MASS
    amu     =  CONSTANT % ATOMIC_MASS_UNIT

    Factor_n  =  FourPi  /  TwoPi ** 3  &
               *  G_F ** 2  /  Pi  *  ( 1  +  3 * g_A ** 2 )
         Q  =  m_n - m_p

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu, Factor_n, Q ) &
      !$OMP private ( N_n, Eta_e, Xi ) &
      !$OMP private ( Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e ) &
      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_n    =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        !--- Energy

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 3             *  Fermi_5_e  &
                    +  3 *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                    +  3 *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                    +  Q ** 3                    *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_J ( iV )  =  Xi
        Chi_J ( iV )  =  Xi / J_eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_eq ( iV )

        !-- Number

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_eq ( iV )

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu, Factor_n, Q ) &
      !$OMP private ( N_n, Eta_e, Xi ) &
      !$OMP private ( Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e ) &
      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_n    =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
                      fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

        !--- Energy

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 3             *  Fermi_5_e  &
                    +  3 *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                    +  3 *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                    +  Q ** 3                    *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_J ( iV )  =  Xi
        Chi_J ( iV )  =  Xi / J_eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_eq ( iV )

        !-- Number

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_eq ( iV )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_E_Bar_Kernel


end submodule Interactions_NM_G__Kernel
