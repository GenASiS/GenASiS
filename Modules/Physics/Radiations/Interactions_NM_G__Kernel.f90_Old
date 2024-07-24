#include "Preprocessor"

submodule ( Interactions_NM_G__Form ) Interactions_NM_G__Kernel

  use Basics 
  
  implicit none

!  real ( KDR ) :: &
!    Pi    =  CONSTANT % PI, &
!    Pi_2  =  CONSTANT % PI ** 2

contains


  module procedure Compute_EA_E_A_Kernel

    !-- Compute_EmissionAbsorption_Electron_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_p, Q, Factor_A, Dlta, &
      N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi, &
      Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

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
      !$OMP shared ( SqrtTiny, Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu ) &
      !$OMP shared ( Factor_p, Q, Factor_A, Dlta ) &
      !$OMP private ( N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi ) &
      !$OMP private ( Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q  ) &
      !$OMP private ( Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp )!&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        ! call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
        Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
        Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
        Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

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

        ! call DFERMI ( 2.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_2_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_3_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_4_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_5_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
        Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
        Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
        Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

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
        Chi_J ( iV )  =  Xi / J_Eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_Eq ( iV )

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
        Chi_N ( iV )  =  Xi / N_Eq ( iV )

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny, Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu ) &
      !$OMP shared ( Factor_p, Q, Factor_A, Dlta ) &
      !$OMP private ( N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi ) &
      !$OMP private ( Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q  ) &
      !$OMP private ( Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp )!&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

        Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

        ! call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
        Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
        Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
        Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

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

        ! call DFERMI ( 2.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_2_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_3_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_4_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_5_e_Qp, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
        Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
        Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
        Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

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
        Chi_J ( iV )  =  Xi / J_Eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_Eq ( iV )

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
        Chi_N ( iV )  =  Xi / N_Eq ( iV )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_E_A_Kernel


  module procedure Compute_EA_E_S_Kernel

    !-- Compute_EmissionAbsorption_Electron_Single_Kernel

    real ( KDR ) :: &
      SqrtTiny, &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_p, Q, Factor_A, Dlta, &
      N_p, N_A, N_p_Z, N_h_N, Qp, Eta_e_Q, Eta_e_Qp, Xi, &
      Fermi_2_e_Q,  Fermi_3_e_Q,  Fermi_4_e_Q,  Fermi_5_e_Q, &
      Fermi_2_e_Qp, Fermi_3_e_Qp, Fermi_4_e_Qp, Fermi_5_e_Qp!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

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

    if ( T ( iV ) == 0.0_KDR ) &
      return

    N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu

    Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

    ! call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    Fermi_2_e_Q  =  Fermi_2 ( Eta_e_Q )
    Fermi_3_e_Q  =  Fermi_3 ( Eta_e_Q )
    Fermi_4_e_Q  =  Fermi_4 ( Eta_e_Q )
    Fermi_5_e_Q  =  Fermi_5 ( Eta_e_Q )

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

    ! call DFERMI ( 2.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_2_e_Qp, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 3.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_3_e_Qp, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 4.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_4_e_Qp, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 5.0_KDR, Eta_e_Qp, 0.0_KDR, Fermi_5_e_Qp, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    Fermi_2_e_Qp  =  Fermi_2 ( Eta_e_Qp )
    Fermi_3_e_Qp  =  Fermi_3 ( Eta_e_Qp )
    Fermi_4_e_Qp  =  Fermi_4 ( Eta_e_Qp )
    Fermi_5_e_Qp  =  Fermi_5 ( Eta_e_Qp )

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
    Chi_J ( iV )  =  Xi / J_Eq ( iV )

    !--- Momentum

     Xi_H ( iV )  =  0.0_KDR
    Chi_H ( iV )  =  Xi / J_Eq ( iV )

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
    Chi_N ( iV )  =  Xi / N_Eq ( iV )

  end procedure Compute_EA_E_S_Kernel


  module procedure Compute_EA_E_Bar_A_Kernel

    !-- Compute_EmissionAbsorption_Electron_Bar_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_n, Q, &
      N_n, Eta_e, Xi, &
      Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
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
      !$OMP private ( Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e ) !&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_n    =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        ! call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e  =  Fermi_2 ( -Eta_e )
        Fermi_3_e  =  Fermi_3 ( -Eta_e )
        Fermi_4_e  =  Fermi_4 ( -Eta_e )
        Fermi_5_e  =  Fermi_5 ( -Eta_e )

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
        Chi_J ( iV )  =  Xi / J_Eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_Eq ( iV )

        !-- Number

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_Eq ( iV )

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Pi, TwoPi, FourPi, G_F, g_A, m_n, m_p, amu, Factor_n, Q ) &
      !$OMP private ( N_n, Eta_e, Xi ) &
      !$OMP private ( Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e ) !&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV ) == 0.0_KDR ) &
          cycle

        N_n    =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        Eta_e  =  Mu_e ( iV )  /  T ( iV )

        ! call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_2_e  =  Fermi_2 ( -Eta_e )
        Fermi_3_e  =  Fermi_3 ( -Eta_e )
        Fermi_4_e  =  Fermi_4 ( -Eta_e )
        Fermi_5_e  =  Fermi_5 ( -Eta_e )

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
        Chi_J ( iV )  =  Xi / J_Eq ( iV )

        !--- Momentum

         Xi_H ( iV )  =  0.0_KDR
        Chi_H ( iV )  =  Xi / J_Eq ( iV )

        !-- Number

        !-- e+ + n  ->  p + nu_e_bar 
        Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
               *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                    +  2 * Q * T ( iV )  *  Fermi_3_e  &
                    +  Q ** 2            *  Fermi_2_e )

        !-- final state blocking
        Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

         Xi_N ( iV )  =  Xi
        Chi_N ( iV )  =  Xi / N_Eq ( iV )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_EA_E_Bar_A_Kernel


  module procedure Compute_EA_E_Bar_S_Kernel

    !-- Compute_EmissionAbsorption_Electron_Bar_Single_Kernel

    real ( KDR ) :: &
      Pi, TwoPi, FourPi, &
      G_F, g_A, m_n, m_p, amu, &
      Factor_n, Q, &
      N_n, Eta_e, Xi, &
      Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
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

      if ( T ( iV ) == 0.0_KDR ) &
        return

      N_n    =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
      Eta_e  =  Mu_e ( iV )  /  T ( iV )

      ! call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
      !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      ! call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
      !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      ! call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
      !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      ! call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
      !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      Fermi_2_e  =  Fermi_2 ( -Eta_e )
      Fermi_3_e  =  Fermi_3 ( -Eta_e )
      Fermi_4_e  =  Fermi_4 ( -Eta_e )
      Fermi_5_e  =  Fermi_5 ( -Eta_e )

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
      Chi_J ( iV )  =  Xi / J_Eq ( iV )

      !--- Momentum

       Xi_H ( iV )  =  0.0_KDR
      Chi_H ( iV )  =  Xi / J_Eq ( iV )

      !-- Number

      !-- e+ + n  ->  p + nu_e_bar 
      Xi  =  Factor_n  *  N_n  *  T ( iV ) ** 3  &
             *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                  +  2 * Q * T ( iV )  *  Fermi_3_e  &
                  +  Q ** 2            *  Fermi_2_e )

      !-- final state blocking
      Xi  =  Xi  *  ( 1.0_KDR  -  F_Ave ( iV ) )

       Xi_N ( iV )  =  Xi
      Chi_N ( iV )  =  Xi / N_Eq ( iV )

  end procedure Compute_EA_E_Bar_S_Kernel


  module procedure Compute_P_A_Kernel

    !-- Compute_Pair_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Pi, G_F, Sin_2_Theta_W, Factor, RhoMin, RhoMax, &
      Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, Xi
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Xi_J )

    Pi             =  CONSTANT % PI
    G_F            =  CONSTANT % FERMI_COUPLING
    Sin_2_Theta_W  =  CONSTANT % SIN_2_WEINBERG

    Factor  =  nSpecies * G_F ** 2  /  ( 9.  *  Pi ** 5 )  &
               *  ( 1.  +  Sign * 4. * Sin_2_Theta_W  &
                        +  8. * Sin_2_Theta_W ** 2 )

    RhoMin  =  0.0_KDR  * UNIT % MASS_DENSITY_CGS
    RhoMax  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Pi, G_F, Sin_2_Theta_W, Factor, RhoMin, RhoMax ) &
      !$OMP private ( Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, Xi )
      do iV = 1, nV
        if (       M ( iV )  *  N ( iV )  >  RhoMin ) then !&
!             .and. M ( iV )  *  N ( iV )  <  RhoMax ) &
!        then

          Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
          Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
          Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
          Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

          !--- Energy

          Xi  =  Factor  *  T ( iV ) ** 9  &
                 *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                            +  Fermi_3_eP * Fermi_4_eM )

          Xi_J ( iV )  =   Xi_J ( iV )  +  Xi
          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_J ( iV )  =  Chi_J ( iV )  +  Xi / J_Eq ( iV )

          !--- Momentum

          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_H ( iV )  =  Chi_H ( iV )  +  Xi / J_Eq ( iV )

          !--- Number

          Xi  =  Factor  *  T ( iV ) ** 8  &
                 *  Fermi_3_eM * Fermi_3_eP

          Xi_N ( iV )  =   Xi_N ( iV )  +  Xi
          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_N ( iV )  =  Chi_N ( iV )  +  Xi / N_Eq ( iV )

        end if
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Pi, G_F, Sin_2_Theta_W, Factor, RhoMin, RhoMax ) &
      !$OMP private ( Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, Xi )
      do iV = 1, nV
        if (       M ( iV )  *  N ( iV )  >  RhoMin ) then !&
!             .and. M ( iV )  *  N ( iV )  <  RhoMax ) &
!        then

          Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
          Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
          Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
          Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

          !--- Energy

          Xi  =  Factor  *  T ( iV ) ** 9  &
                 *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                            +  Fermi_3_eP * Fermi_4_eM )

          Xi_J ( iV )  =   Xi_J ( iV )  +  Xi
          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_J ( iV )  =  Chi_J ( iV )  +  Xi / J_Eq ( iV )

          !--- Momentum

          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_H ( iV )  =  Chi_H ( iV )  +  Xi / J_Eq ( iV )

          !--- Number

          Xi  =  Factor  *  T ( iV ) ** 8  &
                 *  Fermi_3_eM * Fermi_3_eP

          Xi_N ( iV )  =   Xi_N ( iV )  +  Xi
          if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
            Chi_N ( iV )  =  Chi_N ( iV )  +  Xi / N_Eq ( iV )

        end if
      end do !-- iV
      !$OMP end parallel do
    end if
   
  end procedure Compute_P_A_Kernel


  module procedure Compute_P_S_Kernel

    !-- Compute_Pair_Single_Kernel

    real ( KDR ) :: &
      Pi, G_F, Sin_2_Theta_W, Factor, RhoMin, RhoMax, &
      Fermi_3_eM, Fermi_4_eM, Fermi_3_eP, Fermi_4_eP, Xi
          
    Pi             =  CONSTANT % PI
    G_F            =  CONSTANT % FERMI_COUPLING
    Sin_2_Theta_W  =  CONSTANT % SIN_2_WEINBERG

    Factor  =  nSpecies * G_F ** 2  /  ( 9.  *  Pi ** 5 )  &
               *  ( 1.  +  Sign * 4. * Sin_2_Theta_W  &
                        +  8. * Sin_2_Theta_W ** 2 )

    RhoMin  =  0.0_KDR  * UNIT % MASS_DENSITY_CGS
    RhoMax  =  1.0e12_KDR * UNIT % MASS_DENSITY_CGS

   if (       M ( iV )  *  N ( iV )  >  RhoMin ) then!&
!         .and. M ( iV )  *  N ( iV )  <  RhoMax ) &
!    then

      Fermi_3_eM  =  Fermi_3 ( + Mu_e ( iV ) / T ( iV ) )
      Fermi_4_eM  =  Fermi_4 ( + Mu_e ( iV ) / T ( iV ) )
      Fermi_3_eP  =  Fermi_3 ( - Mu_e ( iV ) / T ( iV ) )
      Fermi_4_eP  =  Fermi_4 ( - Mu_e ( iV ) / T ( iV ) )

      !--- Energy

      Xi  =  Factor  *  T ( iV ) ** 9  &
             *  0.5 * (    Fermi_3_eM * Fermi_4_eP  &
                        +  Fermi_3_eP * Fermi_4_eM )

      Xi_J ( iV )  =   Xi_J ( iV )  +  Xi
      if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
        Chi_J ( iV )  =  Chi_J ( iV )  +  Xi / J_Eq ( iV )

      !--- Momentum

      if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
        Chi_H ( iV )  =  Chi_H ( iV )  +  Xi / J_Eq ( iV )

      !--- Number

      Xi  =  Factor  *  T ( iV ) ** 8  &
             *  Fermi_3_eM * Fermi_3_eP

      Xi_N ( iV )  =   Xi_N ( iV )  +  Xi
      if ( M ( iV )  *  N ( iV )  >  RhoMax ) &
        Chi_N ( iV )  =  Chi_N ( iV )  +  Xi / N_Eq ( iV )

    end if
   
  end procedure Compute_P_S_Kernel


  module procedure Compute_S_N_A_A_Kernel

    !-- Compute_Scattering_Nucleons_Nuclei_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      Pi, &
      G_F, Sin_2_Theta_W, g_A, amu, &
      Factor_n, Factor_p, Factor_A, &
      N_p, N_n, N_A, &
      Fermi_3_nu, Fermi_5_nu!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( Chi_H )

    SqrtTiny  =  tiny ( 0.0_KDR )

    Pi             =  CONSTANT % PI
    G_F            =  CONSTANT % FERMI_COUPLING
    Sin_2_Theta_W  =  CONSTANT % SIN_2_WEINBERG
    g_A            =  CONSTANT % NEUTRON_AXIAL_COUPLING
    amu            =  CONSTANT % ATOMIC_MASS_UNIT

    Factor_p  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                 *  ( ( 1. / 2.  -  2. * Sin_2_Theta_W ) ** 2  &
                      +  5. / 4. * g_A ** 2 )

    Factor_n  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                 *  ( 1. / 4.  +  5. / 4. * g_A ** 2 )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny, Pi, G_F, Sin_2_Theta_W, g_A, amu ) &
      !$OMP shared ( Factor_n, Factor_p ) &
      !$OMP private ( Factor_A, N_p, N_n, N_A ) &
      !$OMP private ( Fermi_3_nu, Fermi_5_nu ) !&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        Factor_A  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                     *  (    A ( iV )  *  ( 1. / 2.  -  2. * Sin_2_Theta_W ) &
                          -  Z ( iV )  *  ( 1.  -  2. * Sin_2_Theta_W ) ) ** 2  

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )

        ! call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_5_nu, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

        Chi_H ( iV )  &
          =  Chi_H ( iV )  &
             +  ( Factor_p * N_p  +  Factor_n * N_n  +  Factor_A * N_A )  &
                *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny, Pi, G_F, Sin_2_Theta_W, g_A, amu ) &
      !$OMP shared ( Factor_n, Factor_p ) &
      !$OMP private ( Factor_A, N_p, N_n, N_A ) &
      !$OMP private ( Fermi_3_nu, Fermi_5_nu ) !&
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        Factor_A  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                     *  (    A ( iV )  *  ( 1. / 2.  -  2. * Sin_2_Theta_W ) &
                          -  Z ( iV )  *  ( 1.  -  2. * Sin_2_Theta_W ) ) ** 2  

        N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
        N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
        N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
                /  ( max ( A ( iV ), SqrtTiny ) * amu )

        ! call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 5.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_5_nu, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
        Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

        Chi_H ( iV )  &
          =  Chi_H ( iV )  &
             +  ( Factor_p * N_p  +  Factor_n * N_n  +  Factor_A * N_A )  &
                *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_S_N_A_A_Kernel


  module procedure Compute_S_N_A_S_Kernel

    !-- Compute_Scattering_Nucleons_Nuclei_Single_Kernel

    real ( KDR ) :: &
      SqrtTiny, &
      Pi, &
      G_F, Sin_2_Theta_W, g_A, amu, &
      Factor_n, Factor_p, Factor_A, &
      N_p, N_n, N_A, &
      Fermi_3_nu, Fermi_5_nu!, &
!      fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta

    SqrtTiny  =  tiny ( 0.0_KDR )

    Pi             =  CONSTANT % PI
    G_F            =  CONSTANT % FERMI_COUPLING
    Sin_2_Theta_W  =  CONSTANT % SIN_2_WEINBERG
    g_A            =  CONSTANT % NEUTRON_AXIAL_COUPLING
    amu            =  CONSTANT % ATOMIC_MASS_UNIT

    Factor_p  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                 *  ( ( 1. / 2.  -  2. * Sin_2_Theta_W ) ** 2  &
                      +  5. / 4. * g_A ** 2 )

    Factor_n  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                 *  ( 1. / 4.  +  5. / 4. * g_A ** 2 )

    Factor_A  =  2.  *  G_F ** 2  /  ( 3. * Pi )  &
                 *  (    A ( iV )  *  ( 1. / 2.  -  2. * Sin_2_Theta_W ) &
                      -  Z ( iV )  *  ( 1.  -  2. * Sin_2_Theta_W ) ) ** 2  

    N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  amu
    N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  amu
    N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  &
            /  ( max ( A ( iV ), SqrtTiny ) * amu )

    ! call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 5.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_5_nu, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    Fermi_3_nu  =  Fermi_3 ( Eta_nu ( iV ) )
    Fermi_5_nu  =  Fermi_5 ( Eta_nu ( iV ) )

    Chi_H ( iV )  &
      =  Chi_H ( iV )  &
         +  ( Factor_p * N_p  +  Factor_n * N_n  +  Factor_A * N_A )  &
            *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

  end procedure Compute_S_N_A_S_Kernel


  function Fermi_2 ( Eta ) result ( F_2 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_2

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
    
    associate ( Pi_2 => CONSTANT % PI ** 2 )
    
    if ( Eta  >  0.0_KDR ) then
      F_3  =  Eta**4 / 4.  +  Pi_2 * Eta**2 / 2.  +  12.  -  6. * exp ( -Eta )
    else
      F_3  =  6. * exp ( Eta )
    end if
    
    end associate

  end function Fermi_3


  function Fermi_4 ( Eta ) result ( F_4 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_4
    
    associate ( Pi_2 => CONSTANT % PI ** 2 )
    
    if ( Eta  >  0.0_KDR ) then
      F_4  =  Eta**5 / 5.  +  2. * Pi_2 * Eta**3 / 3.  +  48. * Eta  &
              +  24. * exp ( -Eta )
    else
      F_4  =  24. * exp ( Eta )
    end if
    
    end associate

  end function Fermi_4


  function Fermi_5 ( Eta ) result ( F_5 )

    real ( KDR ), intent ( in ) :: &
      Eta
    real ( KDR ) :: &
      F_5

    associate ( Pi_2 => CONSTANT % PI ** 2 )

    if ( Eta  >  0.0_KDR ) then
      F_5  =  Eta**6 / 6.  +  5. * Pi_2 * Eta**4 / 6.  +  110. * Eta**2  &
              +  240.  -  120. * exp ( -Eta )
    else
      F_5  =  120. * exp ( Eta )
    end if
  
    end associate

  end function Fermi_5


end submodule Interactions_NM_G__Kernel
