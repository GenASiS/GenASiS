#include "Preprocessor"

submodule ( Step_RK_NM_G_1D_C__Form ) Step_RK_NM_G_1D_C__Kernel

  use Basics
  use Interactions_NM_G__Form, &
        only : ComputeInteractions_E_S_Kernel, &
               ComputeInteractions_EB_S_Kernel, ComputeInteractions_HL_S_Kernel
  use Fluid_P_HN__Form, & 
        only : F_P_HN_ComputeFromBalanced => ComputeFromBalanced_S_Kernel

  implicit none
  
contains 


  module procedure SolveKernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iR, &  !-- iRelaxation
      iI, &  !-- iIteration
      nV
    real ( KDR ) :: &
      J_Eq_E_0,  N_Eq_E_0,  &  !-- upon entry
      J_Eq_EB_0, N_Eq_EB_0, &
      J_Eq_X_0,  N_Eq_X_0
    real ( KDR ) :: &
      J_Eq_E_P,  N_Eq_E_P,  &  !-- previous iteration
      J_Eq_EB_P, N_Eq_EB_P, &
      J_Eq_X_P,  N_Eq_X_P
    real ( KDR ) :: &
      E_E_P,  E_E_N,  D_E_P,  D_E_N,  &  !-- previous, new
      E_EB_P, E_EB_N, D_EB_P, D_EB_N, &
      E_X_P,  E_X_N,  D_X_P,  D_X_N
    real ( KDR ) :: &
      dOmega, &
      SqrtTiny !, &
 !     M_Ref, N_Min, E_Min, &
 !     T_Min, Y_Min, Y_Safe

integer ( KDI ) :: &
  iV_Show

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
!    M_Ref   =  F_HN % BaryonMass
!    N_Min   =  F_HN % BaryonDensityMin
!    E_Min   =  F_HN % EnergyDensityMin
!    T_Min   =  F_HN % TemperatureMin
!    Y_Min   =  F_HN % ElectronFractionMin
!    Y_Safe  =  F_HN % ElectronFractionSafe
!    
!    associate &
!      ( F_V     => F_HN % Storage ( iC ) % Value, &
!        EOS     => F_HN % EOS % Table, &
!        T_L_N   => F_HN % EOS % LogDensity, &
!        T_L_T   => F_HN % EOS % LogTemperature, &
!        T_Ye    => F_HN % EOS % ElectronFraction, &
!        E_Shift => F_HN % EOS % EnergyShift, &
!        ia_F_I  => [ F_HN % BARYON_DENSITY_C, &
!                     F_HN % TEMPERATURE, F_HN % ELECTRON_FRACTION ], &
!        ia_F_O  => F_HN % EOS % iaFluidOutput, &
!        ia_E    => F_HN % EOS % iaSelected, &
!        iSolve  => F_HN % ENERGY_DENSITY_C )

    ! dOmega  =  1.0_KDR  /  mRI

    nV  =  size ( ProperCell )

!iV_Show = 88

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP shared ( SqrtTiny, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe ) &
    !$OMP private ( iR, iI ) &
    !$OMP private ( Res_J_Eq_E,  Res_N_Eq_E, &
    !$OMP           Res_J_Eq_EB, Res_N_Eq_EB, &
    !$OMP           Res_J_Eq_X,  Res_N_Eq_X ) &
    !$OMP private ( J_Eq_E_0,  N_Eq_E_0 ) &
    !$OMP private ( J_Eq_EB_0, N_Eq_EB_0 ) &
    !$OMP private ( J_Eq_X_0,  N_Eq_X_0 ) &
    !$OMP private ( J_Eq_E_P,  N_Eq_E_P ) &
    !$OMP private ( J_Eq_EB_P, N_Eq_EB_P ) &
    !$OMP private ( J_Eq_X_P,  N_Eq_X_P ) &
    !$OMP private ( E_E_P,  E_E_N,  D_E_P,  D_E_N ) &
    !$OMP private ( E_EB_P, E_EB_N, D_EB_P, D_EB_N ) &
    !$OMP private ( E_X_P,  E_X_N,  D_X_P,  D_X_N )
    do iV = 1, nV
      if ( ProperCell ( iV ) ) then      

!if ( iV == iV_Show ) &
!  call Show ( iV, '>>> iV' )

        !-- Iterate radiation and fluid energy and number to convergence

        J_Eq_E_0   =  J_Eq_E  ( iV )
        N_Eq_E_0   =  N_Eq_E  ( iV )

        J_Eq_EB_0  =  J_Eq_EB ( iV )
        N_Eq_EB_0  =  N_Eq_EB ( iV )

        J_Eq_X_0   =  J_Eq_X  ( iV )
        N_Eq_X_0   =  N_Eq_X  ( iV )

        ! if ( Omega ( iV )  ==  0.0_KDR ) then
          Omega ( iV )  =  1.0_KDR
        ! else
        !   Omega ( iV )  =  min ( Omega ( iV )  +  dOmega,  1.0_KDR )
        ! end if

        ! iR  =  0
        ! Relaxation: do 

        !   iR  =  iR + 1

        !   if ( iR  >  1 ) &
        !     Omega ( iV )  =  Omega ( iV )  -  dOmega

        !   if ( Omega ( iV )  <  0.99 * dOmega ) then

        !     !-- Solve fails with vanishing relaxation parameter Omegax

        !     !-- Reset

        !     E_E ( iV )  =  E_E_0 ( iV )
        !     D_E ( iV )  =  D_E_0 ( iV )
        !     call R_E % ComputeFromBalanced ( iC, iV )

        !     E_EB ( iV )  =  E_EB_0 ( iV )
        !     D_EB ( iV )  =  D_EB_0 ( iV )
        !     call R_EB % ComputeFromBalanced ( iC, iV )

        !     E_F ( iV )  =  E_F_0 ( iV )
        !     D_F ( iV )  =  D_F_0 ( iV )
        !     call F_HN % ComputeFromBalanced ( iC, iV )

        !     J_Eq_E ( iV )  =  J_Eq_E_0
        !     N_Eq_E ( iV )  =  N_Eq_E_0

        !     J_Eq_EB ( iV )  =  J_Eq_E_0
        !     N_Eq_EB ( iV )  =  N_Eq_E_0

        !     !-- Nullify updates
          
        !     KK_E_E   ( iV )  =  0.0_KDR
        !     KK_E_S_1 ( iV )  =  0.0_KDR
        !     KK_E_S_2 ( iV )  =  0.0_KDR
        !     KK_E_S_3 ( iV )  =  0.0_KDR
        !     KK_E_D   ( iV )  =  0.0_KDR

        !     KK_EB_E   ( iV )  =  0.0_KDR
        !     KK_EB_S_1 ( iV )  =  0.0_KDR
        !     KK_EB_S_2 ( iV )  =  0.0_KDR
        !     KK_EB_S_3 ( iV )  =  0.0_KDR
        !     KK_EB_D   ( iV )  =  0.0_KDR

        !     KK_F_E   ( iV )  =  0.0_KDR
        !     KK_F_S_1 ( iV )  =  0.0_KDR
        !     KK_F_S_2 ( iV )  =  0.0_KDR
        !     KK_F_S_3 ( iV )  =  0.0_KDR
        !     KK_F_D   ( iV )  =  0.0_KDR

        !     !-- Abort solve

        !     Error ( iV )  =  1.0_KDR
        !     exit Relaxation

        !   end if

          iI  =  0
          Implicit: do 

            iI  =  iI + 1

!if ( iV == iV_Show ) &
!  call Show ( iI, '>>> iI' )

            !-- Compute interactions

            !-- Neutrino_E
            call ComputeInteractions_E_S_Kernel &
                   ( Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                     Xi_J_EA_N_E, Xi_J_EA_A_E, Chi_J_EA_N_E, Chi_J_EA_A_E, &
                     Xi_J_P_EP_E, Chi_J_P_EP_E, &
                     Chi_H_S_N_E, Chi_H_S_A_E, &
                     Xi_J_S_EP_E, &
                     Chi_J_S_EP_E, Chi_H_S_EP_E, &
                     J_E, N_E, J_Rd_E, N_Rd_E, &
                     J_Eq_E, N_Eq_E, T_nu_E, Eta_nu_E, T_nu_EB, Eta_nu_EB, &
                     E_Ave_E, F_Ave_E, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )

            !-- Neutrino_EB
            call ComputeInteractions_EB_S_Kernel &
                   ( Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                     Xi_J_EA_N_EB, Xi_J_EA_A_EB, Chi_J_EA_N_EB, Chi_J_EA_A_EB, &
                     Xi_J_P_EP_EB, Chi_J_P_EP_EB, &
                     Chi_H_S_N_EB, Chi_H_S_A_EB, &
                     Xi_J_S_EP_EB, &
                     Chi_J_S_EP_EB, Chi_H_S_EP_EB, &
                     J_EB, N_EB, J_Rd_EB, N_Rd_EB, &
                     J_Eq_EB, N_Eq_EB, T_nu_EB, Eta_nu_EB, T_nu_E, Eta_nu_E, &
                     E_Ave_EB, F_Ave_EB, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )            
            
            !-- Neutrino_X
            call ComputeInteractions_HL_S_Kernel &
                   ( Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                     Xi_J_EA_N_X, Xi_J_EA_A_X, Chi_J_EA_N_X, Chi_J_EA_A_X, &
                     Xi_J_P_EP_X, Chi_J_P_EP_X, &
                     Chi_H_S_N_X, Chi_H_S_A_X, &
                     Xi_J_S_EP_X, &
                     Chi_J_S_EP_X, Chi_H_S_EP_X, &
                     J_X, N_X, J_Rd_X, N_Rd_X, &
                     J_Eq_X, N_Eq_X, T_nu_X, Eta_nu_X, T_nu_X, Eta_nu_X, &
                     E_Ave_X, F_Ave_X, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )            

            !-- For vanishing radiation initial conditions

            if ( J_Eq_E_0  ==  0.0_KDR )  &
              J_Eq_E_0  =  J_Eq_E ( iV )
            if ( N_Eq_E_0  ==  0.0_KDR )  &
              N_Eq_E_0  =  N_Eq_E ( iV )

            if ( J_Eq_EB_0  ==  0.0_KDR )  &
              J_Eq_EB_0  =  J_Eq_EB ( iV )
            if ( N_Eq_EB_0  ==  0.0_KDR )  &
              N_Eq_EB_0  =  N_Eq_EB ( iV )

            if ( J_Eq_X_0  ==  0.0_KDR )  &
              J_Eq_X_0  =  J_Eq_X ( iV )
            if ( N_Eq_X_0  ==  0.0_KDR )  &
              N_Eq_X_0  =  N_Eq_X ( iV )

            !-- Compute radiation energy and number updates

            KK_E_E ( iV )  &
              =  ( Xi_J_E ( iV )  -  Chi_J_E ( iV )  *  E_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_E ( iV ) * AA * dT )
            KK_E_D ( iV )  &
              =  ( Xi_N_E ( iV )  -  Chi_N_E ( iV )  *  D_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_E ( iV ) * AA * dT )

            KK_EB_E ( iV )  &
              =  ( Xi_J_EB ( iV )  -  Chi_J_EB ( iV )  *  E_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_EB ( iV ) * AA * dT )
            KK_EB_D ( iV )  &
              =  ( Xi_N_EB ( iV )  -  Chi_N_EB ( iV )  *  D_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_EB ( iV ) * AA * dT )

            KK_X_E ( iV )  &
              =  ( Xi_J_X ( iV )  -  Chi_J_X ( iV )  *  E_X_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_X ( iV ) * AA * dT )
            KK_X_D ( iV )  &
              =  ( Xi_N_X ( iV )  -  Chi_N_X ( iV )  *  D_X_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_X ( iV ) * AA * dT )

            !-- Previous radiation values

            E_E_P  =  E_E ( iV )
            D_E_P  =  D_E ( iV )

            E_EB_P  =  E_EB ( iV )
            D_EB_P  =  D_EB ( iV )

            E_X_P  =  E_X ( iV )
            D_X_P  =  D_X ( iV )

            !-- New radiation values

            E_E_N  =  E_E_0 ( iV )  +  dT * AA * KK_E_E ( iV )  
            D_E_N  =  D_E_0 ( iV )  +  dT * AA * KK_E_D ( iV )

            E_EB_N  =  E_EB_0 ( iV )  +  dT * AA * KK_EB_E ( iV )  
            D_EB_N  =  D_EB_0 ( iV )  +  dT * AA * KK_EB_D ( iV )

            E_X_N  =  E_X_0 ( iV )  +  dT * AA * KK_X_E ( iV )  
            D_X_N  =  D_X_0 ( iV )  +  dT * AA * KK_X_D ( iV )

            !-- New radiation values with underrelaxation

            E_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_E_P  &
                                      +  Omega ( iV )    *  E_E_N
            D_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_E_P  &
                                      +  Omega ( iV )    *  D_E_N

            E_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_EB_P  &
                                       +  Omega ( iV )    *  E_EB_N
            D_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_EB_P  &
                                       +  Omega ( iV )    *  D_EB_N

            E_X ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_X_P  &
                                      +  Omega ( iV )    *  E_X_N
            D_X ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_X_P  &
                                      +  Omega ( iV )    *  D_X_N

            !-- Adjusted radiation updates

            KK_E_E ( iV )  =  ( E_E ( iV )  -  E_E_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_E_D ( iV )  =  ( D_E ( iV )  -  D_E_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_EB_E ( iV )  =  ( E_EB ( iV )  -  E_EB_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_EB_D ( iV )  =  ( D_EB ( iV )  -  D_EB_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_X_E ( iV )  =  ( E_X ( iV )  -  E_X_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_X_D ( iV )  =  ( D_X ( iV )  -  D_X_0 ( iV ) )  &
                              /  ( dT * AA )

            !-- Fluid updates

            KK_F_E ( iV )  &
              =  - KK_E_E ( iV )  -  KK_EB_E ( iV )  -  KK_X_E ( iV )

            KK_F_D ( iV )  &
              =  - KK_E_D ( iV )  +  KK_EB_D ( iV )

            !-- New fluid values

            E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
            D_F ( iV )  =  D_F_0 ( iV )  +  dT * AA * KK_F_D ( iV )

            !-- If negative radiation density, abort for this value of Omega

            if (      E_E  ( iV )  <  0.0_KDR .or. D_E  ( iV )  <  0.0_KDR  &
                 .or. E_EB ( iV )  <  0.0_KDR .or. D_EB ( iV )  <  0.0_KDR  &
                 .or. E_X  ( iV )  <  0.0_KDR .or. D_X  ( iV )  <  0.0_KDR )  &
            then
              exit Implicit
            end if

            !-- Exit test

            if ( iI  >  1 ) then

              associate &
                ( dJ_Eq_E   =>  Res_J_Eq_E  ( iI ), &
                  dN_Eq_E   =>  Res_N_Eq_E  ( iI ), &
                  dJ_Eq_EB  =>  Res_J_Eq_EB ( iI ), &
                  dN_Eq_EB  =>  Res_N_Eq_EB ( iI ), &
                  dJ_Eq_X   =>  Res_J_Eq_X  ( iI ), &
                  dN_Eq_X   =>  Res_N_Eq_X  ( iI ) )
 
              dJ_Eq_E  =  abs ( J_Eq_E ( iV )  -  J_Eq_E_P )  &
                          /  max ( abs ( J_Eq_E_0 ), SqrtTiny )
              dN_Eq_E  =  abs ( N_Eq_E ( iV )  -  N_Eq_E_P )  &
                          /  max ( abs ( N_Eq_E_0 ), SqrtTiny )

              dJ_Eq_EB  =  abs ( J_Eq_EB ( iV )  -  J_Eq_EB_P )  &
                           /  max ( abs ( J_Eq_EB_0 ), SqrtTiny )
              dN_Eq_EB  =  abs ( N_Eq_EB ( iV )  -  N_Eq_EB_P )  &
                           /  max ( abs ( N_Eq_EB_0 ), SqrtTiny )

              dJ_Eq_X  =  abs ( J_Eq_X ( iV )  -  J_Eq_X_P )  &
                          /  max ( abs ( J_Eq_X_0 ), SqrtTiny )
              dN_Eq_X  =  abs ( N_Eq_X ( iV )  -  N_Eq_X_P )  &
                          /  max ( abs ( N_Eq_X_0 ), SqrtTiny )

              nIterations ( iV )  =  iI
                 Residual ( iV )  =  max ( dJ_Eq_E,  dN_Eq_E, &
                                           dJ_Eq_EB, dN_Eq_EB, &
                                           dJ_Eq_X,  dN_Eq_X )

              if (       dJ_Eq_E   <  Tol .and. dN_Eq_E   <  Tol  &
                   .and. dJ_Eq_EB  <  Tol .and. dN_Eq_EB  <  Tol  &
                   .and. dJ_Eq_X   <  Tol .and. dN_Eq_X   <  Tol )  &
              then
                Error ( iV )  =  0.0_KDR
!                exit Relaxation
                exit Implicit
              else if ( iI  ==  mII ) then
                Error ( iV )  =  1.0_KDR
                exit Implicit
              end if
              end associate !-- dJ_Eq_E, etc.

            end if !-- iI > 1 (exit test)

            !-- Prepare for next iteration

            J_Eq_E_P  =  J_Eq_E ( iV )
            N_Eq_E_P  =  N_Eq_E ( iV )

            J_Eq_EB_P  =  J_Eq_EB ( iV )
            N_Eq_EB_P  =  N_Eq_EB ( iV )

            J_Eq_X_P  =  J_Eq_X ( iV )
            N_Eq_X_P  =  N_Eq_X ( iV )

!            call R_E  % ComputeFromBalanced ( iC, iV )
!            call R_EB % ComputeFromBalanced ( iC, iV )
            
            !call F_HN % ComputeFromBalanced ( iC, iV )
            
            call F_P_HN_ComputeFromBalanced &
                   ( F_V, M_F, N_F, V_F_1, V_F_2, V_F_3, DB_F, E_F, &
                     S_F_1, S_F_2, S_F_3, P_F, T_F, EC_F, YE_F, SS_F, D_F, &
                     Mu_n_F, Mu_p_F, Mu_n_p_F, Mu_e_F, EOS, M_UU_11, &
                     M_UU_22, M_UU_33, T_L_N, T_L_T, T_Ye, M_Ref, N_Min, &
                     E_Min, T_Min, Y_Min, Y_Safe, E_Shift, ia_F_I, ia_F_O, &
                     ia_E, iSolve, iV )
            
!            call Compute_N_V_E_YE_G_S_Kernel &
!                   ( DB_F, S_F_1, S_F_2, S_F_3, E_F, D_F, M_F, &
!                     M_UU_11, M_UU_22, M_UU_33, &
!                     N_Min, E_Min, Y_Min, Y_Safe, iV, &
!                     N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F )
!            call Apply_EOS_Prologue_S_Kernel &
!                   ( M_F, N_F, P_F, T_F, EC_F, YE_F, &
!                     M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, iV )
!                     
!            !-- call CS % EOS % ComputeFromEnergy ( )
!            call F_HN % EOS % ComputeFromEnergy &
!                   ( F_HN % Storage ( iC ), &
!                     iaFluidInput = [ F_HN % BARYON_DENSITY_C, &
!                                      F_HN % TEMPERATURE, &
!                                      F_HN % ELECTRON_FRACTION ], &
!                     iSolve = F_HN % ENERGY_DENSITY_C, iV = iV )
!            call Apply_EOS_Epilogue_S_Kernel &
!                   ( N_F, P_F, T_F, SS_F, EC_F, Mu_n_F, Mu_p_F, Mu_n_p_F, &
!                     Mu_e_F, M_F, iV )

          end do Implicit

        !   !-- Reset, try again with lower relaxation parameter Omega

        !   E_E ( iV )  =  E_E_0 ( iV )
        !   D_E ( iV )  =  D_E_0 ( iV )
        !   call R_E % ComputeFromBalanced ( iC, iV )

        !   E_EB ( iV )  =  E_EB_0 ( iV )
        !   D_EB ( iV )  =  D_EB_0 ( iV )
        !   call R_EB % ComputeFromBalanced ( iC, iV )

        !   E_F ( iV )  =  E_F_0 ( iV )
        !   D_F ( iV )  =  D_F_0 ( iV )
        !   call F_HN % ComputeFromBalanced ( iC, iV )

        !   J_Eq_E ( iV )  =  J_Eq_E_0
        !   N_Eq_E ( iV )  =  N_Eq_E_0

        !   J_Eq_EB ( iV )  =  J_Eq_E_0
        !   N_Eq_EB ( iV )  =  N_Eq_E_0

        ! end do Relaxation

        !--- Momentum update

!        if ( Error ( iV )  ==  0.0_KDR ) then

          KK_E_S_1 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )
          KK_E_S_2 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )
          KK_E_S_3 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )

          KK_EB_S_1 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )
          KK_EB_S_2 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )
          KK_EB_S_3 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )

          KK_X_S_1 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )
          KK_X_S_2 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )
          KK_X_S_3 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )

! if ( iV == iV_Show ) then
!   call Show ( Xi_H_X ( iV ), '>>> Xi_H_X ( iV )' )
!   call Show ( Chi_H_X ( iV ), '>>> Chi_H_X ( iV )' )
!   call Show ( S_X_1_0 ( iV ), '>>> S_X_1_0 ( iV )' )
!   call Show ( KK_X_S_1 ( iV ), '>>> KK_X_S_1' )
! end if

          KK_F_S_1 ( iV )  &
            =  - KK_E_S_1 ( iV )  -  KK_EB_S_1 ( iV )  -  KK_X_S_1 ( iV )
          KK_F_S_2 ( iV )  &
            =  - KK_E_S_2 ( iV )  -  KK_EB_S_2 ( iV )  -  KK_X_S_2 ( iV )
          KK_F_S_3 ( iV )  &
            =  - KK_E_S_3 ( iV )  -  KK_EB_S_3 ( iV )  -  KK_X_S_3 ( iV )

!        end if !-- Error = 0

      else !-- .not. ProperCell

        KK_E_E   ( iV )  =  0.0_KDR
        KK_E_S_1 ( iV )  =  0.0_KDR
        KK_E_S_2 ( iV )  =  0.0_KDR
        KK_E_S_3 ( iV )  =  0.0_KDR
        KK_E_D   ( iV )  =  0.0_KDR

        KK_EB_E   ( iV )  =  0.0_KDR
        KK_EB_S_1 ( iV )  =  0.0_KDR
        KK_EB_S_2 ( iV )  =  0.0_KDR
        KK_EB_S_3 ( iV )  =  0.0_KDR
        KK_EB_D   ( iV )  =  0.0_KDR

        KK_X_E   ( iV )  =  0.0_KDR
        KK_X_S_1 ( iV )  =  0.0_KDR
        KK_X_S_2 ( iV )  =  0.0_KDR
        KK_X_S_3 ( iV )  =  0.0_KDR
        KK_X_D   ( iV )  =  0.0_KDR

        KK_F_E   ( iV )  =  0.0_KDR
        KK_F_S_1 ( iV )  =  0.0_KDR
        KK_F_S_2 ( iV )  =  0.0_KDR
        KK_F_S_3 ( iV )  =  0.0_KDR
        KK_F_D   ( iV )  =  0.0_KDR

      end if !-- ProperCell
    end do !-- iV
    !$OMP end parallel do
    
!    end associate !-- FV

  end procedure SolveKernel


  module procedure SolveKernelDevice

    integer ( KDI ) :: &
      iV, &  !-- iValue
      iR, &  !-- iRelaxation
      iI, &  !-- iIteration
      nV
    real ( KDR ) :: &
      J_Eq_E_0,  N_Eq_E_0,  &  !-- upon entry
      J_Eq_EB_0, N_Eq_EB_0, &
      J_Eq_X_0,  N_Eq_X_0
    real ( KDR ) :: &
      J_Eq_E_P,  N_Eq_E_P,  &  !-- previous iteration
      J_Eq_EB_P, N_Eq_EB_P, &
      J_Eq_X_P,  N_Eq_X_P
    real ( KDR ) :: &
      E_E_P,  E_E_N,  D_E_P,  D_E_N,  &  !-- previous, new
      E_EB_P, E_EB_N, D_EB_P, D_EB_N, &
      E_X_P,  E_X_N,  D_X_P,  D_X_N
    real ( KDR ) :: &
      dOmega, &
      SqrtTiny !, &
 !     M_Ref, N_Min, E_Min, &
 !     T_Min, Y_Min, Y_Safe

integer ( KDI ) :: &
  iV_Show

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
!    M_Ref   =  F_HN % BaryonMass
!    N_Min   =  F_HN % BaryonDensityMin
!    E_Min   =  F_HN % EnergyDensityMin
!    T_Min   =  F_HN % TemperatureMin
!    Y_Min   =  F_HN % ElectronFractionMin
!    Y_Safe  =  F_HN % ElectronFractionSafe
!    
!    associate &
!      ( F_V     => F_HN % Storage ( iC ) % Value, &
!        EOS     => F_HN % EOS % Table, &
!        T_L_N   => F_HN % EOS % LogDensity, &
!        T_L_T   => F_HN % EOS % LogTemperature, &
!        T_Ye    => F_HN % EOS % ElectronFraction, &
!        E_Shift => F_HN % EOS % EnergyShift, &
!        ia_F_I  => [ F_HN % BARYON_DENSITY_C, &
!                     F_HN % TEMPERATURE, F_HN % ELECTRON_FRACTION ], &
!        ia_F_O  => F_HN % EOS % iaFluidOutput, &
!        ia_E    => F_HN % EOS % iaSelected, &
!        iSolve  => F_HN % ENERGY_DENSITY_C )

    ! dOmega  =  1.0_KDR  /  mRI

    nV  =  size ( ProperCell )

!iV_Show = 88

    !$OMP OMP_TARGET_DIRECTIVE parallel do &
    !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
    !$OMP shared ( SqrtTiny, M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe ) &
    !$OMP private ( iR, iI ) &
    !$OMP private ( Res_J_Eq_E,  Res_N_Eq_E, &
    !$OMP           Res_J_Eq_EB, Res_N_Eq_EB, &
    !$OMP           Res_J_Eq_X,  Res_N_Eq_X ) &
    !$OMP private ( J_Eq_E_0,  N_Eq_E_0 ) &
    !$OMP private ( J_Eq_EB_0, N_Eq_EB_0 ) &
    !$OMP private ( J_Eq_X_0,  N_Eq_X_0 ) &
    !$OMP private ( J_Eq_E_P,  N_Eq_E_P ) &
    !$OMP private ( J_Eq_EB_P, N_Eq_EB_P ) &
    !$OMP private ( J_Eq_X_P,  N_Eq_X_P ) &
    !$OMP private ( E_E_P,  E_E_N,  D_E_P,  D_E_N ) &
    !$OMP private ( E_EB_P, E_EB_N, D_EB_P, D_EB_N ) &
    !$OMP private ( E_X_P,  E_X_N,  D_X_P,  D_X_N )
    do iV = 1, nV
      if ( ProperCell ( iV ) ) then      

!if ( iV == iV_Show ) &
!  call Show ( iV, '>>> iV' )

        !-- Iterate radiation and fluid energy and number to convergence

        J_Eq_E_0   =  J_Eq_E  ( iV )
        N_Eq_E_0   =  N_Eq_E  ( iV )

        J_Eq_EB_0  =  J_Eq_EB ( iV )
        N_Eq_EB_0  =  N_Eq_EB ( iV )

        J_Eq_X_0   =  J_Eq_X  ( iV )
        N_Eq_X_0   =  N_Eq_X  ( iV )

        ! if ( Omega ( iV )  ==  0.0_KDR ) then
          Omega ( iV )  =  1.0_KDR
        ! else
        !   Omega ( iV )  =  min ( Omega ( iV )  +  dOmega,  1.0_KDR )
        ! end if

        ! iR  =  0
        ! Relaxation: do 

        !   iR  =  iR + 1

        !   if ( iR  >  1 ) &
        !     Omega ( iV )  =  Omega ( iV )  -  dOmega

        !   if ( Omega ( iV )  <  0.99 * dOmega ) then

        !     !-- Solve fails with vanishing relaxation parameter Omegax

        !     !-- Reset

        !     E_E ( iV )  =  E_E_0 ( iV )
        !     D_E ( iV )  =  D_E_0 ( iV )
        !     call R_E % ComputeFromBalanced ( iC, iV )

        !     E_EB ( iV )  =  E_EB_0 ( iV )
        !     D_EB ( iV )  =  D_EB_0 ( iV )
        !     call R_EB % ComputeFromBalanced ( iC, iV )

        !     E_F ( iV )  =  E_F_0 ( iV )
        !     D_F ( iV )  =  D_F_0 ( iV )
        !     call F_HN % ComputeFromBalanced ( iC, iV )

        !     J_Eq_E ( iV )  =  J_Eq_E_0
        !     N_Eq_E ( iV )  =  N_Eq_E_0

        !     J_Eq_EB ( iV )  =  J_Eq_E_0
        !     N_Eq_EB ( iV )  =  N_Eq_E_0

        !     !-- Nullify updates
          
        !     KK_E_E   ( iV )  =  0.0_KDR
        !     KK_E_S_1 ( iV )  =  0.0_KDR
        !     KK_E_S_2 ( iV )  =  0.0_KDR
        !     KK_E_S_3 ( iV )  =  0.0_KDR
        !     KK_E_D   ( iV )  =  0.0_KDR

        !     KK_EB_E   ( iV )  =  0.0_KDR
        !     KK_EB_S_1 ( iV )  =  0.0_KDR
        !     KK_EB_S_2 ( iV )  =  0.0_KDR
        !     KK_EB_S_3 ( iV )  =  0.0_KDR
        !     KK_EB_D   ( iV )  =  0.0_KDR

        !     KK_F_E   ( iV )  =  0.0_KDR
        !     KK_F_S_1 ( iV )  =  0.0_KDR
        !     KK_F_S_2 ( iV )  =  0.0_KDR
        !     KK_F_S_3 ( iV )  =  0.0_KDR
        !     KK_F_D   ( iV )  =  0.0_KDR

        !     !-- Abort solve

        !     Error ( iV )  =  1.0_KDR
        !     exit Relaxation

        !   end if

          iI  =  0
          Implicit: do 

            iI  =  iI + 1

!if ( iV == iV_Show ) &
!  call Show ( iI, '>>> iI' )

            !-- Compute interactions

            !-- Neutrino_E
            call ComputeInteractions_E_S_Kernel &
                   ( Xi_J_E, Xi_H_E, Xi_N_E, Chi_J_E, Chi_H_E, Chi_N_E, &
                     Xi_J_EA_N_E, Xi_J_EA_A_E, Chi_J_EA_N_E, Chi_J_EA_A_E, &
                     Xi_J_P_EP_E, Chi_J_P_EP_E, &
                     Chi_H_S_N_E, Chi_H_S_A_E, &
                     Xi_J_S_EP_E, &
                     Chi_J_S_EP_E, Chi_H_S_EP_E, &
                     J_E, N_E, J_Rd_E, N_Rd_E, &
                     J_Eq_E, N_Eq_E, T_nu_E, Eta_nu_E, T_nu_EB, Eta_nu_EB, &
                     E_Ave_E, F_Ave_E, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )

            !-- Neutrino_EB
            call ComputeInteractions_EB_S_Kernel &
                   ( Xi_J_EB, Xi_H_EB, Xi_N_EB, Chi_J_EB, Chi_H_EB, Chi_N_EB, &
                     Xi_J_EA_N_EB, Xi_J_EA_A_EB, Chi_J_EA_N_EB, Chi_J_EA_A_EB, &
                     Xi_J_P_EP_EB, Chi_J_P_EP_EB, &
                     Chi_H_S_N_EB, Chi_H_S_A_EB, &
                     Xi_J_S_EP_EB, &
                     Chi_J_S_EP_EB, Chi_H_S_EP_EB, &
                     J_EB, N_EB, J_Rd_EB, N_Rd_EB, &
                     J_Eq_EB, N_Eq_EB, T_nu_EB, Eta_nu_EB, T_nu_E, Eta_nu_E, &
                     E_Ave_EB, F_Ave_EB, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )            
            
            !-- Neutrino_X
            call ComputeInteractions_HL_S_Kernel &
                   ( Xi_J_X, Xi_H_X, Xi_N_X, Chi_J_X, Chi_H_X, Chi_N_X, &
                     Xi_J_EA_N_X, Xi_J_EA_A_X, Chi_J_EA_N_X, Chi_J_EA_A_X, &
                     Xi_J_P_EP_X, Chi_J_P_EP_X, &
                     Chi_H_S_N_X, Chi_H_S_A_X, &
                     Xi_J_S_EP_X, &
                     Chi_J_S_EP_X, Chi_H_S_EP_X, &
                     J_X, N_X, J_Rd_X, N_Rd_X, &
                     J_Eq_X, N_Eq_X, T_nu_X, Eta_nu_X, T_nu_X, Eta_nu_X, &
                     E_Ave_X, F_Ave_X, &
                     M_F, N_F, T_F, X_n_F, X_p_F, X_A_F, Z_F, A_F, &
                     Mu_e_F, Mu_n_p_F, &
                     Rho_DB, iV )            

            !-- For vanishing radiation initial conditions

            if ( J_Eq_E_0  ==  0.0_KDR )  &
              J_Eq_E_0  =  J_Eq_E ( iV )
            if ( N_Eq_E_0  ==  0.0_KDR )  &
              N_Eq_E_0  =  N_Eq_E ( iV )

            if ( J_Eq_EB_0  ==  0.0_KDR )  &
              J_Eq_EB_0  =  J_Eq_EB ( iV )
            if ( N_Eq_EB_0  ==  0.0_KDR )  &
              N_Eq_EB_0  =  N_Eq_EB ( iV )

            if ( J_Eq_X_0  ==  0.0_KDR )  &
              J_Eq_X_0  =  J_Eq_X ( iV )
            if ( N_Eq_X_0  ==  0.0_KDR )  &
              N_Eq_X_0  =  N_Eq_X ( iV )

            !-- Compute radiation energy and number updates

            KK_E_E ( iV )  &
              =  ( Xi_J_E ( iV )  -  Chi_J_E ( iV )  *  E_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_E ( iV ) * AA * dT )
            KK_E_D ( iV )  &
              =  ( Xi_N_E ( iV )  -  Chi_N_E ( iV )  *  D_E_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_E ( iV ) * AA * dT )

            KK_EB_E ( iV )  &
              =  ( Xi_J_EB ( iV )  -  Chi_J_EB ( iV )  *  E_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_EB ( iV ) * AA * dT )
            KK_EB_D ( iV )  &
              =  ( Xi_N_EB ( iV )  -  Chi_N_EB ( iV )  *  D_EB_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_EB ( iV ) * AA * dT )

            KK_X_E ( iV )  &
              =  ( Xi_J_X ( iV )  -  Chi_J_X ( iV )  *  E_X_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_J_X ( iV ) * AA * dT )
            KK_X_D ( iV )  &
              =  ( Xi_N_X ( iV )  -  Chi_N_X ( iV )  *  D_X_0 ( iV ) ) &
                 /  ( 1.0_KDR  +  Chi_N_X ( iV ) * AA * dT )

            !-- Previous radiation values

            E_E_P  =  E_E ( iV )
            D_E_P  =  D_E ( iV )

            E_EB_P  =  E_EB ( iV )
            D_EB_P  =  D_EB ( iV )

            E_X_P  =  E_X ( iV )
            D_X_P  =  D_X ( iV )

            !-- New radiation values

            E_E_N  =  E_E_0 ( iV )  +  dT * AA * KK_E_E ( iV )  
            D_E_N  =  D_E_0 ( iV )  +  dT * AA * KK_E_D ( iV )

            E_EB_N  =  E_EB_0 ( iV )  +  dT * AA * KK_EB_E ( iV )  
            D_EB_N  =  D_EB_0 ( iV )  +  dT * AA * KK_EB_D ( iV )

            E_X_N  =  E_X_0 ( iV )  +  dT * AA * KK_X_E ( iV )  
            D_X_N  =  D_X_0 ( iV )  +  dT * AA * KK_X_D ( iV )

            !-- New radiation values with underrelaxation

            E_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_E_P  &
                                      +  Omega ( iV )    *  E_E_N
            D_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_E_P  &
                                      +  Omega ( iV )    *  D_E_N

            E_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_EB_P  &
                                       +  Omega ( iV )    *  E_EB_N
            D_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_EB_P  &
                                       +  Omega ( iV )    *  D_EB_N

            E_X ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_X_P  &
                                      +  Omega ( iV )    *  E_X_N
            D_X ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_X_P  &
                                      +  Omega ( iV )    *  D_X_N

            !-- Adjusted radiation updates

            KK_E_E ( iV )  =  ( E_E ( iV )  -  E_E_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_E_D ( iV )  =  ( D_E ( iV )  -  D_E_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_EB_E ( iV )  =  ( E_EB ( iV )  -  E_EB_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_EB_D ( iV )  =  ( D_EB ( iV )  -  D_EB_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_X_E ( iV )  =  ( E_X ( iV )  -  E_X_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_X_D ( iV )  =  ( D_X ( iV )  -  D_X_0 ( iV ) )  &
                              /  ( dT * AA )

            !-- Fluid updates

            KK_F_E ( iV )  &
              =  - KK_E_E ( iV )  -  KK_EB_E ( iV )  -  KK_X_E ( iV )

            KK_F_D ( iV )  &
              =  - KK_E_D ( iV )  +  KK_EB_D ( iV )

            !-- New fluid values

            E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
            D_F ( iV )  =  D_F_0 ( iV )  +  dT * AA * KK_F_D ( iV )

            !-- If negative radiation density, abort for this value of Omega

            if (      E_E  ( iV )  <  0.0_KDR .or. D_E  ( iV )  <  0.0_KDR  &
                 .or. E_EB ( iV )  <  0.0_KDR .or. D_EB ( iV )  <  0.0_KDR  &
                 .or. E_X  ( iV )  <  0.0_KDR .or. D_X  ( iV )  <  0.0_KDR )  &
            then
              exit Implicit
            end if

            !-- Exit test

            if ( iI  >  1 ) then

              associate &
                ( dJ_Eq_E   =>  Res_J_Eq_E  ( iI ), &
                  dN_Eq_E   =>  Res_N_Eq_E  ( iI ), &
                  dJ_Eq_EB  =>  Res_J_Eq_EB ( iI ), &
                  dN_Eq_EB  =>  Res_N_Eq_EB ( iI ), &
                  dJ_Eq_X   =>  Res_J_Eq_X  ( iI ), &
                  dN_Eq_X   =>  Res_N_Eq_X  ( iI ) )
 
              dJ_Eq_E  =  abs ( J_Eq_E ( iV )  -  J_Eq_E_P )  &
                          /  max ( abs ( J_Eq_E_0 ), SqrtTiny )
              dN_Eq_E  =  abs ( N_Eq_E ( iV )  -  N_Eq_E_P )  &
                          /  max ( abs ( N_Eq_E_0 ), SqrtTiny )

              dJ_Eq_EB  =  abs ( J_Eq_EB ( iV )  -  J_Eq_EB_P )  &
                           /  max ( abs ( J_Eq_EB_0 ), SqrtTiny )
              dN_Eq_EB  =  abs ( N_Eq_EB ( iV )  -  N_Eq_EB_P )  &
                           /  max ( abs ( N_Eq_EB_0 ), SqrtTiny )

              dJ_Eq_X  =  abs ( J_Eq_X ( iV )  -  J_Eq_X_P )  &
                          /  max ( abs ( J_Eq_X_0 ), SqrtTiny )
              dN_Eq_X  =  abs ( N_Eq_X ( iV )  -  N_Eq_X_P )  &
                          /  max ( abs ( N_Eq_X_0 ), SqrtTiny )

              nIterations ( iV )  =  iI
                 Residual ( iV )  =  max ( dJ_Eq_E,  dN_Eq_E, &
                                           dJ_Eq_EB, dN_Eq_EB, &
                                           dJ_Eq_X,  dN_Eq_X )

              if (       dJ_Eq_E   <  Tol .and. dN_Eq_E   <  Tol  &
                   .and. dJ_Eq_EB  <  Tol .and. dN_Eq_EB  <  Tol  &
                   .and. dJ_Eq_X   <  Tol .and. dN_Eq_X   <  Tol )  &
              then
                Error ( iV )  =  0.0_KDR
!                exit Relaxation
                exit Implicit
              else if ( iI  ==  mII ) then
                Error ( iV )  =  1.0_KDR
                exit Implicit
              end if
              end associate !-- dJ_Eq_E, etc.

            end if !-- iI > 1 (exit test)

            !-- Prepare for next iteration

            J_Eq_E_P  =  J_Eq_E ( iV )
            N_Eq_E_P  =  N_Eq_E ( iV )

            J_Eq_EB_P  =  J_Eq_EB ( iV )
            N_Eq_EB_P  =  N_Eq_EB ( iV )

            J_Eq_X_P  =  J_Eq_X ( iV )
            N_Eq_X_P  =  N_Eq_X ( iV )

!            call R_E  % ComputeFromBalanced ( iC, iV )
!            call R_EB % ComputeFromBalanced ( iC, iV )
            
            !call F_HN % ComputeFromBalanced ( iC, iV )
            
            call F_P_HN_ComputeFromBalanced &
                   ( F_V, M_F, N_F, V_F_1, V_F_2, V_F_3, DB_F, E_F, &
                     S_F_1, S_F_2, S_F_3, P_F, T_F, EC_F, YE_F, SS_F, D_F, &
                     Mu_n_F, Mu_p_F, Mu_n_p_F, Mu_e_F, EOS, M_UU_11, &
                     M_UU_22, M_UU_33, T_L_N, T_L_T, T_Ye, M_Ref, N_Min, &
                     E_Min, T_Min, Y_Min, Y_Safe, E_Shift, ia_F_I, ia_F_O, &
                     ia_E, iSolve, iV )
            
!            call Compute_N_V_E_YE_G_S_Kernel &
!                   ( DB_F, S_F_1, S_F_2, S_F_3, E_F, D_F, M_F, &
!                     M_UU_11, M_UU_22, M_UU_33, &
!                     N_Min, E_Min, Y_Min, Y_Safe, iV, &
!                     N_F, V_F_1, V_F_2, V_F_3, EC_F, YE_F )
!            call Apply_EOS_Prologue_S_Kernel &
!                   ( M_F, N_F, P_F, T_F, EC_F, YE_F, &
!                     M_Ref, N_Min, E_Min, T_Min, Y_Min, Y_Safe, iV )
!                     
!            !-- call CS % EOS % ComputeFromEnergy ( )
!            call F_HN % EOS % ComputeFromEnergy &
!                   ( F_HN % Storage ( iC ), &
!                     iaFluidInput = [ F_HN % BARYON_DENSITY_C, &
!                                      F_HN % TEMPERATURE, &
!                                      F_HN % ELECTRON_FRACTION ], &
!                     iSolve = F_HN % ENERGY_DENSITY_C, iV = iV )
!            call Apply_EOS_Epilogue_S_Kernel &
!                   ( N_F, P_F, T_F, SS_F, EC_F, Mu_n_F, Mu_p_F, Mu_n_p_F, &
!                     Mu_e_F, M_F, iV )

          end do Implicit

        !   !-- Reset, try again with lower relaxation parameter Omega

        !   E_E ( iV )  =  E_E_0 ( iV )
        !   D_E ( iV )  =  D_E_0 ( iV )
        !   call R_E % ComputeFromBalanced ( iC, iV )

        !   E_EB ( iV )  =  E_EB_0 ( iV )
        !   D_EB ( iV )  =  D_EB_0 ( iV )
        !   call R_EB % ComputeFromBalanced ( iC, iV )

        !   E_F ( iV )  =  E_F_0 ( iV )
        !   D_F ( iV )  =  D_F_0 ( iV )
        !   call F_HN % ComputeFromBalanced ( iC, iV )

        !   J_Eq_E ( iV )  =  J_Eq_E_0
        !   N_Eq_E ( iV )  =  N_Eq_E_0

        !   J_Eq_EB ( iV )  =  J_Eq_E_0
        !   N_Eq_EB ( iV )  =  N_Eq_E_0

        ! end do Relaxation

        !--- Momentum update

!        if ( Error ( iV )  ==  0.0_KDR ) then

          KK_E_S_1 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )
          KK_E_S_2 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )
          KK_E_S_3 ( iV )  &
            =  ( Xi_H_E ( iV )  -  Chi_H_E ( iV )  *  S_E_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_E ( iV ) * AA * dT )

          KK_EB_S_1 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )
          KK_EB_S_2 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )
          KK_EB_S_3 ( iV )  &
            =  ( Xi_H_EB ( iV )  -  Chi_H_EB ( iV )  *  S_EB_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_EB ( iV ) * AA * dT )

          KK_X_S_1 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_1_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )
          KK_X_S_2 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_2_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )
          KK_X_S_3 ( iV )  &
            =  ( Xi_H_X ( iV )  -  Chi_H_X ( iV )  *  S_X_3_0 ( iV ) ) &
               /  ( 1.0_KDR  +  Chi_H_X ( iV ) * AA * dT )

! if ( iV == iV_Show ) then
!   call Show ( Xi_H_X ( iV ), '>>> Xi_H_X ( iV )' )
!   call Show ( Chi_H_X ( iV ), '>>> Chi_H_X ( iV )' )
!   call Show ( S_X_1_0 ( iV ), '>>> S_X_1_0 ( iV )' )
!   call Show ( KK_X_S_1 ( iV ), '>>> KK_X_S_1' )
! end if

          KK_F_S_1 ( iV )  &
            =  - KK_E_S_1 ( iV )  -  KK_EB_S_1 ( iV )  -  KK_X_S_1 ( iV )
          KK_F_S_2 ( iV )  &
            =  - KK_E_S_2 ( iV )  -  KK_EB_S_2 ( iV )  -  KK_X_S_2 ( iV )
          KK_F_S_3 ( iV )  &
            =  - KK_E_S_3 ( iV )  -  KK_EB_S_3 ( iV )  -  KK_X_S_3 ( iV )

!        end if !-- Error = 0

      else !-- .not. ProperCell

        KK_E_E   ( iV )  =  0.0_KDR
        KK_E_S_1 ( iV )  =  0.0_KDR
        KK_E_S_2 ( iV )  =  0.0_KDR
        KK_E_S_3 ( iV )  =  0.0_KDR
        KK_E_D   ( iV )  =  0.0_KDR

        KK_EB_E   ( iV )  =  0.0_KDR
        KK_EB_S_1 ( iV )  =  0.0_KDR
        KK_EB_S_2 ( iV )  =  0.0_KDR
        KK_EB_S_3 ( iV )  =  0.0_KDR
        KK_EB_D   ( iV )  =  0.0_KDR

        KK_X_E   ( iV )  =  0.0_KDR
        KK_X_S_1 ( iV )  =  0.0_KDR
        KK_X_S_2 ( iV )  =  0.0_KDR
        KK_X_S_3 ( iV )  =  0.0_KDR
        KK_X_D   ( iV )  =  0.0_KDR

        KK_F_E   ( iV )  =  0.0_KDR
        KK_F_S_1 ( iV )  =  0.0_KDR
        KK_F_S_2 ( iV )  =  0.0_KDR
        KK_F_S_3 ( iV )  =  0.0_KDR
        KK_F_D   ( iV )  =  0.0_KDR

      end if !-- ProperCell
    end do !-- iV
    !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
!    end associate !-- FV

  end procedure SolveKernelDevice


end submodule Step_RK_NM_G_1D_C__Kernel
