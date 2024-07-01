#include "Preprocessor"

submodule ( Step_RK_NM_G_1D_C__Form ) Step_RK_NM_G_1D_C__Kernel

  use Basics

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
      J_Eq_EB_0, N_Eq_EB_0
    real ( KDR ) :: &
      J_Eq_E_P,  N_Eq_E_P,  &  !-- previous iteration
      J_Eq_EB_P, N_Eq_EB_P
    real ( KDR ) :: &
      E_E_P,  E_E_N,  D_E_P,  D_E_N, &  !-- previous, new
      E_EB_P, E_EB_N, D_EB_P, D_EB_N
    real ( KDR ) :: &
      dOmega, &
      SqrtTiny

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    ! dOmega  =  1.0_KDR  /  mRI

    nV  =  size ( ProperCell )

    !$OMP parallel do &
    !$OMP schedule ( OMP_SCHEDULE_HOST ) &
    !$OMP shared ( SqrtTiny ) &
    !$OMP private ( iR, iI ) &
    !$OMP private ( J_Eq_E_0, N_Eq_E_0, J_Eq_EB_0, N_Eq_EB_0 ) &
    !$OMP private ( J_Eq_E_P, N_Eq_E_P, J_Eq_EB_P, N_Eq_EB_P ) &
    !$OMP private ( E_E_P,  E_E_N,  D_E_P,  D_E_N ) &
    !$OMP private ( E_EB_P, E_EB_N, D_EB_P, D_EB_N )
    do iV = 1, nV
      if ( ProperCell ( iV ) ) then      

!call Show ( iV, '>>> iV' )
        !-- Iterate radiation and fluid energy and number to convergence

        J_Eq_E_0   =  J_Eq_E  ( iV )
        N_Eq_E_0   =  N_Eq_E  ( iV )

        J_Eq_EB_0  =  J_Eq_EB ( iV )
        N_Eq_EB_0  =  N_Eq_EB ( iV )

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

!call Show ( iI, '>>> iI' )
            !-- Compute interactions

            call I_E  % Compute ( iC, iV )
            call I_EB % Compute ( iC, iV )

            !-- For vanishing radiation initial conditions

            if ( J_Eq_E_0  ==  0.0_KDR )  &
              J_Eq_E_0  =  J_Eq_E ( iV )
            if ( N_Eq_E_0  ==  0.0_KDR )  &
              N_Eq_E_0  =  N_Eq_E ( iV )

            if ( J_Eq_EB_0  ==  0.0_KDR )  &
              J_Eq_EB_0  =  J_Eq_EB ( iV )
            if ( N_Eq_EB_0  ==  0.0_KDR )  &
              N_Eq_EB_0  =  N_Eq_EB ( iV )

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

            !-- Previous radiation values

            E_E_P  =  E_E ( iV )
            D_E_P  =  D_E ( iV )

            E_EB_P  =  E_EB ( iV )
            D_EB_P  =  D_EB ( iV )

            !-- New radiation values

            E_E_N  =  E_E_0 ( iV )  +  dT * AA * KK_E_E ( iV )  
            D_E_N  =  D_E_0 ( iV )  +  dT * AA * KK_E_D ( iV )

            E_EB_N  =  E_EB_0 ( iV )  +  dT * AA * KK_EB_E ( iV )  
            D_EB_N  =  D_EB_0 ( iV )  +  dT * AA * KK_EB_D ( iV )

            !-- New radiation values with underrelaxation

            E_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_E_P  &
                                      +  Omega ( iV )    *  E_E_N
            D_E ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_E_P  &
                                      +  Omega ( iV )    *  D_E_N

            E_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  E_EB_P  &
                                       +  Omega ( iV )    *  E_EB_N
            D_EB ( iV )  =  ( 1.0_KDR  -  Omega ( iV ) )  *  D_EB_P  &
                                       +  Omega ( iV )    *  D_EB_N

            !-- Adjusted radiation updates

            KK_E_E ( iV )  =  ( E_E ( iV )  -  E_E_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_E_D ( iV )  =  ( D_E ( iV )  -  D_E_0 ( iV ) )  &
                              /  ( dT * AA )

            KK_EB_E ( iV )  =  ( E_EB ( iV )  -  E_EB_0 ( iV ) )  &
                              /  ( dT * AA )
            KK_EB_D ( iV )  =  ( D_EB ( iV )  -  D_EB_0 ( iV ) )  &
                              /  ( dT * AA )

            !-- Fluid updates

            KK_F_E ( iV )  =  - KK_E_E ( iV )  -  KK_EB_E ( iV )
            KK_F_D ( iV )  =  - KK_E_D ( iV )  +  KK_EB_D ( iV )

            !-- New fluid values

            E_F ( iV )  =  E_F_0 ( iV )  +  dT * AA * KK_F_E ( iV )
            D_F ( iV )  =  D_F_0 ( iV )  +  dT * AA * KK_F_D ( iV )

            !-- If negative radiation density, abort for this value of Omega

            if (      E_E  ( iV )  <  0.0_KDR .or. D_E  ( iV )  <  0.0_KDR  &
                 .or. E_EB ( iV )  <  0.0_KDR .or. D_EB ( iV )  <  0.0_KDR )  &
             then
              exit Implicit
            end if

            !-- Exit test

            if ( iI  >  1 ) then

              associate &
                ( dJ_Eq_E   =>  Res_J_Eq_E  ( iI ), &
                  dN_Eq_E   =>  Res_N_Eq_E  ( iI ), &
                  dJ_Eq_EB  =>  Res_J_Eq_EB ( iI ), &
                  dN_Eq_EB  =>  Res_N_Eq_EB ( iI ) )

              dJ_Eq_E  =  abs ( J_Eq_E ( iV )  -  J_Eq_E_P )  &
                          /  max ( abs ( J_Eq_E_0 ), SqrtTiny )
              dN_Eq_E  =  abs ( N_Eq_E ( iV )  -  N_Eq_E_P )  &
                          /  max ( abs ( N_Eq_E_0 ), SqrtTiny )

              dJ_Eq_EB  =  abs ( J_Eq_EB ( iV )  -  J_Eq_EB_P )  &
                           /  max ( abs ( J_Eq_EB_0 ), SqrtTiny )
              dN_Eq_EB  =  abs ( N_Eq_EB ( iV )  -  N_Eq_EB_P )  &
                           /  max ( abs ( N_Eq_EB_0 ), SqrtTiny )

              nIterations ( iV )  =  iI
                 Residual ( iV )  =  max ( dJ_Eq_E,  dN_Eq_E, &
                                           dJ_Eq_EB, dN_Eq_EB )

              if (       dJ_Eq_E   <  Tol .and. dN_Eq_E   <  Tol  &
                   .and. dJ_Eq_EB  <  Tol .and. dN_Eq_EB  <  Tol )  &
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

!            call R_E  % ComputeFromBalanced ( iC, iV )
!            call R_EB % ComputeFromBalanced ( iC, iV )
            call F_HN % ComputeFromBalanced ( iC, iV )

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

          KK_F_S_1 ( iV )  =  - KK_E_S_1 ( iV )  -  KK_EB_S_1 ( iV )
          KK_F_S_2 ( iV )  =  - KK_E_S_2 ( iV )  -  KK_EB_S_2 ( iV )
          KK_F_S_3 ( iV )  =  - KK_E_S_3 ( iV )  -  KK_EB_S_3 ( iV )

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

        KK_F_E   ( iV )  =  0.0_KDR
        KK_F_S_1 ( iV )  =  0.0_KDR
        KK_F_S_2 ( iV )  =  0.0_KDR
        KK_F_S_3 ( iV )  =  0.0_KDR
        KK_F_D   ( iV )  =  0.0_KDR

      end if !-- ProperCell
    end do !-- iV
    !$OMP end parallel do

  end procedure SolveKernel


end submodule Step_RK_NM_G_1D_C__Kernel
