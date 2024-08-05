#include "Preprocessor"

submodule ( NeutrinoMoments_G__Form ) NeutrinoMoments_G__Kernel

  use Basics 
  
  implicit none

  real ( KDR ) :: &
    Pi    =  CONSTANT % PI, &
    Pi_2  =  CONSTANT % PI ** 2

contains


  module procedure Compute_E_S_G_G_Kernel

    !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV = size ( E )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( J ( iV )  <  SqrtTiny ) &
          J ( iV )  =  SqrtTiny
        if ( N ( iV )  <  SqrtTiny ) &
          N ( iV )  =  SqrtTiny

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        if ( H  >  J ( iV ) ) then

          H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
          H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
          H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )
  
          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        end if

        !-- Moment factors ( Minerbo SF )

        FF ( iV )  =  H  /  J ( iV )

        SF ( iV )  =  1.0_KDR / 3.0_KDR &
                      +  2.0_KDR / 3.0_KDR &
                         *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                              * ( 3.0_KDR  -  FF ( iV )  &
                                  +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                         / ( 1.0_KDR / 3.0_KDR )

        !-- FIXME: Add velocity dependence

        E ( iV )  =  J ( iV )

        S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
        S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
        S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

        G ( iV )  =  N ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( J ( iV )  <  SqrtTiny ) &
          J ( iV )  =  SqrtTiny
        if ( N ( iV )  <  SqrtTiny ) &
          N ( iV )  =  SqrtTiny

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        if ( H  >  J ( iV ) ) then

          H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
          H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
          H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )
  
          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        end if

        !-- Moment factors ( Minerbo SF )

        FF ( iV )  =  H  /  J ( iV )

        SF ( iV )  =  1.0_KDR / 3.0_KDR &
                      +  2.0_KDR / 3.0_KDR &
                         *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                              * ( 3.0_KDR  -  FF ( iV )  &
                                  +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                         / ( 1.0_KDR / 3.0_KDR )

        !-- FIXME: Add velocity dependence

        E ( iV )  =  J ( iV )

        S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
        S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
        S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

        G ( iV )  =  N ( iV )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_E_S_G_G_Kernel


  module procedure Compute_J_H_N_G_A_Kernel

    !-- Compute_ComovingEnergy_Momentum_Galileo_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( J )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( E ( iV )  >  0.0_KDR  .and.  N ( iV )  >  0.0_KDR ) then

          ! call ComputeComovingNonlinearSolve &
          !        ( J ( iV ), H_1 ( iV ), H_2 ( iV ), H_3 ( iV ), FF ( iV ), &
          !          SF ( iV ), E ( iV ), S_1 ( iV ), S_2 ( iV ), S_3 ( iV ), &
          !          M_DD_22 ( iV ), M_DD_33 ( iV ), M_UU_22 ( iV ), &
          !          M_UU_33 ( iV ), V_1 ( iV ), V_2 ( iV ), V_3 ( iV ), &
          !          Success, Delta_J_J, Delta_H_H )
          ! if ( .not. Success ) then
          !   call Show ( '>>> ComputeComoving fail', CONSOLE % ERROR )
          !   call Show ( RM % Name, '>>> Species', CONSOLE % ERROR )
          !   call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
          !               CONSOLE % ERROR )
          !   call Show ( iV, '>>> iV', CONSOLE % ERROR )
          !   call Show ( J ( iV ), '>>> J', CONSOLE % ERROR )
          !   call Show ( H_1 ( iV ), '>>> H_1', CONSOLE % ERROR )
          !   call Show ( H_2 ( iV ), '>>> H_2', CONSOLE % ERROR )
          !   call Show ( H_3 ( iV ), '>>> H_3', CONSOLE % ERROR )
          !   call Show ( Delta_J_J, '>>> Delta_J_J', CONSOLE % ERROR )
          !   call Show ( Delta_H_H, '>>> Delta_H_H', CONSOLE % ERROR )
          ! end if

          !-- FIXME: Do solve above

          J ( iV )  =  max ( E ( iV ), SqrtTiny )
          
          H_1 ( iV )  =  M_UU_11 ( iV )  *  S_1 ( iV )
          H_2 ( iV )  =  M_UU_22 ( iV )  *  S_2 ( iV )
          H_3 ( iV )  =  M_UU_33 ( iV )  *  S_3 ( iV )          

          N ( iV )  =  max ( G ( iV ), SqrtTiny )
          
          !-- Moment factors ( Minerbo SF )

          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

          FF ( iV )  =  H  /  J ( iV )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        else

          J   ( iV )  =  SqrtTiny
          H_1 ( iV )  =  0.0_KDR
          H_2 ( iV )  =  0.0_KDR
          H_3 ( iV )  =  0.0_KDR
          N   ( iV )  =  SqrtTiny
          E   ( iV )  =  SqrtTiny
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          G   ( iV )  =  SqrtTiny
          FF  ( iV )  =  0.0_KDR
          SF  ( iV )  =  1.0_KDR / 3.0_KDR

          cycle 

        end if

        SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                         / ( 1.0_KDR / 3.0_KDR )

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
        
        if ( H  >  J ( iV ) ) then

          H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
          H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
          H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )

          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

          !-- Moment factors ( Minerbo SF )

          FF ( iV )  =  H  /  J ( iV )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

          SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                           / ( 1.0_KDR / 3.0_KDR )

          !-- FIXME: Add velocity dependence

          E ( iV )  =  J ( iV )

          S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
          S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
          S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

          G ( iV )  =  N ( iV )

        end if

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( E ( iV )  >  0.0_KDR  .and.  N ( iV )  >  0.0_KDR ) then

          ! call ComputeComovingNonlinearSolve &
          !        ( J ( iV ), H_1 ( iV ), H_2 ( iV ), H_3 ( iV ), FF ( iV ), &
          !          SF ( iV ), E ( iV ), S_1 ( iV ), S_2 ( iV ), S_3 ( iV ), &
          !          M_DD_22 ( iV ), M_DD_33 ( iV ), M_UU_22 ( iV ), &
          !          M_UU_33 ( iV ), V_1 ( iV ), V_2 ( iV ), V_3 ( iV ), &
          !          Success, Delta_J_J, Delta_H_H )
          ! if ( .not. Success ) then
          !   call Show ( '>>> ComputeComoving fail', CONSOLE % ERROR )
          !   call Show ( RM % Name, '>>> Species', CONSOLE % ERROR )
          !   call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
          !               CONSOLE % ERROR )
          !   call Show ( iV, '>>> iV', CONSOLE % ERROR )
          !   call Show ( J ( iV ), '>>> J', CONSOLE % ERROR )
          !   call Show ( H_1 ( iV ), '>>> H_1', CONSOLE % ERROR )
          !   call Show ( H_2 ( iV ), '>>> H_2', CONSOLE % ERROR )
          !   call Show ( H_3 ( iV ), '>>> H_3', CONSOLE % ERROR )
          !   call Show ( Delta_J_J, '>>> Delta_J_J', CONSOLE % ERROR )
          !   call Show ( Delta_H_H, '>>> Delta_H_H', CONSOLE % ERROR )
          ! end if

          !-- FIXME: Do solve above

          J ( iV )  =  max ( E ( iV ), SqrtTiny )
          
          H_1 ( iV )  =  M_UU_11 ( iV )  *  S_1 ( iV )
          H_2 ( iV )  =  M_UU_22 ( iV )  *  S_2 ( iV )
          H_3 ( iV )  =  M_UU_33 ( iV )  *  S_3 ( iV )          

          N ( iV )  =  max ( G ( iV ), SqrtTiny )
          
          !-- Moment factors ( Minerbo SF )

          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

          FF ( iV )  =  H  /  J ( iV )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        else

          J   ( iV )  =  SqrtTiny
          H_1 ( iV )  =  0.0_KDR
          H_2 ( iV )  =  0.0_KDR
          H_3 ( iV )  =  0.0_KDR
          N   ( iV )  =  SqrtTiny
          E   ( iV )  =  SqrtTiny
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          G   ( iV )  =  SqrtTiny
          FF  ( iV )  =  0.0_KDR
          SF  ( iV )  =  1.0_KDR / 3.0_KDR

          cycle 

        end if

        SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                         / ( 1.0_KDR / 3.0_KDR )

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
        
        if ( H  >  J ( iV ) ) then

          H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
          H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
          H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )

          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

          !-- Moment factors ( Minerbo SF )

          FF ( iV )  =  H  /  J ( iV )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

          SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                           / ( 1.0_KDR / 3.0_KDR )

          !-- FIXME: Add velocity dependence

          E ( iV )  =  J ( iV )

          S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
          S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
          S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

          G ( iV )  =  N ( iV )

        end if

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_J_H_N_G_A_Kernel


  module procedure Compute_J_H_N_G_S_Kernel

    !-- Compute_ComovingEnergy_Momentum_Galileo_Single_Kernel

    real ( KDR ) :: &
      H, &
      SqrtTiny

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( E ( iV )  >  0.0_KDR  .and.  N ( iV )  >  0.0_KDR ) then

      ! call ComputeComovingNonlinearSolve &
      !        ( J ( iV ), H_1 ( iV ), H_2 ( iV ), H_3 ( iV ), FF ( iV ), &
      !          SF ( iV ), E ( iV ), S_1 ( iV ), S_2 ( iV ), S_3 ( iV ), &
      !          M_DD_22 ( iV ), M_DD_33 ( iV ), M_UU_22 ( iV ), &
      !          M_UU_33 ( iV ), V_1 ( iV ), V_2 ( iV ), V_3 ( iV ), &
      !          Success, Delta_J_J, Delta_H_H )
      ! if ( .not. Success ) then
      !   call Show ( '>>> ComputeComoving fail', CONSOLE % ERROR )
      !   call Show ( RM % Name, '>>> Species', CONSOLE % ERROR )
      !   call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', &
      !               CONSOLE % ERROR )
      !   call Show ( iV, '>>> iV', CONSOLE % ERROR )
      !   call Show ( J ( iV ), '>>> J', CONSOLE % ERROR )
      !   call Show ( H_1 ( iV ), '>>> H_1', CONSOLE % ERROR )
      !   call Show ( H_2 ( iV ), '>>> H_2', CONSOLE % ERROR )
      !   call Show ( H_3 ( iV ), '>>> H_3', CONSOLE % ERROR )
      !   call Show ( Delta_J_J, '>>> Delta_J_J', CONSOLE % ERROR )
      !   call Show ( Delta_H_H, '>>> Delta_H_H', CONSOLE % ERROR )
      ! end if

      !-- FIXME: Do solve above

      J ( iV )  =  max ( E ( iV ), SqrtTiny )
      
      H_1 ( iV )  =  M_UU_11 ( iV )  *  S_1 ( iV )
      H_2 ( iV )  =  M_UU_22 ( iV )  *  S_2 ( iV )
      H_3 ( iV )  =  M_UU_33 ( iV )  *  S_3 ( iV )          

      N ( iV )  =  max ( G ( iV ), SqrtTiny )
      
      !-- Moment factors ( Minerbo SF )

      H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                   +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                   +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

      FF ( iV )  =  H  /  J ( iV )

      SF ( iV )  =  1.0_KDR / 3.0_KDR &
                    +  2.0_KDR / 3.0_KDR &
                       *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                            * ( 3.0_KDR  -  FF ( iV )  &
                                +  3.0_KDR  *  FF ( iV ) ** 2 ) )

    else

      J   ( iV )  =  SqrtTiny
      H_1 ( iV )  =  0.0_KDR
      H_2 ( iV )  =  0.0_KDR
      H_3 ( iV )  =  0.0_KDR
      N   ( iV )  =  SqrtTiny
      E   ( iV )  =  SqrtTiny
      S_1 ( iV )  =  0.0_KDR
      S_2 ( iV )  =  0.0_KDR
      S_3 ( iV )  =  0.0_KDR
      G   ( iV )  =  SqrtTiny
      FF  ( iV )  =  0.0_KDR
      SF  ( iV )  =  1.0_KDR / 3.0_KDR

      return

    end if

    SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                     / ( 1.0_KDR / 3.0_KDR )

    H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                 +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                 +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
    
    if ( H  >  J ( iV ) ) then

      H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
      H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
      H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )

      H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                   +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                   +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

      !-- Moment factors ( Minerbo SF )

      FF ( iV )  =  H  /  J ( iV )

      SF ( iV )  =  1.0_KDR / 3.0_KDR &
                    +  2.0_KDR / 3.0_KDR &
                       *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                            * ( 3.0_KDR  -  FF ( iV )  &
                                +  3.0_KDR  *  FF ( iV ) ** 2 ) )

      SF_RD ( iV )  =  abs ( SF ( iV )  -  1.0_KDR / 3.0_KDR )  &
                       / ( 1.0_KDR / 3.0_KDR )

      !-- FIXME: Add velocity dependence

      E ( iV )  =  J ( iV )

      S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
      S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
      S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

      G ( iV )  =  N ( iV )

    end if

  end procedure Compute_J_H_N_G_S_Kernel


  module procedure Compute_SP_A_Kernel

    !-- Compute_SpectralParameters_All_xKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
!      SixPi_2, EightPi_2, &
!      OnePlusEpsilon, &
!      Factor_ND, Factor_ED_1, Factor_ED_2, &
      Factor_ND, &
!      LHS, &
!      Eta_ND, Eta_ED, Eta, &
      Eta_ND, EtaMax, &
      F_2, F_3!, &
 !     fdeta, fdeta2, &
 !     fdtheta, fdtheta2, &
 !     fdetadtheta
    logical ( KDL ) :: &
      UseDevice, &
      Bracket, &
      Converge
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J )

    ! OnePlusEpsilon  =  1.0_KDR  +  10.0_KDR * epsilon ( 0.0_KDR )

    !   SixPi_2  =  6.0_KDR  *  Pi_2
    ! EightPi_2  =  8.0_KDR  *  Pi_2

    Factor_ND    =  Pi ** ( - 2.0_KDR / 3.0_KDR )  /  3.0_KDR
 !   Factor_ED_1  =  EightPi_2 ** ( 1.0_KDR / 4.0_KDR )  &
 !                   /  SixPi_2 ** ( 1.0_KDR / 3.0_KDR )
 !   Factor_ED_2  =  6.0_KDR  /  Pi ** 2

    EtaMax  =  25.

    if ( UseDevice ) then
  !     !$OMP OMP_TARGET_DIRECTIVE parallel do &
  !     !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
  !     !$OMP shared ( a )
  !     do iV = 1, nV
  !       T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
  !     end do
  !     !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
!      !$OMP shared ( OnePlusEpsilon, Factor_ND, Factor_ED_1, Factor_ED_2 ) &
!      !$OMP private ( LHS, Eta_ND, Eta_ED, Eta, F_2, F_3 ) &
!      !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta ) &
!      !$OMP private ( Success )
      !$OMP shared ( Factor_ND, EtaMax ) &
      !$OMP private ( Eta_ND, F_2, F_3, Bracket, Converge )
      do iV = 1, nV

        if ( J ( iV )  <=  0.0_KDR  .or.  N ( iV )  <=  0.0_KDR ) &
          cycle

        ! LHS  =  J ( iV ) ** ( 1.0_KDR / 4.0_KDR )  &
        !         *  N ( iV ) ** ( - 1.0_KDR / 3.0_KDR )
        ! LHS  =  max ( LHS, OnePlusEpsilon / Factor_ED_1 )

        ! Eta_ND  =  - 3.0_KDR  * log ( Factor_ND  *  LHS ** 4 )
        ! Eta_ED  =  ( Factor_ED_2 * ( Factor_ED_1 * LHS  -  1.0_KDR ) ) &
        !            ** ( - 0.5_KDR )

        ! if ( Eta_ND < -10.0_KDR ) then
        !   Eta_R ( iV )  =  Eta_ND
        ! else if ( Eta_ED > 50.0_KDR ) then
        !   Eta_R ( iV )  =  Eta_ED
        ! else
        !   Eta  =  max ( Eta_ND, Eta_R ( iV ) )
        !   call SolveSecant ( LHS, 0.99 * Eta, Eta, Success, Eta )
        !   if ( Success ) then
        !     Eta_R ( iV )  =  Eta
        !   else
        !     !-- crude last resort
        !     if ( Eta_ND < 1.0_KDR ) then
        !       Eta_R ( iV )  =  Eta_ND
        !     else
        !       Eta_R ( iV )  =  Eta_ED
        !     end if
        !   end if
        ! end if

        Eta_ND  =  - 3.0_KDR  &
                     * log ( Factor_ND  *  J ( iV ) &
                             *  N ( iV ) ** ( - 4.0_KDR / 3.0_KDR ) )

        if ( Eta_ND  <=  0.0_KDR ) then
          Eta_R ( iV )  =  Eta_ND
        else
          Eta_R ( iV )  =  max ( Eta_R ( iV ), Eta_ND )
          call SolveEtaBisection &
                 ( Eta_R ( iV ), J ( iV ), N ( iV ), EtaMax, iV, &
                   Bracket, Converge )
          Eta_R ( iV )  =  min ( Eta_R ( iV ), EtaMax )
          if ( .not. Bracket .or. .not. Converge ) then
    !        call Show ( Eta_0, 'Eta_0', CONSOLE % ERROR )
    !        call Show ( Eta_ND, 'Eta_ND', CONSOLE % ERROR )
    !        Eta_R ( iV )  =  Eta_0
    !        Eta_R ( iV )  =  0.0_KDR
            Eta_R ( iV )  =  EtaMax
    !        call PROGRAM_HEADER % Abort ( )
          end if
        end if

        ! call DFERMI ( 2.0_KDR, Eta_R ( iV ), 0.0_KDR, F_2, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_R ( iV ), 0.0_KDR, F_3, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        F_2  =  Fermi_2 ( Eta_R ( iV ) )
        F_3  =  Fermi_3 ( Eta_R ( iV ) )

        T_R  ( iV )  =  J ( iV )  /  N ( iV )  *  F_2 / F_3

        E_Ave ( iV )  &
          =  J ( iV )  /  N ( iV )
        F_Ave ( iV )  &
          =  1.0_KDR &
             /  ( exp ( E_Ave ( iV ) / T_R ( iV )  -  Eta_R ( iV ) )  &
                  +  1.0_KDR )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_SP_A_Kernel


  module procedure Compute_SP_S_Kernel

    !-- Compute_SpectralParameters_Singlex_Kernel

    real ( KDR ) :: &
!      SixPi_2, EightPi_2, &
!      OnePlusEpsilon, &
!      Factor_ND, Factor_ED_1, Factor_ED_2, &
      Factor_ND, &
!      LHS, &
!      Eta_ND, Eta_ED, Eta, &
      Eta_ND, EtaMax, &
      F_2, F_3!, &
 !     fdeta, fdeta2, &
 !     fdtheta, fdtheta2, &
 !     fdetadtheta
    logical ( KDL ) :: &
      Bracket, &
      Converge
          
    ! OnePlusEpsilon  =  1.0_KDR  +  10.0_KDR * epsilon ( 0.0_KDR )

    !   SixPi_2  =  6.0_KDR  *  Pi_2
    ! EightPi_2  =  8.0_KDR  *  Pi_2

    Factor_ND    =  Pi ** ( - 2.0_KDR / 3.0_KDR )  /  3.0_KDR
 !   Factor_ED_1  =  EightPi_2 ** ( 1.0_KDR / 4.0_KDR )  &
 !                   /  SixPi_2 ** ( 1.0_KDR / 3.0_KDR )
 !   Factor_ED_2  =  6.0_KDR  /  Pi ** 2

    EtaMax  =  25.

    if ( J ( iV )  <=  0.0_KDR  .or.  N ( iV )  <=  0.0_KDR ) &
      return

    ! LHS  =  J ( iV ) ** ( 1.0_KDR / 4.0_KDR )  &
    !         *  N ( iV ) ** ( - 1.0_KDR / 3.0_KDR )
    ! LHS  =  max ( LHS, OnePlusEpsilon / Factor_ED_1 )

    ! Eta_ND  =  - 3.0_KDR  * log ( Factor_ND  *  LHS ** 4 )
    ! Eta_ED  =  ( Factor_ED_2 * ( Factor_ED_1 * LHS  -  1.0_KDR ) ) &
    !            ** ( - 0.5_KDR )

    ! if ( Eta_ND < -10.0_KDR ) then
    !   Eta_R ( iV )  =  Eta_ND
    ! else if ( Eta_ED > 50.0_KDR ) then
    !   Eta_R ( iV )  =  Eta_ED
    ! else
    !   Eta  =  max ( Eta_ND, Eta_R ( iV ) )
    !   call SolveSecant ( LHS, 0.99 * Eta, Eta, Success, Eta )
    !   if ( Success ) then
    !     Eta_R ( iV )  =  Eta
    !   else
    !     !-- crude last resort
    !     if ( Eta_ND < 1.0_KDR ) then
    !       Eta_R ( iV )  =  Eta_ND
    !     else
    !       Eta_R ( iV )  =  Eta_ED
    !     end if
    !   end if
    ! end if

    Eta_ND  =  - 3.0_KDR  &
                 * log ( Factor_ND  *  J ( iV ) &
                         *  N ( iV ) ** ( - 4.0_KDR / 3.0_KDR ) )

    if ( Eta_ND  <=  0.0_KDR ) then
      Eta_R ( iV )  =  Eta_ND
    else
      Eta_R ( iV )  =  max ( Eta_R ( iV ), Eta_ND )
      call SolveEtaBisection &
             ( Eta_R ( iV ), J ( iV ), N ( iV ), EtaMax, iV, Bracket, Converge )
      Eta_R ( iV )  =  min ( Eta_R ( iV ), EtaMax )
      if ( .not. Bracket .or. .not. Converge ) then
!        call Show ( Eta_0, 'Eta_0', CONSOLE % ERROR )
!        call Show ( Eta_ND, 'Eta_ND', CONSOLE % ERROR )
!        Eta_R ( iV )  =  Eta_0
!        Eta_R ( iV )  =  0.0_KDR
        Eta_R ( iV )  =  EtaMax
!        call PROGRAM_HEADER % Abort ( )
      end if
    end if

    ! call DFERMI ( 2.0_KDR, Eta_R ( iV ), 0.0_KDR, F_2, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 3.0_KDR, Eta_R ( iV ), 0.0_KDR, F_3, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    F_2  =  Fermi_2 ( Eta_R ( iV ) )
    F_3  =  Fermi_3 ( Eta_R ( iV ) )

    T_R  ( iV )  =  J ( iV )  /  N ( iV )  *  F_2 / F_3

    E_Ave ( iV )  &
      =  J ( iV )  /  N ( iV )
    F_Ave ( iV )  &
      =  1.0_KDR &
         /  ( exp ( E_Ave ( iV ) / T_R ( iV )  -  Eta_R ( iV ) )  &
              +  1.0_KDR )

  end procedure Compute_SP_S_Kernel


  module procedure Compute_Eq_E_A_Kernel

    !-- Compute_Equilibrium_E_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      TwoPi, FourPi, & 
      Factor_J_N, &
      Eta_Eq, &
      F_2_Eq, F_3_Eq!, &
!      fdeta, fdeta2, &
!      fdtheta, fdtheta2, &
!      fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

     TwoPi  =  2.0_KDR  *  Pi
    FourPi  =  4.0_KDR  *  Pi

    Factor_J_N   =  FourPi  /  TwoPi ** 3

    if ( UseDevice ) then
    !   !$OMP OMP_TARGET_DIRECTIVE parallel do &
    !   !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
    !   !$OMP shared ( a )
    !   do iV = 1, nV
    !     J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
    !   end do
    !   !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny, Factor_J_N ) &
      !$OMP private ( Eta_Eq, F_2_Eq, F_3_Eq ) !&
 !     !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV )  <=  0.0_KDR ) &
          cycle

        Eta_Eq  =  Sign  *  ( Mu_E ( iV )  -  Mu_NP ( iV ) )  /  T ( iV )
        
        ! call DFERMI ( 2.0_KDR, Eta_Eq, 0.0_KDR, F_2_Eq, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_Eq, 0.0_KDR, F_3_Eq, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        F_2_Eq  =  Fermi_2 ( Eta_Eq )
        F_3_Eq  =  Fermi_3 ( Eta_Eq )

        N_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 3  *  F_2_Eq
        J_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 4  *  F_3_Eq

        N_RD  ( iV )  =  abs ( N ( iV )  -  N_Eq ( iV ) )  &
                         /  max ( SqrtTiny, N_Eq ( iV ) )
        J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                         /  max ( SqrtTiny, J_Eq ( iV ) )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_Eq_E_A_Kernel


  module procedure Compute_Eq_E_S_Kernel

    !-- Compute_Equilibrium_E_Single_Kernel

    real ( KDR ) :: &
      SqrtTiny, &
      TwoPi, FourPi, & 
      Factor_J_N, &
      Eta_Eq, &
      F_2_Eq, F_3_Eq!, &
!      fdeta, fdeta2, &
!      fdtheta, fdtheta2, &
!      fdetadtheta

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

     TwoPi  =  2.0_KDR  *  Pi
    FourPi  =  4.0_KDR  *  Pi

    Factor_J_N   =  FourPi  /  TwoPi ** 3

    if ( T ( iV )  <=  0.0_KDR ) &
      return

    Eta_Eq  =  Sign  *  ( Mu_E ( iV )  -  Mu_NP ( iV ) )  /  T ( iV )

! Eta_Eq  =  min ( Eta_Eq,   50.0_KDR )
! Eta_Eq  =  max ( Eta_Eq, - 50.0_KDR )
! call Show ( Mu_E ( iV ), '>>> Mu_E' )
! call Show ( Mu_NP ( iV ), '>>> Mu_NP' )
! call Show ( T ( iV ), '>>> T' )
! call Show ( Eta_Eq, '>>> Eta_Eq' )
    
    ! call DFERMI ( 2.0_KDR, Eta_Eq, 0.0_KDR, F_2_Eq, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 3.0_KDR, Eta_Eq, 0.0_KDR, F_3_Eq, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    F_2_Eq  =  Fermi_2 ( Eta_Eq )
    F_3_Eq  =  Fermi_3 ( Eta_Eq )

!call Show ( F_2_Eq, '>>> F_2_Eq' )
!call Show ( F_3_Eq, '>>> F_3_Eq' )

    N_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 3  *  F_2_Eq
    J_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 4  *  F_3_Eq

!call Show ( N_Eq ( iV ), '>>> N_Eq' )
!call Show ( J_Eq ( iV ), '>>> J_Eq' )

    N_RD  ( iV )  =  abs ( N ( iV )  -  N_Eq ( iV ) )  &
                     /  max ( SqrtTiny, N_Eq ( iV ) )
    J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                     /  max ( SqrtTiny, J_Eq ( iV ) )

  end procedure Compute_Eq_E_S_Kernel


  module procedure Compute_Eq_HL_A_Kernel

    !-- Compute_Equilibrium_HL_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      TwoPi, FourPi, & 
      Factor_J_N, &
      F_2_Eq, F_3_Eq!, &
!      fdeta, fdeta2, &
!      fdtheta, fdtheta2, &
!      fdetadtheta
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

     TwoPi  =  2.0_KDR  *  Pi
    FourPi  =  4.0_KDR  *  Pi

    Factor_J_N   =  nSpecies  *  FourPi  /  TwoPi ** 3

    if ( UseDevice ) then
    !   !$OMP OMP_TARGET_DIRECTIVE parallel do &
    !   !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
    !   !$OMP shared ( a )
    !   do iV = 1, nV
    !     J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
    !   end do
    !   !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny, Factor_J_N ) &
      !$OMP private ( F_2_Eq, F_3_Eq ) !&
 !     !$OMP private ( fdeta, fdeta2, fdtheta, fdtheta2, fdetadtheta )
      do iV = 1, nV

        if ( T ( iV )  <=  0.0_KDR ) &
          cycle

        ! call DFERMI ( 2.0_KDR, Eta_Eq, 0.0_KDR, F_2_Eq, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        ! call DFERMI ( 3.0_KDR, Eta_Eq, 0.0_KDR, F_3_Eq, &
        !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
        F_2_Eq  =  Fermi_2 ( 0.0_KDR )
        F_3_Eq  =  Fermi_3 ( 0.0_KDR )

        N_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 3  *  F_2_Eq
        J_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 4  *  F_3_Eq

        N_RD  ( iV )  =  abs ( N ( iV )  -  N_Eq ( iV ) )  &
                         /  max ( SqrtTiny, N_Eq ( iV ) )
        J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                         /  max ( SqrtTiny, J_Eq ( iV ) )

      end do
      !$OMP end parallel do
    end if

  end procedure Compute_Eq_HL_A_Kernel


  module procedure Compute_Eq_HL_S_Kernel

    !-- Compute_Equilibrium_HL_Single_Kernel

    real ( KDR ) :: &
      SqrtTiny, &
      TwoPi, FourPi, & 
      Factor_J_N, &
      F_2_Eq, F_3_Eq!, &
!      fdeta, fdeta2, &
!      fdtheta, fdtheta2, &
!      fdetadtheta

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

     TwoPi  =  2.0_KDR  *  Pi
    FourPi  =  4.0_KDR  *  Pi

    Factor_J_N   =  nSpecies  *  FourPi  /  TwoPi ** 3

    if ( T ( iV )  <=  0.0_KDR ) &
      return

! Eta_Eq  =  min ( Eta_Eq,   50.0_KDR )
! Eta_Eq  =  max ( Eta_Eq, - 50.0_KDR )
! call Show ( Mu_E ( iV ), '>>> Mu_E' )
! call Show ( Mu_NP ( iV ), '>>> Mu_NP' )
! call Show ( T ( iV ), '>>> T' )
! call Show ( Eta_Eq, '>>> Eta_Eq' )
    
    ! call DFERMI ( 2.0_KDR, Eta_Eq, 0.0_KDR, F_2_Eq, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    ! call DFERMI ( 3.0_KDR, Eta_Eq, 0.0_KDR, F_3_Eq, &
    !               fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
    F_2_Eq  =  Fermi_2 ( 0.0_KDR )
    F_3_Eq  =  Fermi_3 ( 0.0_KDR )

!call Show ( F_2_Eq, '>>> F_2_Eq' )
!call Show ( F_3_Eq, '>>> F_3_Eq' )

    N_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 3  *  F_2_Eq
    J_Eq ( iV )  =  Factor_J_N  *  T ( iV ) ** 4  *  F_3_Eq

!call Show ( N_Eq ( iV ), '>>> N_Eq' )
!call Show ( J_Eq ( iV ), '>>> J_Eq' )

    N_RD  ( iV )  =  abs ( N ( iV )  -  N_Eq ( iV ) )  &
                     /  max ( SqrtTiny, N_Eq ( iV ) )
    J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                     /  max ( SqrtTiny, J_Eq ( iV ) )

  end procedure Compute_Eq_HL_S_Kernel


!   subroutine SolveSecant ( LHS, Guess_1, Guess_2, Success, Root )
  
!     real ( KDR ), intent ( in ) :: &
!       LHS, &
!       Guess_1, &
!       Guess_2
!     logical ( KDL ), intent ( out ) :: &
!       Success
!     real ( KDR ), intent ( out ) :: &
!       Root
      
!     integer ( KDI ) :: &
!       iIteration, &
!       MaxIterations
!     real ( KDR ) :: &
!       X, &
!       X_0, X_1, &
!       Y_0, Y_1, &
!       RequestedAccuracy, &
!       Accuracy, &
!       AbsolutePrecision, &
!       RelativePrecision
      
!     MaxIterations      =  20
!     RequestedAccuracy  =  1.0e-10_KDR !epsilon ( 1.0_KDR ) * 10.0_KDR 

!     Success  =  .false.
      
!     X_0  =  Guess_1
!     X_1  =  Guess_2
    
! !    call EvaluateZero ( LHS, X_0, Y_0 )
! !    call EvaluateZero ( LHS, X_1, Y_1 )
    
!     do iIteration  =  1, MaxIterations
      
!       Root  =  X_1
      
!       Accuracy           =  abs ( Y_1 )
!       AbsolutePrecision  =  abs ( X_1 - X_0 )
!       RelativePrecision  =  abs ( ( X_1 - X_0 ) )  &
!                             /  max ( abs ( X_1 ), tiny ( 0.0_KDR ) )

!       if (      Accuracy           <=  RequestedAccuracy   &
!            .or. AbsolutePrecision  <=  RequestedAccuracy   &
!            .or. RelativePrecision  <=  RequestedAccuracy ) &
!       then
!         Success = .true. 
!         exit
!       end if
      
!       if ( Y_1 == Y_0 ) &
!         exit

!       X = X_1 - Y_1 * ( X_1 - X_0 ) / ( Y_1 - Y_0 )
      
!       X_0 = X_1
!       Y_0 = Y_1
      
!       X_1 = X
! !      call EvaluateZero ( LHS, X_1, Y_1 )

!     end do
    
!   end subroutine SolveSecant
  

  subroutine SolveEtaBisection ( Eta, J, N, EtaMax, iV, Bracket, Converge )

    real ( KDR ), intent ( inout ) :: &
      Eta
    real ( KDR ), intent ( in ) :: &
      J, N
    real ( KDR ), intent ( in ) :: &
      EtaMax
    integer ( KDI ), intent ( in ) :: &
      iV
    logical ( KDL ), intent ( out ) :: &
      Bracket, &
      Converge

    integer ( KDI ) :: &
      iIteration, &
      MaxIterations
    real ( KDR ) :: &
      dX, &
      X_0, X_1, X_M, &
      Y_0, Y_1, Y_M, &
      Factor, &
      Tolerance, &
      AbsolutePrecision, &
      RelativePrecision
      
    MaxIterations  =  50
    Tolerance      =  1.0e-9_KDR !epsilon ( 1.0_KDR ) * 10.0_KDR 
    Factor         =  1.6_KDR

    X_0  =  0.0
    X_1  =  2.0 * Eta
    
    Y_0  =  ZeroEta ( J, N, X_0 )
    Y_1  =  ZeroEta ( J, N, X_1 )

    !-- Bracket root

    Bracket  =  .false.

! if ( iV == 3 ) then
!   call Show ( '>>> Bracketing Eta' ) 
!   call Show ( [ X_0, Y_0 ], '>>> [ X_0, Y_0 ]' )
!   call Show ( [ X_1, Y_1 ], '>>> [ X_1, Y_1 ]' )
! end if
      
    do iIteration  =  1, MaxIterations
      !-- Only move outer bound
      X_1  =  Factor * X_1
      Y_1  =  ZeroEta ( J, N, X_1 )
! if ( iV == 3 ) then
!   call Show ( iIteration, '>>> iIteration' )
!   call Show ( [ X_1, Y_1 ], '>>> [ X_1, Y_1 ]' )
! end if      
      if ( Y_0 * Y_1  <  0.0_KDR ) then
        Bracket  =  .true.
        exit
      else if ( X_1  >  EtaMax ) then
        exit
      end if
    end do

    if ( .not. Bracket ) then
      ! call Show ( 'Failed to bracket Eta', CONSOLE % ERROR )
      ! call Show ( iV, 'iV', CONSOLE % ERROR )
      ! call Show ( J, 'J', CONSOLE % ERROR )
      ! call Show ( N, 'N', CONSOLE % ERROR )
      ! call Show ( Eta, 'Eta', CONSOLE % ERROR )
!      call PROGRAM_HEADER % Abort ( )
      return
    end if

    !-- Find root

    Converge  =  .false.

    !-- Orient the search
    if ( Y_0  <  0.0_KDR ) then
      Eta  =  X_0
       dX  =  X_1 - X_0
    else !-- Y_1  <  0.0
      Eta  =  X_1
       dX  =  X_0 - X_1
    end if

    do iIteration  =  1, MaxIterations

      AbsolutePrecision  =  abs ( dX )
      RelativePrecision  =  abs ( dX )  &
                            /  max ( abs ( Eta ), abs ( Eta + dX ) )

!      if (      AbsolutePrecision  <=  Tolerance   &
!           .or. RelativePrecision  <=  Tolerance ) &
!      then
      if ( RelativePrecision  <=  Tolerance ) then
! if ( iV == 3 ) then
!   call Show ( '>>> Converged' )
!   call Show ( iV, '>>> iV' )
!   call Show ( iIteration, '>>> iIteration' )
!   call Show ( Eta, '>>> Eta' )
! end if      
        Converge = .true. 
        exit
      end if

       dX  =  0.5_KDR * dX
      X_M  =  Eta + dX
      Y_M  =  ZeroEta ( J, N, X_M )
      if ( Y_M  <  0.0_KDR ) &
        Eta  =  X_M
      
    end do
    
    if ( .not. Converge ) then
      ! call Show ( 'Failed to converge Eta', CONSOLE % ERROR )
      ! call Show ( iV, 'iV', CONSOLE % ERROR )
      ! call Show ( J, 'J', CONSOLE % ERROR )
      ! call Show ( N, 'N', CONSOLE % ERROR )
      ! call Show ( Eta, 'Eta', CONSOLE % ERROR )
      ! call Show ( AbsolutePrecision, 'AbsolutePrecision', CONSOLE % ERROR )
      ! call Show ( RelativePrecision, 'RelativePrecision', CONSOLE % ERROR )
!      call PROGRAM_HEADER % Abort ( )
      return
    end if

  end subroutine SolveEtaBisection


  ! subroutine EvaluateZero ( LHS, EtaIn, Result )

  !   real ( KDR ), intent ( in ) :: &
  !     LHS, &
  !     EtaIn
  !   real ( KDR ), intent ( out ) :: &
  !     Result

  !   real ( KDR ) :: &
  !     Factor, &
  !     Eta, &
  !     F_2, F_3, &
  !     fdeta, fdeta2, &
  !     fdtheta, fdtheta2, &
  !     fdetadtheta

  !   Factor  =  ( 2  *  Pi ** 2 ) ** ( 1. / 12. )

  !   Eta  =  min ( EtaIn, log ( huge ( 1.0_KDR ) ) - 10.0_KDR )

  !   call DFERMI ( 2.0_KDR, Eta, 0.0_KDR, F_2, &
  !                 fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
  !   call DFERMI ( 3.0_KDR, Eta, 0.0_KDR, F_3, &
  !                 fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

  !   F_2  =  max ( F_2, tiny ( 0.0_KDR ) )

  !   Result  =  Factor  *  F_3 ** ( 1.0_KDR / 4.0_KDR )  &
  !                 *  F_2 ** ( - 1.0_KDR / 3.0_KDR ) &
  !              -  LHS

  ! end subroutine EvaluateZero


  function ZeroEta ( J, N, Eta ) result ( ZE )

    real ( KDR ), intent ( in ) :: &
      J, N, Eta
    real ( KDR ) :: &
      ZE

    real ( KDR ) :: &
      FourThirds, &
      Factor, &
      F_2, F_3

    FourThirds  =  4.0_KDR / 3.0_KDR
    Factor      =  ( 2 * Pi_2 ) ** ( 1.0_KDR / 3.0_KDR )

    F_2  =  Fermi_2 ( Eta )
    F_3  =  Fermi_3 ( Eta )

    ZE  =  J  *  F_2 ** FourThirds   -   Factor  *  N ** FourThirds  *  F_3
   
  end function ZeroEta


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

    if ( Eta  >  0.0_KDR ) then
      F_3  =  Eta**4 / 4.  +  Pi_2 * Eta**2 / 2.  +  12.  -  6. * exp ( -Eta )
    else
      F_3  =  6. * exp ( Eta )
    end if

  end function Fermi_3


end submodule NeutrinoMoments_G__Kernel
