#include "Preprocessor"

submodule ( RadiationMoments_BM__Form ) RadiationMoments_BM__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_E_S_G_Kernel

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

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_E_S_G_Kernel


  module procedure Compute_J_H_G_A_Kernel

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

        if ( E ( iV )  >  0.0_KDR ) then

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
          E   ( iV )  =  SqrtTiny
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
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

        end if

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( E ( iV )  >  0.0_KDR ) then

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
          E   ( iV )  =  SqrtTiny
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
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

        end if

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_J_H_G_A_Kernel


  module procedure Compute_J_H_G_S_Kernel

    !-- Compute_ComovingEnergy_Momentum_Galileo_Single_Kernel

    real ( KDR ) :: &
      H, &
      SqrtTiny
      
#ifdef ENABLE_OMP_OFFLOAD  
    !$OMP declare target
#endif

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( E ( iV )  >  0.0_KDR ) then

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
      E   ( iV )  =  SqrtTiny
      S_1 ( iV )  =  0.0_KDR
      S_2 ( iV )  =  0.0_KDR
      S_3 ( iV )  =  0.0_KDR
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

    end if

  end procedure Compute_J_H_G_S_Kernel


  module procedure Compute_ES_G_Kernel

    !-- Compute_EigenspeedSet_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H, &
      H_Dim_H, &
      c_D, &
      SF, SFP, &
      D, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( EF_P )

    SqrtTiny   =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( H, H_Dim_H, c_D, SF, SFP, D ) &
      !$OMP shared ( SqrtTiny )
      do iV = 1, nV
        ! !-- Speed of light streaming

        ! EF_P ( iV )  =    sqrt ( M_UU_Dim ( iV ) )  *  c
        ! EF_M ( iV )  =  - sqrt ( M_UU_Dim ( iV ) )  *  c

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
        H  =  max ( H, SqrtTiny )

        H_Dim_H  =  H_Dim ( iV )  /  H
        
        ! !-- A convex combination of streaming and diffusion speeds
        ! !        is my own hack, not actual eigenvalues.

        ! c_D  =  sqrt ( M_UU_Dim ( iV ) )  *  c  /  sqrt ( 3.0_KDR )

        ! EF_P ( iV )  =  FF ( iV )**2  *  H_Dim_H  *  c  &
        !                 +  ( 1.0_KDR  -  FF ( iV )**2 )  *  c_D
        ! EF_M ( iV )  =  FF ( iV )**2  *  H_Dim_H  *  c  &
        !                 -  ( 1.0_KDR  -  FF ( iV )**2 )  *  c_D

        !-- 1D eigenvalues from Endeve et al. 2017

        SF  =  1.0_KDR / 3.0_KDR  &
               +  2.0_KDR / 15.0_KDR  &
                  *  FF ( iV ) ** 2 &
                     *  ( 3.0_KDR  -  FF ( iV )  +  3.0_KDR  *  FF ( iV ) ** 2 )

        SFP  =  2.0_KDR / 5.0_KDR  *  FF ( iV )  &
                *  ( 2.0_KDR  -  FF ( iV )  +  4.0_KDR  *  FF ( iV ) ** 2 )

        D  =  max ( ( SFP  -  2.0_KDR  *  FF ( iV ) ) ** 2  &
                    +  4.0_KDR * ( SF  -  FF ( iV ) ** 2 ), &
                    tiny ( 0.0_KDR ) )

        EF_P ( iV )  =  ( H_Dim_H * SFP  &
                          +  sqrt ( M_UU_Dim ( iV )  * D ) )  /  2.0_KDR 
        EF_M ( iV )  =  ( H_Dim_H * SFP  &
                          -  sqrt ( M_UU_Dim ( iV )  * D ) )  /  2.0_KDR 

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( H, H_Dim_H, c_D, SF, SFP, D ) &
      !$OMP shared ( SqrtTiny )
      do iV = 1, nV

        ! !-- Speed of light streaming

        ! EF_P ( iV )  =    sqrt ( M_UU_Dim ( iV ) )  *  c
        ! EF_M ( iV )  =  - sqrt ( M_UU_Dim ( iV ) )  *  c

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )
        H  =  max ( H, SqrtTiny )

        H_Dim_H  =  H_Dim ( iV )  /  H
        
        ! !-- A convex combination of streaming and diffusion speeds
        ! !        is my own hack, not actual eigenvalues.

        ! c_D  =  sqrt ( M_UU_Dim ( iV ) )  *  c  /  sqrt ( 3.0_KDR )

        ! EF_P ( iV )  =  FF ( iV )**2  *  H_Dim_H  *  c  &
        !                 +  ( 1.0_KDR  -  FF ( iV )**2 )  *  c_D
        ! EF_M ( iV )  =  FF ( iV )**2  *  H_Dim_H  *  c  &
        !                 -  ( 1.0_KDR  -  FF ( iV )**2 )  *  c_D

        !-- 1D eigenvalues from Endeve et al. 2017

        SF  =  1.0_KDR / 3.0_KDR  &
               +  2.0_KDR / 15.0_KDR  &
                  *  FF ( iV ) ** 2 &
                     *  ( 3.0_KDR  -  FF ( iV )  +  3.0_KDR  *  FF ( iV ) ** 2 )

        SFP  =  2.0_KDR / 5.0_KDR  *  FF ( iV )  &
                *  ( 2.0_KDR  -  FF ( iV )  +  4.0_KDR  *  FF ( iV ) ** 2 )

        D  =  max ( ( SFP  -  2.0_KDR  *  FF ( iV ) ) ** 2  &
                    +  4.0_KDR * ( SF  -  FF ( iV ) ** 2 ), &
                    tiny ( 0.0_KDR ) )

        EF_P ( iV )  =  ( H_Dim_H * SFP  &
                          +  sqrt ( M_UU_Dim ( iV )  * D ) )  /  2.0_KDR 
        EF_M ( iV )  =  ( H_Dim_H * SFP  &
                          -  sqrt ( M_UU_Dim ( iV )  * D ) )  /  2.0_KDR 

      end do
      !$OMP end parallel do

    end if

  end procedure Compute_ES_G_Kernel


end submodule RadiationMoments_BM__Kernel
