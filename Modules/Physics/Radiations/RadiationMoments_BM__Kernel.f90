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
      H
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( E )

    if ( UseDevice ) then

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( J ( iV )  <  0.0_KDR ) &
          J ( iV )  =  0.0_KDR

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

        FF ( iV )  =  H  /  max ( J ( iV ), tiny ( 0.0_KDR ) )

        SF ( iV )  =  1.0_KDR / 3.0_KDR &
                      +  2.0_KDR / 3.0_KDR &
                         *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                              * ( 3.0_KDR  -  FF ( iV )  &
                                  +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        !-- FIXME: Add velocity dependence

        E ( iV )  =  J ( iV )

        S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
        S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
        S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_E_S_G_Kernel


  module procedure Compute_J_H_G_Kernel

    !-- Compute_ComovingEnergy_Momentum_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( J )

    if ( UseDevice ) then

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
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

          J ( iV )  =  E ( iV )
          
          H_1 ( iV )  =  M_UU_11 ( iV )  *  S_1 ( iV )
          H_2 ( iV )  =  M_UU_22 ( iV )  *  S_2 ( iV )
          H_3 ( iV )  =  M_UU_33 ( iV )  *  S_3 ( iV )          

          !-- Moment factors ( Minerbo SF )

          FF ( iV )  =  H  /  max ( J ( iV ), tiny ( 0.0_KDR ) )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

        else

          J   ( iV )  =  0.0_KDR
          H_1 ( iV )  =  0.0_KDR
          H_2 ( iV )  =  0.0_KDR
          H_3 ( iV )  =  0.0_KDR
          E   ( iV )  =  0.0_KDR
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          FF  ( iV )  =  0.0_KDR
          SF  ( iV )  =  0.0_KDR

          cycle 

        end if

        if ( J ( iV )  <  0.0_KDR ) then

          J   ( iV )  =  0.0_KDR
          H_1 ( iV )  =  0.0_KDR
          H_2 ( iV )  =  0.0_KDR
          H_3 ( iV )  =  0.0_KDR
          E   ( iV )  =  0.0_KDR
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          FF  ( iV )  =  0.0_KDR
          SF  ( iV )  =  0.0_KDR

          cycle

        end if

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

          FF ( iV )  =  H  /  max ( J ( iV ), tiny ( 0.0_KDR ) )

          SF ( iV )  =  1.0_KDR / 3.0_KDR &
                        +  2.0_KDR / 3.0_KDR &
                           *  ( FF ( iV ) ** 2  /  5.0_KDR  &
                                * ( 3.0_KDR  -  FF ( iV )  &
                                    +  3.0_KDR  *  FF ( iV ) ** 2 ) )

          !-- FIXME: Add velocity dependence

          E ( iV )  =  J ( iV )

          S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
          S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
          S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

        end if

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_J_H_G_Kernel


  module procedure Compute_ES_G_Kernel

    !-- Compute_EigenspeedSet_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( EF_P )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        EF_P ( iV )  =  + sqrt ( M_UU_Dim ( iV ) )  *  c 
        EF_M ( iV )  =  - sqrt ( M_UU_Dim ( iV ) )  *  c
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        EF_P ( iV )  =  + sqrt ( M_UU_Dim ( iV ) )  *  c 
        EF_M ( iV )  =  - sqrt ( M_UU_Dim ( iV ) )  *  c
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_ES_G_Kernel


end submodule RadiationMoments_BM__Kernel
