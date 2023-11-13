#include "Preprocessor"

submodule ( DivergencePart_NM_G__Form ) DivergencePart_NM_G__Kernel

  use Basics

  implicit none

contains 


  module procedure Compute_FS_G_Kernel

    !-- Compute_FluxSet_Galileo_Kernel

    integer :: &
      Delta_1, Delta_2, Delta_3
    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H_Sq, &
      K_U_Dim_D_1, K_U_Dim_D_2, K_U_Dim_D_3, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( F_E )
    
    Delta_1  =  int (    real ( 1 + iDim  -  abs ( 1 - iDim ) )  &
                      /  real ( 1 + iDim  +  abs ( 1 - iDim ) ) )

    Delta_2  =  int (    real ( 2 + iDim  -  abs ( 2 - iDim ) )  &
                      /  real ( 2 + iDim  +  abs ( 2 - iDim ) ) )

    Delta_3  =  int (    real ( 3 + iDim  -  abs ( 3 - iDim ) )  &
                      /  real ( 3 + iDim  +  abs ( 3 - iDim ) ) )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3, SqrtTiny ) &
      !$OMP private ( H_Sq, K_U_Dim_D_1, K_U_Dim_D_2, K_U_Dim_D_3 )
      do iV = 1, nV

        !-- Energy

        F_E ( iV )  =  M_UU_Dim ( iV )  *  S_Dim ( iV ) 

        !-- Momentum

        H_Sq  =     M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                 +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                 +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2
        H_Sq  =  max ( H_Sq, SqrtTiny )

        K_U_Dim_D_1  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_1  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_11 ( iV ) * H_1 ( iV )  /  H_Sq )  &
             *  J ( iV )

        K_U_Dim_D_2  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_2  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_22 ( iV ) * H_2 ( iV )  /  H_Sq )  &
             *  J ( iV )

        K_U_Dim_D_3  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_3  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_33 ( iV ) * H_3 ( iV )  /  H_Sq )  &
             *  J ( iV )

        !-- FIXME: Add velocity dependence

        F_S_1 ( iV )  =  K_U_Dim_D_1
        F_S_2 ( iV )  =  K_U_Dim_D_2
        F_S_3 ( iV )  =  K_U_Dim_D_3

        !-- Number
        !-- FIXME: Add velocity dependence

        F_D ( iV )  =  ( H_Dim ( iV )  /  max ( J ( iV ), SqrtTiny ) )  &
                       *  N ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3, SqrtTiny ) &
      !$OMP private ( H_Sq, K_U_Dim_D_1, K_U_Dim_D_2, K_U_Dim_D_3 )
      do iV = 1, nV

        !-- Energy

        F_E ( iV )  =  M_UU_Dim ( iV )  *  S_Dim ( iV ) 

        !-- Momentum

        H_Sq  =     M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                 +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                 +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2
        H_Sq  =  max ( H_Sq, SqrtTiny )

        K_U_Dim_D_1  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_1  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_11 ( iV ) * H_1 ( iV )  /  H_Sq )  &
             *  J ( iV )

        K_U_Dim_D_2  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_2  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_22 ( iV ) * H_2 ( iV )  /  H_Sq )  &
             *  J ( iV )

        K_U_Dim_D_3  &
          =  0.5_KDR * ( ( 1.0_KDR  -  SF ( iV ) )  *  Delta_3  &
                         +  ( 3.0_KDR * SF ( iV )  -  1.0_KDR )  &
                            *  H_Dim ( iV )  &
                            *  M_DD_33 ( iV ) * H_3 ( iV )  /  H_Sq )  &
             *  J ( iV )

        !-- FIXME: Add velocity dependence

        F_S_1 ( iV )  =  K_U_Dim_D_1
        F_S_2 ( iV )  =  K_U_Dim_D_2
        F_S_3 ( iV )  =  K_U_Dim_D_3

        !-- Number
        !-- FIXME: Add velocity dependence

        F_D ( iV )  =  ( H_Dim ( iV )  /  max ( J ( iV ), SqrtTiny ) )  &
                       *  N ( iV )

      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_FS_G_Kernel


end submodule DivergencePart_NM_G__Kernel
