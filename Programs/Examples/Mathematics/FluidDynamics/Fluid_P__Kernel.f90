#include "Preprocessor"

submodule ( Fluid_P__Template ) Fluid_P__Kernel

  use Basics
  
  implicit none
  
contains


  module procedure ComputeConservedEnergyKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( G )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues 
        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeConservedEnergyKernel


  module procedure ComputeInternalEnergyKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( E )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        E ( iV )  =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                              +  S_2 ( iV ) * V_2 ( iV ) &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( E ( iV ) < 0.0_KDR ) then
          E ( iV ) = 0.0_KDR
          G ( iV ) = E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                              +  S_2 ( iV ) * V_2 ( iV ) &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
        end if
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        E ( iV )  =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                              +  S_2 ( iV ) * V_2 ( iV ) &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
      end do
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( E ( iV ) < 0.0_KDR ) then
          E ( iV ) = 0.0_KDR
          G ( iV ) = E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                              +  S_2 ( iV ) * V_2 ( iV ) &
                                              +  S_3 ( iV ) * V_3 ( iV ) )
        end if
      end do
      !$OMP end parallel do
    
    end if

  end procedure ComputeInternalEnergyKernel


  module procedure ComputeEigenspeedsFluidKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( FEP_1 )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
        else
          CS ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( P ( iV ) > 0.0_KDR ) then
          MN ( iV ) = sqrt ( (    S_1 ( iV ) * V_1 ( iV )  &
                               +  S_2 ( iV ) * V_2 ( iV )  &
                               +  S_3 ( iV ) * V_3 ( iV ) ) &
                             / ( Gamma ( iV ) * P ( iV ) ) )
        else
          MN ( iV ) = 0.0_KDR
        end if
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        FEP_1 ( iV )  =  V_1 ( iV )  +  CS ( iV ) 
        FEP_2 ( iV )  =  V_2 ( iV )  +  sqrt ( M_UU_22 ( iV ) ) * CS ( iV ) 
        FEP_3 ( iV )  =  V_3 ( iV )  +  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
        FEM_1 ( iV )  =  V_1 ( iV )  -  CS ( iV )
        FEM_2 ( iV )  =  V_2 ( iV )  -  sqrt ( M_UU_22 ( iV ) ) * CS ( iV )
        FEM_3 ( iV )  =  V_3 ( iV )  -  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else 
    
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then
          CS ( iV ) = sqrt ( Gamma ( iV ) * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
        else
          CS ( iV ) = 0.0_KDR
        end if 
      end do
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( P ( iV ) > 0.0_KDR ) then
          MN ( iV ) = sqrt ( (    S_1 ( iV ) * V_1 ( iV )  &
                               +  S_2 ( iV ) * V_2 ( iV )  &
                               +  S_3 ( iV ) * V_3 ( iV ) ) &
                             / ( Gamma ( iV ) * P ( iV ) ) )
        else
          MN ( iV ) = 0.0_KDR
        end if 
      end do
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        FEP_1 ( iV )  =  V_1 ( iV )  +  CS ( iV ) 
        FEP_2 ( iV )  =  V_2 ( iV )  +  sqrt ( M_UU_22 ( iV ) ) * CS ( iV ) 
        FEP_3 ( iV )  =  V_3 ( iV )  +  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
        FEM_1 ( iV )  =  V_1 ( iV )  -  CS ( iV )
        FEM_2 ( iV )  =  V_2 ( iV )  -  sqrt ( M_UU_22 ( iV ) ) * CS ( iV )
        FEM_3 ( iV )  =  V_3 ( iV )  -  sqrt ( M_UU_33 ( iV ) ) * CS ( iV )
      end do
      !$OMP end parallel do
    
    end if

  end procedure ComputeEigenspeedsFluidKernel



  module procedure ComputeCenterSpeedKernel 

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      D_Numerator, &
      S_Numerator, &
      D_Numerator_Inv
    logical ( KDL ) :: &
      UseDevice 

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( AC_I )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) &
      !$OMP& private ( iV, D_Numerator, S_Numerator, D_Numerator_Inv )
      do iV = 1, nValues

        D_Numerator &
          =  AP_I ( iV ) * D_IR ( iV )  +  AM_I ( iV ) * D_IL ( iV ) &
             -  F_D_IR ( iV )  +  F_D_IL ( iV )

        S_Numerator &
          =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
             -  F_S_IR ( iV )  +  F_S_IL ( iV )

        D_Numerator_Inv  &
          =  max ( D_Numerator, 0.0_KDR )  &
             /  max ( D_Numerator ** 2, tiny ( 0.0_KDR ) )

        AC_I ( iV )  &
          =  M_UU ( iV )  *  S_Numerator  *  D_Numerator_Inv

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& private ( iV, D_Numerator, S_Numerator, D_Numerator_Inv )
      do iV = 1, nValues

        D_Numerator &
          =  AP_I ( iV ) * D_IR ( iV )  +  AM_I ( iV ) * D_IL ( iV ) &
             -  F_D_IR ( iV )  +  F_D_IL ( iV )

        S_Numerator &
          =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
             -  F_S_IR ( iV )  +  F_S_IL ( iV )

        D_Numerator_Inv  &
          =  max ( D_Numerator, 0.0_KDR )  &
             /  max ( D_Numerator ** 2, tiny ( 0.0_KDR ) )

        AC_I ( iV )  &
          =  M_UU ( iV )  *  S_Numerator  *  D_Numerator_Inv

      end do !-- iV
      !$OMP  end parallel do
      
    end if

  end procedure ComputeCenterSpeedKernel


  module procedure ComputeCenterStatesKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      AM_VL, &
      AM_AC, &
      AM_AC_Inv, &
      AP_VR, &
      AP_AC, &
      AP_AC_Inv, &
      SqrtTiny
    real ( KDR ), dimension ( : ), pointer :: &
      V_Dim_ICL, &
      V_Dim_ICR, &
      S_Dim_ICL, &
      S_Dim_ICR
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nValues = size ( AC_I )
    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )
    
    select case ( iDim )
    case ( 1 )
      V_Dim_ICL => V_1_ICL
      V_Dim_ICR => V_1_ICR
      S_Dim_ICL => S_1_ICL
      S_Dim_ICR => S_1_ICR
    case ( 2 )
      V_Dim_ICL => V_2_ICL
      V_Dim_ICR => V_2_ICR
      S_Dim_ICL => S_2_ICL
      S_Dim_ICR => S_2_ICR
    case ( 3 )
      V_Dim_ICL => V_3_ICL
      V_Dim_ICR => V_3_ICR
      S_Dim_ICL => S_3_ICL
      S_Dim_ICR => S_3_ICR
    end select 

    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        V_Dim_ICL ( iV )  =  AC_I ( iV )
        V_Dim_ICR ( iV )  =  AC_I ( iV )

        AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
        D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

        S_1_ICL ( iV )  =  D_ICL ( iV )  *  V_1_ICL ( iV )
        S_1_ICR ( iV )  =  D_ICR ( iV )  *  V_1_ICR ( iV )

        S_2_ICL ( iV )  =  M_DD_22 ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
        S_2_ICR ( iV )  =  M_DD_22 ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

        S_3_ICL ( iV )  =  M_DD_33 ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
        S_3_ICR ( iV )  =  M_DD_33 ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

        P_ICL ( iV )  =  P_IL ( iV )  +  S_Dim_IL  ( iV ) * AM_VL &
                                      -  S_Dim_ICL ( iV ) * AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_Dim_IR  ( iV ) * AP_VR &
                                      +  S_Dim_ICR ( iV ) * AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV ) * AM_VL &
                           +  V_Dim_IL ( iV ) * P_IL ( iV ) &
                           -  AC_I ( iV ) * P_ICL ( iV ) ) &
                         * AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV ) * AP_VR &
                           -  V_Dim_IR ( iV ) * P_IR ( iV ) &
                           +  AC_I ( iV ) * P_ICR ( iV ) ) &
                         * AP_AC_Inv

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        V_Dim_ICL ( iV )  =  AC_I ( iV )
        V_Dim_ICR ( iV )  =  AC_I ( iV )

        AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
        D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

        S_1_ICL ( iV )  =  D_ICL ( iV )  *  V_1_ICL ( iV )
        S_1_ICR ( iV )  =  D_ICR ( iV )  *  V_1_ICR ( iV )

        S_2_ICL ( iV )  =  M_DD_22 ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
        S_2_ICR ( iV )  =  M_DD_22 ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

        S_3_ICL ( iV )  =  M_DD_33 ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
        S_3_ICR ( iV )  =  M_DD_33 ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

        P_ICL ( iV )  =  P_IL ( iV )  +  S_Dim_IL  ( iV ) * AM_VL &
                                      -  S_Dim_ICL ( iV ) * AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_Dim_IR  ( iV ) * AP_VR &
                                      +  S_Dim_ICR ( iV ) * AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV ) * AM_VL &
                           +  V_Dim_IL ( iV ) * P_IL ( iV ) &
                           -  AC_I ( iV ) * P_ICL ( iV ) ) &
                         * AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV ) * AP_VR &
                           -  V_Dim_IR ( iV ) * P_IR ( iV ) &
                           +  AC_I ( iV ) * P_ICR ( iV ) ) &
                         * AP_AC_Inv

      end do !-- iV
      !$OMP end parallel do
      
    end if

  end procedure ComputeCenterStatesKernel


  module procedure ComputeFluxes_HLLC_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( F_I )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        !-- If flagged for diffusive flux, leave HLL flux in place
        if ( DF_I ( iV ) > 0.0_KDR ) &
          cycle

        if ( AP_I ( iV ) /= 0.0_KDR .and. AM_I ( iV ) /= 0.0_KDR ) then
          !-- Use the appropriate center state flux
          if ( AC_I ( iV ) >= 0.0_KDR ) then
            F_I ( iV )  =  F_ICL ( iV )
          else !-- AC < 0
            F_I ( iV )  =  F_ICR ( iV )
          end if !-- AC >= 0
        else  !-- AP or AM == 0
          !-- Leave HLL flux in place (which is upwind in this case)
        end if !-- AP and AM /= 0

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP parallel do private ( iV ) 
      do iV = 1, nValues

        !-- If flagged for diffusive flux, leave HLL flux in place
        if ( DF_I ( iV ) > 0.0_KDR ) &
          cycle

        if ( AP_I ( iV ) /= 0.0_KDR .and. AM_I ( iV ) /= 0.0_KDR ) then
          !-- Use the appropriate center state flux
          if ( AC_I ( iV ) >= 0.0_KDR ) then
            F_I ( iV )  =  F_ICL ( iV )
          else !-- AC < 0
            F_I ( iV )  =  F_ICR ( iV )
          end if !-- AC >= 0
        else  !-- AP or AM == 0
          !-- Leave HLL flux in place (which is upwind in this case)
        end if !-- AP and AM /= 0

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure ComputeFluxes_HLLC_Kernel


  module procedure ComputeRawFluxesTemplate_P_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( F_S_Dim )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, nValues
        F_S_Dim ( iV ) = F_S_Dim ( iV )  +  P ( iV )
        F_G     ( iV ) =     ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        F_S_Dim ( iV ) = F_S_Dim ( iV )  +  P ( iV )
        F_G     ( iV ) =     ( G ( iV )  +  P ( iV ) ) * V_Dim ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeRawFluxesTemplate_P_Kernel


end submodule Fluid_P__Kernel
