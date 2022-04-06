#include "Preprocessor"

submodule ( RiemannSolver_HLLC_P__Form ) RiemannSolver_HLLC_P__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeCenterSpeedKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      MD_Numerator, &
      S_Numerator, &
      MD_Numerator_Inv
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
  
    nV  =  size ( AC_I )
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( MD_Numerator, S_Numerator, MD_Numerator_Inv )
      do iV = 1, nV

        MD_Numerator &
          =  AP_I ( iV ) * M_IR ( iV ) * D_IR ( iV )  &
             +  AM_I ( iV ) * M_IL ( iV ) * D_IL ( iV ) &
             -  F_D_IR ( iV )  +  F_D_IL ( iV )

        S_Numerator &
          =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
             -  F_S_IR ( iV )  +  F_S_IL ( iV )

        MD_Numerator_Inv  &
          =  max ( MD_Numerator, 0.0_KDR )  &
             /  max ( MD_Numerator ** 2, tiny ( 0.0_KDR ) )

        AC_I ( iV )  &
          =  M_UU ( iV )  *  S_Numerator  *  MD_Numerator_Inv

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( MD_Numerator, S_Numerator, MD_Numerator_Inv )
      do iV = 1, nV

        MD_Numerator &
          =  AP_I ( iV ) * M_IR ( iV ) * D_IR ( iV )  &
             +  AM_I ( iV ) * M_IL ( iV ) * D_IL ( iV ) &
             -  F_D_IR ( iV )  +  F_D_IL ( iV )

        S_Numerator &
          =  AP_I ( iV ) * S_IR ( iV )  +  AM_I ( iV ) * S_IL ( iV ) &
             -  F_S_IR ( iV )  +  F_S_IL ( iV )

        MD_Numerator_Inv  &
          =  max ( MD_Numerator, 0.0_KDR )  &
             /  max ( MD_Numerator ** 2, tiny ( 0.0_KDR ) )

        AC_I ( iV )  &
          =  M_UU ( iV )  *  S_Numerator  *  MD_Numerator_Inv

      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure ComputeCenterSpeedKernel


  module procedure ComputeCenterStatesKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      AM_VL, &
      AM_AC, &
      AM_AC_Inv, &
      AP_VR, &
      AP_AC, &
      AP_AC_Inv, &
      SqrtTiny
    logical ( KDL ) :: &  
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( AC_I )
    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
            
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( AM_VL, AM_AC, AM_AC_Inv, AP_VR, AP_AC, AP_AC_Inv ) &
      !$OMP firstprivate ( SqrtTiny )
      do iV  =  1,  nV

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        V_D_ICL ( iV )  =  AC_I ( iV )
        V_D_ICR ( iV )  =  AC_I ( iV )

        AM_VL     =  AM_I ( iV )  +  V_D_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_D_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        M_ICL ( iV )  =  M_IL ( iV )
        M_ICR ( iV )  =  M_IR ( iV )

        D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
        D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

        S_1_ICL ( iV )  =  M_DD_11 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_1_ICL ( iV )
        S_1_ICR ( iV )  =  M_DD_11 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_1_ICR ( iV )

        S_2_ICL ( iV )  =  M_DD_22 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
        S_2_ICR ( iV )  =  M_DD_22 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

        S_3_ICL ( iV )  =  M_DD_33 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
        S_3_ICR ( iV )  =  M_DD_33 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

        P_ICL ( iV )  =  P_IL ( iV )  +  S_D_IL  ( iV ) * AM_VL &
                                      -  S_D_ICL ( iV ) * AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_D_IR  ( iV ) * AP_VR &
                                      +  S_D_ICR ( iV ) * AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV ) * AM_VL &
                           +  V_D_IL ( iV ) * P_IL ( iV ) &
                           -  AC_I ( iV ) * P_ICL ( iV ) ) &
                         * AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV ) * AP_VR &
                           -  V_D_IR ( iV ) * P_IR ( iV ) &
                           +  AC_I ( iV ) * P_ICR ( iV ) ) &
                         * AP_AC_Inv

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else 
      
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( AM_VL, AM_AC, AM_AC_Inv, AP_VR, AP_AC, AP_AC_Inv ) &
      !$OMP firstprivate ( SqrtTiny )
      do iV  =  1,  nV

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        V_D_ICL ( iV )  =  AC_I ( iV )
        V_D_ICR ( iV )  =  AC_I ( iV )

        AM_VL     =  AM_I ( iV )  +  V_D_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_D_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        M_ICL ( iV )  =  M_IL ( iV )
        M_ICR ( iV )  =  M_IR ( iV )

        D_ICL ( iV )  =  D_IL ( iV ) * AM_VL * AM_AC_Inv
        D_ICR ( iV )  =  D_IR ( iV ) * AP_VR * AP_AC_Inv

        S_1_ICL ( iV )  =  M_DD_11 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_1_ICL ( iV )
        S_1_ICR ( iV )  =  M_DD_11 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_1_ICR ( iV )

        S_2_ICL ( iV )  =  M_DD_22 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_2_ICL ( iV )
        S_2_ICR ( iV )  =  M_DD_22 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_2_ICR ( iV )

        S_3_ICL ( iV )  =  M_DD_33 ( iV )  &
                           *  M_ICL ( iV )  *  D_ICL ( iV )  *  V_3_ICL ( iV )
        S_3_ICR ( iV )  =  M_DD_33 ( iV )  &
                           *  M_ICR ( iV )  *  D_ICR ( iV )  *  V_3_ICR ( iV )

        P_ICL ( iV )  =  P_IL ( iV )  +  S_D_IL  ( iV ) * AM_VL &
                                      -  S_D_ICL ( iV ) * AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_D_IR  ( iV ) * AP_VR &
                                      +  S_D_ICR ( iV ) * AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV ) * AM_VL &
                           +  V_D_IL ( iV ) * P_IL ( iV ) &
                           -  AC_I ( iV ) * P_ICL ( iV ) ) &
                         * AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV ) * AP_VR &
                           -  V_D_IR ( iV ) * P_IR ( iV ) &
                           +  AC_I ( iV ) * P_ICR ( iV ) ) &
                         * AP_AC_Inv

      end do !-- iV
      !$OMP end parallel do

    end if
    
  end procedure ComputeCenterStatesKernel


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      iF, &
      iF_F, &
      nV, &
      nF
    logical ( KDL ) :: &  
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV  =  size ( RSV, dim = 1 )
    nF  =  size ( iaFluxes )
    
    associate &
      ( F_I  => RSV, &
        AP_I => RSV ( :, iAP ), &
        AM_I => RSV ( :, iAM ), &
        AC_I => RSV ( :, iAC ) )
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iF_F )
      do iF  =  1,  nF
        do iV  =  1,  nV

        !   !-- If flagged for diffusive flux, leave HLL flux in place
        !   if ( DF_I ( iV ) > 0.0_KDR ) &
        !     cycle

          iF_F  =  iaFluxes ( iF )

          if ( AP_I ( iV )  /=  0.0_KDR .and. AM_I ( iV )  /=  0.0_KDR ) then
            !-- Use the appropriate center state flux
            if ( AC_I ( iV )  >=  0.0_KDR ) then
              F_I ( iV, iF_F )  =  F_ICL ( iV, iF_F )
            else !-- AC < 0
              F_I ( iV, iF_F )  =  F_ICR ( iV, iF_F )
            end if !-- AC >= 0
          else  !-- AP or AM == 0
            !-- Leave HLL flux in place (which is upwind in this case)
          end if !-- AP and AM /= 0

        end do !-- iV
      end do !-- iF
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iF_F )
      do iF  =  1,  nF
        do iV  =  1,  nV

        !   !-- If flagged for diffusive flux, leave HLL flux in place
        !   if ( DF_I ( iV ) > 0.0_KDR ) &
        !     cycle

          iF_F  =  iaFluxes ( iF )

          if ( AP_I ( iV )  /=  0.0_KDR .and. AM_I ( iV )  /=  0.0_KDR ) then
            !-- Use the appropriate center state flux
            if ( AC_I ( iV )  >=  0.0_KDR ) then
              F_I ( iV, iF_F )  =  F_ICL ( iV, iF_F )
            else !-- AC < 0
              F_I ( iV, iF_F )  =  F_ICR ( iV, iF_F )
            end if !-- AC >= 0
          else  !-- AP or AM == 0
            !-- Leave HLL flux in place (which is upwind in this case)
          end if !-- AP and AM /= 0

        end do !-- iV
      end do !-- iF
      !$OMP end parallel do
    
    end if

    end associate   !-- F_I, AP_I, AM_I

  end procedure ComputeKernel


end submodule RiemannSolver_HLLC_P__Kernel
