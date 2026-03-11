#include "Preprocessor"

submodule ( RiemannSolver_HLLC_P_HN__Form ) RiemannSolver_HLLC_P_HN__Kernel
  
  use Basics
  
  implicit none
  
contains


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
      V_D_IL,  V_D_IR, &
      S_D_IL,  S_D_IR, &
      S_D_ICL, S_D_ICR, &
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
      !$OMP private ( V_D_IL, V_D_IR, S_D_IL, S_D_IR, S_D_ICL, S_D_ICR ) &
      !$OMP firstprivate ( SqrtTiny )
      do iV  =  1,  nV

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        select case ( iD )
        case ( 1 )
          V_D_IL          =  V_1_IL ( iV )
          V_D_IR          =  V_1_IR ( iV )
          V_1_ICL ( iV )  =  AC_I ( iV )
          V_1_ICR ( iV )  =  AC_I ( iV )
        case ( 2 )
          V_D_IL          =  V_2_IL ( iV )
          V_D_IR          =  V_2_IR ( iV )
          V_2_ICL ( iV )  =  AC_I ( iV )
          V_2_ICR ( iV )  =  AC_I ( iV )
        case ( 3 )
          V_D_IL          =  V_3_IL ( iV )
          V_D_IR          =  V_3_IR ( iV )
          V_3_ICL ( iV )  =  AC_I ( iV )
          V_3_ICR ( iV )  =  AC_I ( iV )
        end select !-- iD

        AM_VL     =  AM_I ( iV )  +  V_D_IL
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_D_IR
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

        select case ( iD )
        case ( 1 )
          S_D_IL   =  S_1_IL ( iV )
          S_D_IR   =  S_1_IR ( iV )
          S_D_ICL  =  S_1_ICL ( iV )
          S_D_ICR  =  S_1_ICR ( iV )
        case ( 2 )
          S_D_IL   =  S_2_IL ( iV )
          S_D_IR   =  S_2_IR ( iV )
          S_D_ICL  =  S_2_ICL ( iV )
          S_D_ICR  =  S_2_ICR ( iV )
        case ( 3 )
          S_D_IL   =  S_3_IL ( iV )
          S_D_IR   =  S_3_IR ( iV )
          S_D_ICL  =  S_3_ICL ( iV )
          S_D_ICR  =  S_3_ICR ( iV )
        end select !-- iD

        P_ICL ( iV )  =  P_IL ( iV )  +  S_D_IL   *  AM_VL  &
                                      -  S_D_ICL  *  AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_D_IR   *  AP_VR  &
                                      +  S_D_ICR  *  AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV )  *  AM_VL &
                           +  V_D_IL  *  P_IL ( iV ) &
                           -  AC_I ( iV )  *  P_ICL ( iV ) ) &
                         *  AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV )  *  AP_VR &
                           -  V_D_IR  *  P_IR ( iV ) &
                           +  AC_I ( iV )  *  P_ICR ( iV ) ) &
                         *  AP_AC_Inv

        DS_ICL ( iV )  =  DS_IL ( iV ) * AM_VL * AM_AC_Inv
        DS_ICR ( iV )  =  DS_IR ( iV ) * AP_VR * AP_AC_Inv

        DE_ICL ( iV )  =  DE_IL ( iV ) * AM_VL * AM_AC_Inv
        DE_ICR ( iV )  =  DE_IR ( iV ) * AP_VR * AP_AC_Inv

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else 
      
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( AM_VL, AM_AC, AM_AC_Inv, AP_VR, AP_AC, AP_AC_Inv ) &
      !$OMP private ( V_D_IL, V_D_IR, S_D_IL, S_D_IR, S_D_ICL, S_D_ICR ) &
      !$OMP firstprivate ( SqrtTiny )
      do iV  =  1,  nV

        V_1_ICL ( iV )  =  V_1_IL ( iV )
        V_1_ICR ( iV )  =  V_1_IR ( iV )

        V_2_ICL ( iV )  =  V_2_IL ( iV )
        V_2_ICR ( iV )  =  V_2_IR ( iV )

        V_3_ICL ( iV )  =  V_3_IL ( iV )
        V_3_ICR ( iV )  =  V_3_IR ( iV )

        select case ( iD )
        case ( 1 )
          V_D_IL          =  V_1_IL ( iV )
          V_D_IR          =  V_1_IR ( iV )
          V_1_ICL ( iV )  =  AC_I ( iV )
          V_1_ICR ( iV )  =  AC_I ( iV )
        case ( 2 )
          V_D_IL          =  V_2_IL ( iV )
          V_D_IR          =  V_2_IR ( iV )
          V_2_ICL ( iV )  =  AC_I ( iV )
          V_2_ICR ( iV )  =  AC_I ( iV )
        case ( 3 )
          V_D_IL          =  V_3_IL ( iV )
          V_D_IR          =  V_3_IR ( iV )
          V_3_ICL ( iV )  =  AC_I ( iV )
          V_3_ICR ( iV )  =  AC_I ( iV )
        end select !-- iD

        AM_VL     =  AM_I ( iV )  +  V_D_IL
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_D_IR
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

        select case ( iD )
        case ( 1 )
          S_D_IL   =  S_1_IL ( iV )
          S_D_IR   =  S_1_IR ( iV )
          S_D_ICL  =  S_1_ICL ( iV )
          S_D_ICR  =  S_1_ICR ( iV )
        case ( 2 )
          S_D_IL   =  S_2_IL ( iV )
          S_D_IR   =  S_2_IR ( iV )
          S_D_ICL  =  S_2_ICL ( iV )
          S_D_ICR  =  S_2_ICR ( iV )
        case ( 3 )
          S_D_IL   =  S_3_IL ( iV )
          S_D_IR   =  S_3_IR ( iV )
          S_D_ICL  =  S_3_ICL ( iV )
          S_D_ICR  =  S_3_ICR ( iV )
        end select !-- iD

        P_ICL ( iV )  =  P_IL ( iV )  +  S_D_IL   *  AM_VL  &
                                      -  S_D_ICL  *  AM_AC
        P_ICR ( iV )  =  P_IR ( iV )  -  S_D_IR   *  AP_VR  &
                                      +  S_D_ICR  *  AP_AC

        G_ICL ( iV )  =  ( G_IL ( iV )  *  AM_VL &
                           +  V_D_IL  *  P_IL ( iV ) &
                           -  AC_I ( iV )  *  P_ICL ( iV ) ) &
                         *  AM_AC_Inv
        G_ICR ( iV )  =  ( G_IR ( iV )  *  AP_VR &
                           -  V_D_IR  *  P_IR ( iV ) &
                           +  AC_I ( iV )  *  P_ICR ( iV ) ) &
                         *  AP_AC_Inv

        DS_ICL ( iV )  =  DS_IL ( iV ) * AM_VL * AM_AC_Inv
        DS_ICR ( iV )  =  DS_IR ( iV ) * AP_VR * AP_AC_Inv

        DE_ICL ( iV )  =  DE_IL ( iV ) * AM_VL * AM_AC_Inv
        DE_ICR ( iV )  =  DE_IR ( iV ) * AP_VR * AP_AC_Inv

      end do !-- iV
      !$OMP end parallel do

    end if
    
  end procedure ComputeCenterStatesKernel


end submodule RiemannSolver_HLLC_P_HN__Kernel
