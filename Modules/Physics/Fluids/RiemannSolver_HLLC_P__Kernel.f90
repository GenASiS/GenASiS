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


end submodule RiemannSolver_HLLC_P__Kernel
