#include "Preprocessor"

submodule ( Geometry_F__Form ) Geometry_F__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure Compute_FV_R_Kernel

    !-- Compute_FiniteVolume_Rectangular_Kernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      dX, dY, dZ

    !$OMP parallel do private ( dX, dY, dZ )
    do iV = oV + 1, oV + nV

      dX  =  W_1 ( iV )
      dY  =  W_2 ( iV )
      dZ  =  W_3 ( iV )

      select case ( nD )
      case ( 1 )
        A_I_1 ( iV )  =  1.0_KDR
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX
            V ( iV )  =  dX
      case ( 2 )
        A_I_1 ( iV )  =  dY
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX * dY
            V ( iV )  =  dX * dY
      case ( 3 )
        A_I_1 ( iV )  =  dY * dZ
        A_I_2 ( iV )  =  dZ * dX
        A_I_3 ( iV )  =  dX * dY
            V ( iV )  =  dX * dY * dZ
      end select

    end do
    !$OMP end parallel do

  end procedure Compute_FV_R_Kernel


  module procedure Compute_M_R_Kernel

    !-- Compute_Metric_Rectangular_Kernel

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = oV + 1, oV + nV
        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_11 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST ) 
      do iV = oV + 1, oV + nV
        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_11 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP end parallel do

    end if

  end procedure Compute_M_R_Kernel


end submodule Geometry_F__Kernel
