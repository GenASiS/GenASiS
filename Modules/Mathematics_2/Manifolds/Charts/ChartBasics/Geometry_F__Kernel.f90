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


  module procedure Compute_FV_C_Kernel

    !-- Compute_FiniteVolume_Cylindrical_Kernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Pi, &
      R_I, R_O, &
      dZ, dPh

    Pi  =  CONSTANT % PI

    !$OMP parallel do private ( R_I, R_O, dZ, dPh )
    do iV = oV + 1, oV + nV

        dZ  =  W_2 ( iV )
       dPh  =  W_3 ( iV )

      R_I  =  RP_I ( iV )
      R_O  =  RP_I ( iV )  +  W_1 ( iV )

      select case ( nD )
      case ( 1 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I  
        A_I_2 ( iV )  =  Pi  *  ( R_O ** 2  -  R_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( R_O ** 2  -  R_I ** 2 )
        V ( iV )      =  Pi  *  ( R_O ** 2  -  R_I ** 2 )  
      case ( 2 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I  *  dZ  
        A_I_2 ( iV )  =  Pi  *  ( R_O ** 2  -  R_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( R_O ** 2  -  R_I ** 2 ) * dZ
        V ( iV )      =  Pi  *  ( R_O ** 2  -  R_I ** 2 )  *  dZ
      case ( 3 )
        A_I_1 ( iV )  =  R_I * dZ * dPh  
        A_I_2 ( iV )  =  0.5_KDR  *  ( R_O ** 2  -  R_I ** 2 )  *  dPh
        A_I_3 ( iV )  =  0.5_KDR  *  ( R_O ** 2  -  R_I ** 2 )  *  dZ
        V ( iV )      =  0.5_KDR  *  ( R_O ** 2  -  R_I ** 2 )  *  dZ * dPh
      end select

    end do
    !$OMP end parallel do

  end procedure Compute_FV_C_Kernel


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


  module procedure Compute_M_C_Kernel

    !-- Compute_Metric_Cylindrical_Kernel

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = oV + 1, oV + nV

        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 

        M_UU_11 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        if ( abs ( RP ( iV ) )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = oV + 1, oV + nV

        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 

        M_UU_11 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        if ( abs ( RP ( iV ) )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end parallel do

    end if

  end procedure Compute_M_C_Kernel


end submodule Geometry_F__Kernel
