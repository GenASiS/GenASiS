#include "Preprocessor"

submodule ( Geometry_F_C__Form ) Geometry_F_C__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure Compute_FV_R_Kernel

    !-- Compute_FiniteVolume_Rectangular_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues
    real ( KDR ) :: &
      dX, dY, dZ

    nV  =  size ( V )

    !$OMP parallel do private ( dX, dY, dZ )
    do iV = 1, nV

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
      iV, &  !-- iValue
      nV     !-- nValues
    real ( KDR ) :: &
      Pi, &
      RP_I, RP_O, &
      dZ, dPh

    Pi  =  CONSTANT % PI

    nV  =  size ( V )

    !$OMP parallel do private ( RP_I, RP_O, dZ, dPh )
    do iV = 1, nV

      RP_I  =  E_I_1 ( iV )
      RP_O  =  E_I_1 ( iV )  +  W_1 ( iV )

        dZ  =  W_2 ( iV )
       dPh  =  W_3 ( iV )

      select case ( nD )
      case ( 1 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  
      case ( 2 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  *  dZ  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) * dZ
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ
      case ( 3 )
        A_I_1 ( iV )  =  RP_I * dZ * dPh  
        A_I_2 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dPh
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ
        V ( iV )      =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ * dPh
      end select

    end do
    !$OMP end parallel do

  end procedure Compute_FV_C_Kernel


  module procedure Compute_FV_S_Kernel

    !-- Compute_FiniteVolume_Spherical_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues
    real ( KDR ) :: &
      Pi, &
      R_I, R_O, &
      Th_I, Th_O, &
      dPh

    Pi  =  CONSTANT % PI

    nV  =  size ( V )

    !$OMP parallel do private ( R_I, R_O, Th_I, Th_O, dPh )
    do iV = 1, nV

      R_I  =  E_I_1 ( iV )
      R_O  =  E_I_1 ( iV )  +  W_1 ( iV )

      Th_I  =  E_I_2 ( iV )
      Th_O  =  E_I_2 ( iV )  +  W_2 ( iV )

      dPh  =  W_3 ( iV )

      select case ( nD )
      case ( 1 )
        A_I_1 ( iV )  =  4.0_KDR  *  Pi  *  R_I ** 2
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_3 ( iV )  =  2.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )
        V ( iV )      =  4.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
      case ( 2 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I ** 2  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  sin ( Th_I )
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        V ( iV )      =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
      case ( 3 )
        A_I_1 ( iV )  =  R_I ** 2  *  ( cos ( Th_I )  -  cos ( Th_O ) )  * dPh 
        A_I_2 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  sin ( Th_I )  *  dPh
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        V ( iV )      =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )  *  dPh
      end select

    end do
    !$OMP end parallel do

  end procedure Compute_FV_S_Kernel


  module procedure Compute_M_R_Kernel

    !-- Compute_Metric_Rectangular_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV  =  size ( M_DD_11 )

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
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
      do iV = 1, nV
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
      iV, &  !-- iValue
      nV     !-- nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV  =  size ( M_DD_11 )

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV

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

      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

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


  module procedure Compute_M_S_Kernel

    !-- Compute_Metric_Spherical_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues
    real ( KDR ) :: &
      Sin_Th
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV  =  size ( M_DD_11 )

    if ( UseDevice ) then      
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( Sin_Th )
      do iV = 1, nV

        select case ( nD )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2

        M_UU_11 ( iV )  =  1.0_KDR
        if ( abs ( R ( iV )  *  Sin_Th )  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST ) private ( Sin_Th )
      do iV = 1, nV

        select case ( nD )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_11 ( iV )  =  1.0_KDR
        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2

        M_UU_11 ( iV )  =  1.0_KDR
        if ( abs ( R ( iV )  *  Sin_Th )  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end parallel do
      
    end if

  end procedure Compute_M_S_Kernel


end submodule Geometry_F_C__Kernel
