#include "Preprocessor"

submodule ( Geometry_F__Form ) Geometry_F__Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure Compute_FV_R_Kernel

    !-- Compute_FiniteVolume_Rectangular_Kernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV     !-- nValues
    real ( KDR ) :: &
      dX, dY, dZ, &
      X_I, X_O, &
      Y_I, Y_O, &
      Z_I, Z_O

    nV  =  size ( V )

    !$OMP parallel do private ( dX, dY, dZ, X_I, X_O, Y_I, Y_O, Z_I, Z_O )
    do iV = 1, nV

      dX  =  W_1 ( iV )
      dY  =  W_2 ( iV )
      dZ  =  W_3 ( iV )

      X_I  =  E_I_1 ( iV )
      X_O  =  E_I_1 ( iV )  +  W_1 ( iV )

      Y_I  =  E_I_2 ( iV )
      Y_O  =  E_I_2 ( iV )  +  W_2 ( iV )

      Z_I  =  E_I_3 ( iV )
      Z_O  =  E_I_3 ( iV )  +  W_3 ( iV )

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

      AV_1_1 ( iV )  =  ( X_I  +  X_O )  /  2.0_KDR
!      AV_2_1 ( iV )  =  ( X_I ** 2  +  X_I * X_O  +  X_O ** 2 )  /  3.0_KDR
      AV_2_1 ( iV )  =  ( ( X_I  +  X_O )  /  2.0_KDR ) ** 2

      AV_1_2 ( iV )  =  ( Y_I  +  Y_O )  /  2.0_KDR
!      AV_2_2 ( iV )  =  ( Y_I ** 2  +  Y_I * Y_O  +  Y_O ** 2 )  /  3.0_KDR
      AV_2_2 ( iV )  =  ( ( Y_I  +  Y_O )  /  2.0_KDR ) ** 2

      AV_1_3 ( iV )  =  ( Z_I  +  Z_O )  /  2.0_KDR
!      AV_2_3 ( iV )  =  ( Z_I ** 2  +  Z_I * Z_O  +  Z_O ** 2 )  /  3.0_KDR
      AV_2_3 ( iV )  =  ( ( Z_I  +  Z_O )  /  2.0_KDR ) ** 2

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
      dZ, dPh, &
      RP_I, RP_O, &
       Z_I,  Z_O, &
      Ph_I, Ph_O

    Pi  =  CONSTANT % PI

    nV  =  size ( V )

    !$OMP parallel do private ( dZ, dPh, RP_I, RP_O, Z_I, Z_O, Ph_I, Ph_O )
    do iV = 1, nV

        dZ  =  W_2 ( iV )
       dPh  =  W_3 ( iV )

      RP_I  =  E_I_1 ( iV )
      RP_O  =  E_I_1 ( iV )  +  W_1 ( iV )

      Z_I  =  E_I_2 ( iV )
      Z_O  =  E_I_2 ( iV )  +  W_2 ( iV )

      Ph_I  =  E_I_3 ( iV )
      Ph_O  =  E_I_3 ( iV )  +  W_3 ( iV )

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

      AV_1_1 ( iV )  =  2.0_KDR  *  ( RP_O ** 3  -  RP_I ** 3 )  &
                        /  ( 3.0_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) )
!      AV_2_1 ( iV )  =  2.0_KDR  *  ( RP_O ** 4  -  RP_I ** 4 )  &
!                        /  ( 4.0_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) )
      AV_2_1 ( iV )  =  ( 2.0_KDR  *  ( RP_O ** 3  -  RP_I ** 3 )  &
                          /  ( 3.0_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) ) ) ** 2

      AV_1_2 ( iV )  =    ( Z_I  +  Z_O )  /  2.0_KDR
!      AV_2_2 ( iV )  =  ( Z_I ** 2  +  Z_I * Z_O  +  Z_O ** 2 )  /  3.0_KDR
      AV_2_2 ( iV )  =  ( ( Z_I  +  Z_O )  /  2.0_KDR ) ** 2

      AV_1_3 ( iV )  =    ( Ph_I  +  Ph_O )  /  2.0_KDR
!      AV_2_3 ( iV )  =  ( Ph_I ** 2  +  Ph_I * Ph_O  +  Ph_O ** 2 )  /  3.0_KDR
      AV_2_3 ( iV )  =  ( ( Ph_I  +  Ph_O )  /  2.0_KDR ) ** 2

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
       R_I,  R_O, &
      Th_I, Th_O, &
      Ph_I, Ph_O, &
       S_I,  S_O, &
       C_I,  C_O, &
      dPh

    Pi  =  CONSTANT % PI

    nV  =  size ( V )

    !$OMP parallel do private ( R_I, R_O, Th_I, Th_O, Ph_I, Ph_O, &
    !$OMP                       S_I, S_O, C_I, C_O, dPh )
    do iV = 1, nV

      R_I  =  E_I_1 ( iV )
      R_O  =  E_I_1 ( iV )  +  W_1 ( iV )

      Th_I  =  E_I_2 ( iV )
      Th_O  =  E_I_2 ( iV )  +  W_2 ( iV )

      Ph_I  =  E_I_3 ( iV )
      Ph_O  =  E_I_3 ( iV )  +  W_3 ( iV )

      S_I  =  sin ( Th_I )
      S_O  =  sin ( Th_O )

      C_I  =  cos ( Th_I )
      C_O  =  cos ( Th_O )

      dPh  =  W_3 ( iV )

      select case ( nD )
      case ( 1 )
        A_I_1 ( iV )  =  4.0_KDR  *  Pi  *  R_I ** 2
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_3 ( iV )  =  2.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )
        V ( iV )      =  4.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
      case ( 2 )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I ** 2  *  ( C_I  -  C_O )
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  S_I
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( C_I  -  C_O )
        V ( iV )      =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( C_I  -  C_O )
      case ( 3 )
        A_I_1 ( iV )  =  R_I ** 2  *  ( C_I  -  C_O )  * dPh 
        A_I_2 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  S_I  *  dPh
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( C_I  -  C_O )
        V ( iV )      =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( C_I  -  C_O )  *  dPh
      end select

!      AV_1_1 ( iV )  =    ( R_I  +  R_O )  /  2.0_KDR
!      AV_2_1 ( iV )  =  ( ( R_I  +  R_O )  /  2.0_KDR ) ** 2

      AV_1_1 ( iV )  =  3.0_KDR  *  ( R_O ** 4  -  R_I ** 4 )  &
                        /  ( 4.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) )
!      AV_2_1 ( iV )  =  3.0_KDR  *  ( R_O ** 5  -  R_I ** 5 )  &
!                        /  ( 5.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) )
      AV_2_1 ( iV )  =  ( 3.0_KDR  *  ( R_O ** 4  -  R_I ** 4 )  &
                          /  ( 4.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) ) ) ** 2

      if ( nD  >  1 ) then

!        AV_1_2 ( iV )  =    ( Th_I  +  Th_O )  /  2.0_KDR
!        AV_2_2 ( iV )  =  ( ( Th_I  +  Th_O )  /  2.0_KDR ) ** 2

        AV_1_2 ( iV )  =  ( S_O  -  S_I  +  Th_I * C_I  -  Th_O * C_O )  &
                          /  ( C_I  -  C_O )
        ! AV_2_2 ( iV )  =  ( 2.0_KDR * S_O  -  2.0_KDR * S_I  &
        !                     +  ( Th_I ** 2  -  2.0_KDR ) * C_I  &
        !                     -  ( Th_O ** 2  -  2.0_KDR ) * C_O )  &
        !                   /  ( C_I  -  C_O )
        AV_2_2 ( iV )  =  ( ( S_O  -  S_I  +  Th_I * C_I  -  Th_O * C_O )  &
                            /  ( C_I  -  C_O ) ) ** 2

      end if

      if ( nD  >  2 ) then

        AV_1_3 ( iV )  =    ( Ph_I  +  Ph_O )  /  2.0_KDR
!        AV_2_3 ( iV )  =  ( Ph_I ** 2  +  Ph_I * Ph_O  +  Ph_O ** 2 )  &
!                          /  3.0_KDR
        AV_2_3 ( iV )  =  ( ( Ph_I  +  Ph_O )  /  2.0_KDR ) ** 2

      end if

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


end submodule Geometry_F__Kernel
