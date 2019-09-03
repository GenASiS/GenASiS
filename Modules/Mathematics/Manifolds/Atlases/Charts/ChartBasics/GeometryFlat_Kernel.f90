#include "Preprocessor"

submodule ( GeometryFlat_Form ) GeometryFlat_Kernel

  use Basics
  
  implicit none
  
contains 


  module procedure SetFiniteVolumeRectangularKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      dX, dY, dZ

    !$OMP parallel do private ( iV, dX, dY, dZ )
    do iV = oValue + 1, oValue + nValues

      dX  =  W_L_1 ( iV )  +  W_R_1 ( iV )
      dY  =  W_L_2 ( iV )  +  W_R_2 ( iV )
      dZ  =  W_L_3 ( iV )  +  W_R_3 ( iV )

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  dX
        A_I_1 ( iV )  =  1.0_KDR
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX
      case ( 2 )
        V ( iV )      =  dX * dY
        A_I_1 ( iV )  =  dY
        A_I_2 ( iV )  =  dX
        A_I_3 ( iV )  =  dX * dY
      case ( 3 )
        V ( iV )      =  dX * dY * dZ
        A_I_1 ( iV )  =  dY * dZ
        A_I_2 ( iV )  =  dZ * dX
        A_I_3 ( iV )  =  dX * dY
      end select

    end do
    !$OMP end parallel do

  end procedure SetFiniteVolumeRectangularKernel


  module procedure SetFiniteVolumeCylindricalKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Pi, &
      RP_I, RP_O, &
      dZ, dPh

    Pi  =  CONSTANT % PI

    !$OMP parallel do private ( iV, RP_I, RP_O, dZ, dPh )
    do iV = oValue + 1, oValue + nValues

      RP_I  =  RP_C ( iV )  -  W_L_1 ( iV )
      RP_O  =  RP_C ( iV )  +  W_R_1 ( iV )

        dZ  =  W_L_2 ( iV )  +  W_R_2 ( iV ) 
       dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )
      case ( 2 )
        V ( iV )      =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  RP_I  *  dZ  
        A_I_2 ( iV )  =  Pi  *  ( RP_O ** 2  -  RP_I ** 2 )
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) * dZ
      case ( 3 )
        V ( iV )      =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dZ * dPh
        A_I_1 ( iV )  =  RP_I * dZ * dPh  
        A_I_2 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 )  *  dPh
        A_I_3 ( iV )  =  0.5_KDR  *  ( RP_O ** 2  -  RP_I ** 2 ) * dZ
      end select

    end do
    !$OMP end parallel do

  end procedure SetFiniteVolumeCylindricalKernel


  module procedure SetFiniteVolumeSphericalKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Pi, &
      R_I, R_O, &
      Th_I, Th_O, &
      dPh

    Pi  =  CONSTANT % PI

    !$OMP parallel do private ( iV, R_I, R_O, Th_I, Th_O, dPh )
    do iV = oValue + 1, oValue + nValues

      R_I  =  R_C ( iV )  -  W_L_1 ( iV )
      R_O  =  R_C ( iV )  +  W_R_1 ( iV )

      Th_I  =  Th_C ( iV )  -  W_L_2 ( iV )
      Th_O  =  Th_C ( iV )  +  W_R_2 ( iV )

      dPh  =  W_L_3 ( iV )  +  W_R_3 ( iV ) 

      select case ( nDimensions )
      case ( 1 )
        V ( iV )      =  4.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_1 ( iV )  =  4.0_KDR  *  Pi  *  R_I ** 2
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 )
        A_I_3 ( iV )  =  2.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )
      case ( 2 )
        V ( iV )      =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        A_I_1 ( iV )  =  2.0_KDR  *  Pi  *  R_I ** 2  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
        A_I_2 ( iV )  =  2.0_KDR / 3.0_KDR * Pi *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  sin ( Th_I )
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
      case ( 3 )
        V ( iV )      =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 ) &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )  *  dPh
        A_I_1 ( iV )  =  R_I ** 2  *  ( cos ( Th_I )  -  cos ( Th_O ) )  * dPh 
        A_I_2 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  sin ( Th_I )  *  dPh
        A_I_3 ( iV )  =  1.0_KDR / 3.0_KDR  *  ( R_O ** 3  -  R_I ** 3 )  &
                         *  ( cos ( Th_I )  -  cos ( Th_O ) )
      end select

    end do
    !$OMP end parallel do

  end procedure SetFiniteVolumeSphericalKernel


  module procedure SetMetricRectangularKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  1.0_KDR
        M_UU_22 ( iV )  =  1.0_KDR
        M_UU_33 ( iV )  =  1.0_KDR
      end do
      !$OMP end parallel do

    end if

  end procedure SetMetricRectangularKernel


  module procedure SetMetricCylindricalKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 
        M_UU_22 ( iV )  =  1.0_KDR
        if ( RP ( iV )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = oValue + 1, oValue + nValues
        M_DD_22 ( iV )  =  1.0_KDR
        M_DD_33 ( iV )  =  RP ( iV ) ** 2 
        M_UU_22 ( iV )  =  1.0_KDR
        if ( RP ( iV )  >  0.0_KDR ) then
          M_UU_33 ( iV )  =  1.0_KDR  /  RP ( iV ) ** 2
        else
          M_UU_33 ( iV )  =  0.0_KDR
        end if
      end do
      !$OMP end parallel do

    end if

  end procedure SetMetricCylindricalKernel


  module procedure SetMetricSphericalKernel

    integer ( KDI ) :: &
      iV  !-- iValue
    real ( KDR ) :: &
      Sin_Th
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then      
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, Sin_Th )
      do iV = oValue + 1, oValue + nValues

        select case ( nDimensions )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2
        if ( R ( iV )  *  Sin_Th  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV, Sin_Th )
      do iV = oValue + 1, oValue + nValues

        select case ( nDimensions )
        case ( 1 )
          Sin_Th  =  1.0_KDR
        case ( 2 )
          Sin_Th  =  sin ( Th ( iV ) )
        case ( 3 )
          Sin_Th  =  sin ( Th ( iV ) )
        end select

        M_DD_22 ( iV )  =  R ( iV ) ** 2
        M_DD_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** 2
        if ( R ( iV )  *  Sin_Th  >  0.0_KDR ) then
          M_UU_22 ( iV )  =  R ( iV ) ** ( -2 )
          M_UU_33 ( iV )  =  ( R ( iV )  *  Sin_Th ) ** ( -2 )
        else
          M_UU_22 ( iV )  =  0.0_KDR
          M_UU_33 ( iV )  =  0.0_KDR
        end if

      end do
      !$OMP end parallel do
      
    end if

  end procedure SetMetricSphericalKernel
  

  module procedure ComputeEdgesKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV = size ( X )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do private ( iV ) &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, nV
        X_I ( iV )  =  X ( iV )  -  dX_L ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else      
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nV
        X_I ( iV )  =  X ( iV )  -  dX_L ( iV )
      end do
      !$OMP end parallel do
    end if
    
  end procedure ComputeEdgesKernel


end submodule GeometryFlat_Kernel
