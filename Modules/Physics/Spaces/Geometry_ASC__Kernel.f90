#include "Preprocessor"

submodule ( Geometry_ASC__Form ) Geometry_ASC__Kernel

  use Basics
  implicit none
  
contains


  module procedure ComputeGravityUniformKernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    logical ( KDL ) :: &
      UseDevice
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( Phi )

    select case ( nDimensions )
    case ( 1 )
    
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  X ( iV )
          GradPhi_1 ( iV )  =  Acceleration
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST )
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  X ( iV )
          GradPhi_1 ( iV )  =  Acceleration
        end do
        !$OMP  end parallel do
      end if

    case ( 2 )
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  Y ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  Acceleration
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST )
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  Y ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  Acceleration
        end do
        !$OMP  end parallel do
      end if

    case ( 3 )
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  Z ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  0.0_KDR
          GradPhi_3 ( iV )  =  Acceleration
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST )        
        do iV = 1, nValues
          Phi ( iV )        =  Acceleration  *  Z ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  0.0_KDR
          GradPhi_3 ( iV )  =  Acceleration
        end do
        !$OMP  end parallel do
      end if

    end select !-- nDimensions

  end procedure ComputeGravityUniformKernel


  module procedure ComputeGravityCentralMassKernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( Phi )
    
    if ( UseDevice ) then
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues
        Phi ( iV )        =  - G * M  /  X_1 ( iV )
        GradPhi_1 ( iV )  =  G * M  /  X_1 ( iV ) ** 2
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues
        Phi ( iV )        =  - G * M  /  X_1 ( iV )
        GradPhi_1 ( iV )  =  G * M  /  X_1 ( iV ) ** 2
      end do
      !$OMP  end parallel do
    
    end if

  end procedure ComputeGravityCentralMassKernel


  module procedure ComputeGravitySource

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      FourPi_G
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    FourPi_G  =  4.0_KDR  *  CONSTANT % PI  *  G

    nValues = size ( S )
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& shared ( FourPI_G )
      do iV = 1, nValues
        S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& shared ( FourPI_G )
      do iV = 1, nValues
        S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
      end do
      !$OMP  end parallel do
    
    end if

  end procedure ComputeGravitySource


end submodule Geometry_ASC__Kernel
