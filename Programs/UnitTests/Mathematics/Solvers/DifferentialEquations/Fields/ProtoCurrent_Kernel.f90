#include "Preprocessor"

submodule ( ProtoCurrent_Form ) ProtoCurrent_Kernel
  
  use Basics
  
  implicit none
  
contains

  module procedure ComputeConservedDensityKernel
  
    integer ( KDI ) :: &
      iV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
     
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV )
      do iV = 1, size ( N )
        D ( iV ) = N ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV )
      do iV = 1, size ( N )
        D ( iV ) = N ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeConservedDensityKernel


  module procedure ComputeComovingDensityKernel

    integer ( KDI ) :: &
      iV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
     
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV )
      do iV = 1, size ( N )
        N ( iV ) = D ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV )
      do iV = 1, size ( N )
        N ( iV ) = D ( iV )
      end do
      !$OMP end parallel do
    end if


  end procedure ComputeComovingDensityKernel


  module procedure ComputeEigenspeedsKernel

    integer ( KDI ) :: &
      iV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV )
      do iV = 1, size ( FEP_1 )
        FEP_1 ( iV ) = V_1 ( iV )
        FEP_2 ( iV ) = V_2 ( iV )
        FEP_3 ( iV ) = V_3 ( iV )
        FEM_1 ( iV ) = V_1 ( iV )
        FEM_2 ( iV ) = V_2 ( iV )
        FEM_3 ( iV ) = V_3 ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV )
      do iV = 1, size ( FEP_1 )
        FEP_1 ( iV ) = V_1 ( iV )
        FEP_2 ( iV ) = V_2 ( iV )
        FEP_3 ( iV ) = V_3 ( iV )
        FEM_1 ( iV ) = V_1 ( iV )
        FEM_2 ( iV ) = V_2 ( iV )
        FEM_3 ( iV ) = V_3 ( iV )
      end do
      !$OMP end parallel do
    end if
    
  end procedure ComputeEigenspeedsKernel


  module procedure ComputeVelocityKernel
    
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Abs_K
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    Abs_K = sqrt ( dot_product ( K, K ) )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV )
      do iV = 1, size ( V_1 )  
        V_1 ( iV ) = V * K ( 1 ) / Abs_K
        V_2 ( iV ) = V * K ( 2 ) / Abs_K
        V_3 ( iV ) = V * K ( 3 ) / Abs_K
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV )
      do iV = 1, size ( V_1 )  
        V_1 ( iV ) = V * K ( 1 ) / Abs_K
        V_2 ( iV ) = V * K ( 2 ) / Abs_K
        V_3 ( iV ) = V * K ( 3 ) / Abs_K
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeVelocityKernel


  module procedure ComputeRawFluxesKernel
  
    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      Abs_K
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV )
      do iV = 1, size ( F_D )
        F_D ( iV )  =  D ( iV ) * V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP& private ( iV )
      do iV = 1, size ( F_D )
        F_D ( iV )  =  D ( iV ) * V_Dim ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeRawFluxesKernel


end submodule ProtoCurrent_Kernel
