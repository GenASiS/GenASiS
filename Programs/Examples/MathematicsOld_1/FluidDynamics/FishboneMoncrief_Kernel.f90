#include "Preprocessor"

submodule ( FishboneMoncrief_Form ) FishboneMoncrief_Kernel

  use Basics

  implicit none

contains


  module procedure ApplySourcesKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ), dimension ( size ( N ) ) :: &
      F_1
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( N )
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nV
        F_1 ( iV )  =  - G * M * N ( iV )  /  R ( iV ) ** 2
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KV_M_1 ( iV )  =  KV_M_1 ( iV )  +  dT * F_1 ( iV )
        KV_E   ( iV )  =  KV_E   ( iV )  +  dT * F_1 ( iV ) * V_1 ( iV ) 
        SV_M_1 ( iV )  =  SV_M_1 ( iV )  +  F_1 ( iV ) * Weight_RK
        SV_E   ( iV )  =  SV_E   ( iV )  +  F_1 ( iV ) * V_1 ( iV ) * Weight_RK
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP  parallel do schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = 1, nV
        F_1 ( iV )  =  - G * M * N ( iV )  /  R ( iV ) ** 2
      end do
      !$OMP  end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        KV_M_1 ( iV )  =  KV_M_1 ( iV )  +  dT * F_1 ( iV )
        KV_E   ( iV )  =  KV_E   ( iV )  +  dT * F_1 ( iV ) * V_1 ( iV ) 
        SV_M_1 ( iV )  =  SV_M_1 ( iV )  +  F_1 ( iV ) * Weight_RK
        SV_E   ( iV )  =  SV_E   ( iV )  +  F_1 ( iV ) * V_1 ( iV ) * Weight_RK
      end do
      !$OMP end parallel do
      
    end if

  end procedure ApplySourcesKernel


end submodule FishboneMoncrief_Kernel
