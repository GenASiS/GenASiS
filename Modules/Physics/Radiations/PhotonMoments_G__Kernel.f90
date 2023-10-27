#include "Preprocessor"

submodule ( PhotonMoments_G__Form ) PhotonMoments_G__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_SP_Kernel

    !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        TP ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
        TE ( iV )  =  T ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        TP ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
        TE ( iV )  =  T ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_SP_Kernel


end submodule PhotonMoments_G__Kernel
