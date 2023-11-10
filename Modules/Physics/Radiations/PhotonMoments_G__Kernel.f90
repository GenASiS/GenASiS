#include "Preprocessor"

submodule ( PhotonMoments_G__Form ) PhotonMoments_G__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_Eq_Kernel

    !-- Compute_J_Eq_Photons_Grey_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      a
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( a )
      do iV = 1, nV
        J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( a )
      do iV = 1, nV
        J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_Eq_Kernel


  module procedure Compute_SP_Kernel

    !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      a
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J )

    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( a )
      do iV = 1, nV
        T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( a )
      do iV = 1, nV
        T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_SP_Kernel


end submodule PhotonMoments_G__Kernel
