#include "Preprocessor"

submodule ( RadiationMoments_BM__Form ) RadiationMoments_BM__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_E_S_G_Kernel

    !-- Compute_BalancedEnergy_Momentum_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      H
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( E )

    if ( UseDevice ) then

    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( H )
      do iV = 1, nV

        if ( J ( iV )  <  0.0_KDR ) &
          J ( iV )  =  0.0_KDR

        H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                     +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                     +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        if ( H  >  J ( iV ) ) then

          H_1 ( iV )  =  ( H_1 ( iV )  /  H )  *  J ( iV )
          H_2 ( iV )  =  ( H_2 ( iV )  /  H )  *  J ( iV )
          H_3 ( iV )  =  ( H_3 ( iV )  /  H )  *  J ( iV )
  
          H  =  sqrt (    M_DD_11 ( iV )  *  H_1 ( iV ) ** 2  &
                       +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                       +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 )

        end if

        !-- FIXME: Add velocity dependence

        E ( iV )  =  J ( iV )

        S_1 ( iV )  =  M_DD_11 ( iV )  *  H_1 ( iV )
        S_2 ( iV )  =  M_DD_22 ( iV )  *  H_2 ( iV )
        S_3 ( iV )  =  M_DD_33 ( iV )  *  H_3 ( iV )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_E_S_G_Kernel


end submodule RadiationMoments_BM__Kernel
