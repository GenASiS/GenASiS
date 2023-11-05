#include "Preprocessor"

submodule ( Interactions_MWV_3__Form ) Interactions_MWV_3__Kernel

  use Basics 
  
  implicit none

contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      S, &
      S_EQ
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP private ( S, S_Eq )
      do iV = 1, nV
         Xi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_Eq ( iV )
        Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
        Chi_H ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP private ( S, S_Eq )
      do iV = 1, nV

        S     =  1.0_KDR  -  Ratio_P  *  k_B  *  TP ( iV )  /  E_Max
        S_Eq  =  1.0_KDR  -  Ratio_P  *  k_B  *   T ( iV )  /  E_Max

         Xi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  * S_Eq * J_Eq ( iV ) &
                         *  ( T ( iV ) / T_0 ) ** 1.5_KDR  
        Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  S  & 
                         *  ( T ( iV ) / T_0 ) ** 1.5_KDR  
        Chi_H ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  S  &
                         *  ( T ( iV ) / T_0 ) ** 1.5_KDR  

      end do
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule Interactions_MWV_3__Kernel
