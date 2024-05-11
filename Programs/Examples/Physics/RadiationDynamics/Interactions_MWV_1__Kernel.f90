#include "Preprocessor"

submodule ( Interactions_MWV_1__Form ) Interactions_MWV_1__Kernel

  use Basics 
  
  implicit none

contains


  module procedure ComputeAllKernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
         Xi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_Eq ( iV )
        Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
        Chi_H ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
         Xi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_Eq ( iV )
        Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
        Chi_H ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeAllKernel


  module procedure ComputeSingleKernel

     Xi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )  *  J_Eq ( iV )
    Chi_J ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )
    Chi_H ( iV )  =  Kappa  *  M ( iV )  *  N ( iV )

  end procedure ComputeSingleKernel


end submodule Interactions_MWV_1__Kernel
