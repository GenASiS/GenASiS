#include "Preprocessor"

submodule ( Interactions_C__Form ) Interactions_C__Kernel

  use Basics 
  
  implicit none

contains


  module procedure ComputeKernel

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
      
    nV  =  size ( J_EQ )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
         Xi_J ( iV )  =  Kappa_A  *  J_EQ ( iV )
        Chi_J ( iV )  =  Kappa_A
        Chi_H ( iV )  =  Kappa_A
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
         Xi_J ( iV )  =  Kappa_A  *  J_EQ ( iV )
        Chi_J ( iV )  =  Kappa_A
        Chi_H ( iV )  =  Kappa_A
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule Interactions_C__Kernel
