#include "Preprocessor"

submodule ( Gravitation_N_SG__Form ) Gravitation_N_SG__Kernel

  use Basics
  implicit none
  
contains


  module procedure ComputeSourceKernel

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

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( FourPI_G )
      do iV = 1, nValues
        S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( FourPI_G )
      do iV = 1, nValues
        S ( iV )  =  FourPi_G  *  M ( iV )  *  N ( iV )
      end do
      !$OMP end parallel do
    
    end if

  end procedure ComputeSourceKernel


end submodule Gravitation_N_SG__Kernel
