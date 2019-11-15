#include "Preprocessor"

submodule ( RayleighTaylor_Form ) RayleighTaylor_Kernel

  use Basics

  implicit none

contains

  
  module procedure ApplySourcesKernel
  
    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( N )
    
    !call Show ( UseDevice, '==== UseDevice ====')
    !call Show ( OnDevice ( KVM ), '==== onDevice KVM ====')
    !call Show ( OnDevice ( N ), '==== onDevice N ====')
    !call Show ( OnDevice ( VY ), '==== onDevice VY ====')
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, size ( N )
        KVM ( iV ) = KVM ( iV )  -  dT * N ( iV ) * A
        KVE ( iV ) = KVE ( iV )  -  dT * N ( iV ) * A * VY ( iV )
        SVM ( iV ) = SVM ( iV )  -  B * N ( iV ) * A 
        SVE ( iV ) = SVE ( iV )  -  B * N ( iV ) * A * VY ( iV )
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, size ( N )
        KVM ( iV ) = KVM ( iV )  -  dT * N ( iV ) * A
        KVE ( iV ) = KVE ( iV )  -  dT * N ( iV ) * A * VY ( iV )
        SVM ( iV ) = SVM ( iV )  -  B * N ( iV ) * A 
        SVE ( iV ) = SVE ( iV )  -  B * N ( iV ) * A * VY ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ApplySourcesKernel  


end submodule RayleighTaylor_Kernel
