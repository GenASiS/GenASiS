#include "Preprocessor"

submodule ( Difference_Form ) Difference_Kernel

  use Basics 
  
  implicit none
  
contains

  
  module procedure ComputeChart_SL_Kernel

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      iaS, &
      iaVS, &   
      lV, uV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
!    dV_I = V - cshift ( V, shift = -1, dim = iD )

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    lV ( iD ) = oV

    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD )
    
    iaS = 0
    iaS ( iD ) = -1
    
    if ( UseDevice ) then
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_I ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else 
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dV_I ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  -  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if
       
  end procedure ComputeChart_SL_Kernel


end submodule Difference_Kernel
