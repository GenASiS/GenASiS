#include "Preprocessor"

submodule ( Chart_SL__Template ) Chart_SL__Kernel

  use Basics
  
  implicit none
  
contains

  
  module procedure CopyBoundaryKernel

    integer ( KDI ) :: &
      iV, jV, kV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                oBE ( 2 )  +  dBE ( 2 ) * jV, &
                oBE ( 3 )  +  dBE ( 3 ) * kV ) &
              = V ( oBI ( 1 )  +  dBI ( 1 ) * iV, &
                    oBI ( 2 )  +  dBI ( 2 ) * jV, &
                    oBI ( 3 )  +  dBI ( 3 ) * kV )
          end do 
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
      
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                oBE ( 2 )  +  dBE ( 2 ) * jV, &
                oBE ( 3 )  +  dBE ( 3 ) * kV ) &
              = V ( oBI ( 1 )  +  dBI ( 1 ) * iV, &
                    oBI ( 2 )  +  dBI ( 2 ) * jV, &
                    oBI ( 3 )  +  dBI ( 3 ) * kV )
          end do 
        end do
      end do
      !$OMP end parallel do
    
    end if
      
  end procedure CopyBoundaryKernel

  
  module procedure ReverseBoundaryKernel

    integer ( KDI ) :: &
      iV, jV, kV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                oBE ( 2 )  +  dBE ( 2 ) * jV, &
                oBE ( 3 )  +  dBE ( 3 ) * kV ) &
              = - V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                      oBE ( 2 )  +  dBE ( 2 ) * jV, &
                      oBE ( 3 )  +  dBE ( 3 ) * kV )
          end do 
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                oBE ( 2 )  +  dBE ( 2 ) * jV, &
                oBE ( 3 )  +  dBE ( 3 ) * kV ) &
              = - V ( oBE ( 1 )  +  dBE ( 1 ) * iV, &
                      oBE ( 2 )  +  dBE ( 2 ) * jV, &
                      oBE ( 3 )  +  dBE ( 3 ) * kV )
          end do 
        end do
      end do
      !$OMP end parallel do
    
    end if

  end procedure ReverseBoundaryKernel


end submodule Chart_SL__Kernel
