#include "Preprocessor"

submodule ( IncrementDivergence_FV__Form ) IncrementDivergence_FV__Kernel

  use Basics
  
  implicit none
  
contains

  module procedure ComputeReconstruction_CSL_Kernel

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

!    V_IL  =  cshift ( V  +  0.5_KDR * dX * dVdX, shift = -1, dim = iD )

!    V_IR  =  V  -  0.5_KDR * dX * dVdX

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
    uV ( iD ) = size ( V, dim = iD ) - oV + 1 
      
    iaS = 0
    iaS ( iD ) = -1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            V_IL ( iV, jV, kV )  &
              =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                 +  dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                    *  dVdX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
            
            V_IR ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  &
                 -  dX_L ( iV, jV, kV )  *  dVdX ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
              
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            V_IL ( iV, jV, kV )  &
              =  V ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )  &
                 +  dX_R ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                    *  dVdX ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) )
                    
            V_IR ( iV, jV, kV )  &
              =  V ( iV, jV, kV )  &
                 -  dX_L ( iV, jV, kV )  *  dVdX ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if
        
  end procedure ComputeReconstruction_CSL_Kernel


  module procedure ComputeLogDerivative_CSL_Kernel

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

    lV = 1
    where ( shape ( V ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( V ) > 1 )
      uV = shape ( V ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dLVdX ( iV, jV, kV ) &
              =  (    A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                   -  A_I ( iV, jV, kV )  ) &
                 /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else
        
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dLVdX ( iV, jV, kV ) &
              =  (    A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                   -  A_I ( iV, jV, kV )  ) &
                 /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
      
    end if
      
  end procedure ComputeLogDerivative_CSL_Kernel

  
  module procedure ComputeIncrement_CSL_Kernel
    
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

!    dU  =  dU  +  dT * ( VJ_I * F_I  &
!                         -  cshift ( VJ_I * F_I, shift = +1, dim = iD ) ) &
!                       / ( VJ * dX )

    lV = 1
    where ( shape ( dU ) > 1 )
      lV = oV + 1
    end where
    
    uV = 1
    where ( shape ( dU ) > 1 )
      uV = shape ( dU ) - oV
    end where
      
    iaS = 0
    iaS ( iD ) = +1
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dU ( iV, jV, kV ) &
              = dU ( iV, jV, kV )  &
                +  dT * (    A_I ( iV, jV, kV ) &
                             *  F_I ( iV, jV, kV ) &
                          -  A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                             *  F_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                        /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else
        
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV, iaVS )
      do kV = lV ( 3 ), uV ( 3 ) 
        do jV = lV ( 2 ), uV ( 2 )
          do iV = lV ( 1 ), uV ( 1 )

            iaVS = [ iV, jV, kV ] + iaS

            dU ( iV, jV, kV ) &
              = dU ( iV, jV, kV )  &
                +  dT * (    A_I ( iV, jV, kV ) &
                             *  F_I ( iV, jV, kV ) &
                          -  A_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) &
                             *  F_I ( iaVS ( 1 ), iaVS ( 2 ), iaVS ( 3 ) ) ) &
                        /  V ( iV, jV, kV )

          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if
    
  end procedure ComputeIncrement_CSL_Kernel
  
  
  module procedure RecordBoundaryFluence_CSL_Kernel

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
            BF ( iV, jV, kV ) &
              =  BF ( iV, jV, kV ) &
                 +  Factor &
                    *   F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

    else 
    
      !$OMP  parallel do collapse ( 3 ) &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV, jV, kV ) 
      do kV = 1, nB ( 3 )
        do jV = 1, nB ( 2 )
          do iV = 1, nB ( 1 )
            BF ( iV, jV, kV ) &
              =  BF ( iV, jV, kV ) &
                 +  Factor &
                    *   F ( oB ( 1 ) + iV, oB ( 2 ) + jV, oB ( 3 ) + kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP  end parallel do
    
    end if

  end procedure RecordBoundaryFluence_CSL_Kernel


end submodule IncrementDivergence_FV__Kernel
