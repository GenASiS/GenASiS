#include "Preprocessor"

submodule ( RiemannSolver_HLL__Form ) RiemannSolver_HLL__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      iF, &
      nV, &
      nF
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( F_I, dim = 1 )
    nF  =  size ( F_I, dim = 2 )
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV  =  1,  nV
        AP_I ( iV )  =  max ( 0.0_KDR, + EP_IL ( iV ), + EP_IR ( iV ) )
        AM_I ( iV )  =  max ( 0.0_KDR, - EM_IL ( iV ), - EM_IR ( iV ) )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP firstprivate ( SqrtTiny )
      do iF  =  1,  nF
        do iV  =  1,  nV

          F_I ( iV, iF ) &
            =  (    AP_I ( iV )  *  F_IL ( iV, iF ) &
                 +  AM_I ( iV )  *  F_IR ( iV, iF ) &
                 -  AP_I ( iV )  *  AM_I ( iV ) &
                    *  ( U_IR ( iV, iF )  -  U_IL ( iV, iF ) ) ) &
               /  max ( AP_I ( iV )  +  AM_I ( iV ),  SqrtTiny )

        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV
        AP_I ( iV )  =  max ( 0.0_KDR, + EP_IL ( iV ), + EP_IR ( iV ) )
        AM_I ( iV )  =  max ( 0.0_KDR, - EM_IL ( iV ), - EM_IR ( iV ) )
      end do
      !$OMP end parallel do
    
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP firstprivate ( SqrtTiny )
      do iF  =  1,  nF
        do iV  =  1,  nV

          F_I ( iV, iF ) &
            =  (    AP_I ( iV )  *  F_IL ( iV, iF ) &
                 +  AM_I ( iV )  *  F_IR ( iV, iF ) &
                 -  AP_I ( iV )  *  AM_I ( iV ) &
                    *  ( U_IR ( iV, iF )  -  U_IL ( iV, iF ) ) ) &
               /  max ( AP_I ( iV )  +  AM_I ( iV ),  SqrtTiny )

        end do
      end do
      !$OMP end parallel do
    
    end if

  end procedure ComputeKernel


end submodule RiemannSolver_HLL__Kernel
