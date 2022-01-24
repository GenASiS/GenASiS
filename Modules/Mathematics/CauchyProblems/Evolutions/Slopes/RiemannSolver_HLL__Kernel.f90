#include "Preprocessor"

submodule ( RiemannSolver_HLL__Form ) RiemannSolver_HLL__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      iF, &
      iF_F, &
      iF_B, &
      nV, &
      nF
    real ( KDR ) :: &
      SqrtTiny!, &
      ! dF_L, dF_R, &
      !  F_L,  F_R, &
      ! Sign, Scale
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( RSV, dim = 1 )
    nF  =  size ( iaFluxes )
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )
    
    associate &
      ( F_I  => RSV, &
        AP_I => RSV ( :, iAP ), &
        AM_I => RSV ( :, iAM ) )
    
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
      !$OMP private ( iF_B, iF_F ) firstprivate ( SqrtTiny )
      do iF  =  1,  nF
        do iV  =  1,  nV

          iF_B  =  iaBalanced ( iF )
          iF_F  =  iaFluxes ( iF )

          F_I ( iV, iF_F ) &
            =  (    AP_I ( iV )  *  F_IL ( iV, iF ) &
                 +  AM_I ( iV )  *  F_IR ( iV, iF ) &
                 -  AP_I ( iV )  *  AM_I ( iV ) &
                    *  ( U_IR ( iV, iF_B )  -  U_IL ( iV, iF_B ) ) ) &
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
      !$OMP private ( iF_B, iF_F ) firstprivate ( SqrtTiny )
      do iF  =  1,  nF
        do iV  =  1,  nV

          iF_B  =  iaBalanced ( iF )
          iF_F  =  iaFluxes ( iF )

          F_I ( iV, iF_F ) &
            =  (    AP_I ( iV )  *  F_IL ( iV, iF ) &
                 +  AM_I ( iV )  *  F_IR ( iV, iF ) &
                 -  AP_I ( iV )  *  AM_I ( iV ) &
                    *  ( U_IR ( iV, iF_B )  -  U_IL ( iV, iF_B ) ) ) &
               /  max ( AP_I ( iV )  +  AM_I ( iV ),  SqrtTiny )

        end do
      end do
      !$OMP end parallel do
    
! do iF  =  1,  nF
!   do iV  =  1,  nV
!     dF_L   =  F_I  ( iV, iF )  -  F_IL ( iV, iF )
!     dF_R   =  F_IR ( iV, iF )  -  F_I  ( iV, iF )
!      F_L   =  F_IL ( iV, iF )
!      F_R   =  F_IR ( iV, iF )
!     Sign   =  dF_L * dF_R
!     Scale  =  abs ( Sign ) / max ( abs ( F_L * F_R ), SqrtTiny )
!     if ( Sign  <  0.0_KDR .and. Scale > 1.0e-6) then
!       call Show ( [ iF, iV ], '>>> iF, iV' )
!       call Show ( [ F_IL ( iV, iF ), F_I ( iV, iF ), F_IR ( iV, iF ) ], &
!                   '>>> F_IL, F_I, F_IR' )
!       call Show ( Sign, '>>> Sign' )
!       call Show ( Scale, '>>> Scale' )
!     end if
!   end do !-- iV
! end do !-- iF

    end if
    
    end associate   !-- F_I, AP_I, AM_I

  end procedure ComputeKernel


end submodule RiemannSolver_HLL__Kernel
