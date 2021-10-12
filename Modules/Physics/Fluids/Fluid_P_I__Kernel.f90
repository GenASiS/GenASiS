#include "Preprocessor"

submodule ( Fluid_P_I__Form ) Fluid_P_I__Kernel

  use Basics

  implicit none

contains


  module procedure Apply_EOS_I_T_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      SqrtHuge
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    SqrtHuge = sqrt ( huge ( 1.0_KDR ) )

    nValues = size ( P )
    
    if ( UseDevice ) then
    
      !-- FIXME: Performance tuning, test if breaking this into multiple loops
      !          work better
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
      !$OMP firstprivate ( SqrtHuge, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        E ( iV )  =  C_V  *  N ( iV )  *  T ( iV )

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then

          SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                      *  ( N0 / N ( iV ) ) ** Gamma ) 

          CS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

        else
          SB ( iV )  =  - SqrtHuge
          CS ( iV )  =    0.0_KDR
        end if

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP firstprivate ( SqrtHuge, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        E ( iV )  =  C_V  *  N ( iV )  *  T ( iV )

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        if ( N ( iV ) > 0.0_KDR .and. P ( iV ) > 0.0_KDR ) then

          SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                      *  ( N0 / N ( iV ) ) ** Gamma ) 

          CS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

        else
          SB ( iV )  =  - SqrtHuge
          CS ( iV )  =    0.0_KDR
        end if

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Apply_EOS_I_T_Kernel


end submodule Fluid_P_I__Kernel
