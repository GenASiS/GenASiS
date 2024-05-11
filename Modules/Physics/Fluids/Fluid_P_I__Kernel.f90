#include "Preprocessor"

submodule ( Fluid_P_I__Form ) Fluid_P_I__Kernel

  use Basics

  implicit none

contains


  module procedure Apply_EOS_I_T_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( P )
    
    if ( UseDevice ) then
    
      !-- FIXME: Performance tuning, test if breaking this into multiple loops
      !          work better
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
      !$OMP firstprivate ( M_Ref, N_Min, T_Min, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        M ( iV )  =  M_Ref

        if ( N ( iV )  <  N_Min ) &
          N ( iV )  =  N_Min
        if ( T ( iV )  <  T_Min ) &
          T ( iV )  =  T_Min

        E ( iV )  =  C_V  *  N ( iV )  *  T ( iV )

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        SS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP firstprivate ( M_Ref, N_Min, T_Min, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        M ( iV )  =  M_Ref

        if ( N ( iV )  <  N_Min ) &
          N ( iV )  =  N_Min
        if ( T ( iV )  <  T_Min ) &
          T ( iV )  =  T_Min

        E ( iV )  =  C_V  *  N ( iV )  *  T ( iV )

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        SS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Apply_EOS_I_T_Kernel


  module procedure Apply_EOS_I_E_A_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( P )

    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
      !$OMP firstprivate ( M_Ref, N_Min, E_Min, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        M ( iV )  =  M_Ref

        if ( N ( iV )  <  N_Min ) &
          N ( iV )  =  N_Min
        if ( E ( iV )  <  E_Min ) &
          E ( iV )  =  E_Min

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        T ( iV )  =  E ( iV )  /  ( C_V  *  N ( iV ) )

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        SS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP firstprivate ( M_Ref, N_Min, E_Min, Gamma, C_V, N0, P0 )
      do iV = 1, nValues

        M ( iV )  =  M_Ref

        if ( N ( iV )  <  N_Min ) &
          N ( iV )  =  N_Min
        if ( E ( iV )  <  E_Min ) &
          E ( iV )  =  E_Min

        P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

        T ( iV )  =  E ( iV )  /  ( C_V  *  N ( iV ) )

        SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                    *  ( N0 / N ( iV ) ) ** Gamma ) 

        SS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Apply_EOS_I_E_A_Kernel


  module procedure Apply_EOS_I_E_S_Kernel

    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    if ( UseDevice ) then
      
    else

      M ( iV )  =  M_Ref

      if ( N ( iV )  <  N_Min ) &
        N ( iV )  =  N_Min
      if ( E ( iV )  <  E_Min ) &
        E ( iV )  =  E_Min

      P ( iV )  =  ( Gamma - 1.0_KDR )  *  E ( iV ) 

      T ( iV )  =  E ( iV )  /  ( C_V  *  N ( iV ) )

      SB ( iV )  =  C_V  *  log ( P ( iV ) / P0  &
                                  *  ( N0 / N ( iV ) ) ** Gamma ) 

      SS ( iV )  =  sqrt ( Gamma * P ( iV ) / ( M ( iV ) * N ( iV ) ) )

    end if

  end procedure Apply_EOS_I_E_S_Kernel


end submodule Fluid_P_I__Kernel
