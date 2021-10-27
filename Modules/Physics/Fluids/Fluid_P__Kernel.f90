#include "Preprocessor"

submodule ( Fluid_P__Form ) Fluid_P__Kernel

  use Basics

  implicit none

contains


  module procedure Compute_D_S_G_G_Kernel
 	 
    !-- Compute_DensityB_Momentum_EnergyB_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( D )

    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP firstprivate ( N_Min, E_Min )
      do iV = 1, nV

        if ( N ( iV )  <=  N_Min  .or.  E ( iV )  <=  E_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
          E   ( iV )  =  E_Min
        end if

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP firstprivate ( N_Min, E_Min )
      do iV = 1, nV

        if ( N ( iV )  <=  N_Min  .or.  E ( iV )  <=  E_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
        end if

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_D_S_G_G_Kernel 	 	 


  module procedure Compute_N_V_E_G_Kernel

    !-- Compute_DensityC_Velocity_EnergyC_Galileo

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( N )
    
    if ( UseDevice ) then

      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP firstprivate ( N_Min, E_Min )
      do iV = 1, nV

        if ( D ( iV )  <=  N_Min  .or.  G ( iV )  <=  E_Min ) then
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          G   ( iV )  =  0.0_KDR
        end if

        N ( iV )    =  D ( iV )

        V_1 ( iV )  =  M_UU_11 ( iV )  &
                       *  S_1 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_2 ( iV )  =  M_UU_22 ( iV )  &
                       *  S_2 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_3 ( iV )  =  M_UU_33 ( iV )  &
                       *  S_3 ( iV )  /  ( M ( iV )  *  D ( iV ) )

        E ( iV )    =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                                +  S_2 ( iV ) * V_2 ( iV ) &
                                                +  S_3 ( iV ) * V_3 ( iV ) )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP firstprivate ( N_Min, E_Min )
      do iV = 1, nV

        if ( D ( iV )  <=  N_Min  .or.  G ( iV )  <=  E_Min ) then
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
        end if

        N ( iV )    =  D ( iV )

        V_1 ( iV )  =  M_UU_11 ( iV )  &
                       *  S_1 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_2 ( iV )  =  M_UU_22 ( iV )  &
                       *  S_2 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_3 ( iV )  =  M_UU_33 ( iV )  &
                       *  S_3 ( iV )  /  ( M ( iV )  *  D ( iV ) )

        E ( iV )    =  G ( iV )  -  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV ) &
                                                +  S_2 ( iV ) * V_2 ( iV ) &
                                                +  S_3 ( iV ) * V_3 ( iV ) )

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Compute_N_V_E_G_Kernel


  module procedure Compute_FS_G_Kernel

    !-- Compute_FluxSet_Galileo_Kernel

    integer :: &
      Delta_1, Delta_2, Delta_3
    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( F_D )
    
    Delta_1  =  int (    real ( 1 + iDim  -  abs ( 1 - iDim ) )  &
                      /  real ( 1 + iDim  +  abs ( 1 - iDim ) ) )

    Delta_2  =  int (    real ( 2 + iDim  -  abs ( 2 - iDim ) )  &
                      /  real ( 2 + iDim  +  abs ( 2 - iDim ) ) )

    Delta_3  =  int (    real ( 3 + iDim  -  abs ( 3 - iDim ) )  &
                      /  real ( 3 + iDim  +  abs ( 3 - iDim ) ) )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3 )
      do iV = 1, nV

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )  +  Delta_1  *  P ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )  +  Delta_2  *  P ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )  +  Delta_3  *  P ( iV )

        F_G ( iV )  =  ( G ( iV )  +  P ( iV ) )  *  V_Dim ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3 )
      do iV = 1, nV

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )  +  Delta_1  *  P ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )  +  Delta_2  *  P ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )  +  Delta_3  *  P ( iV )

        F_G ( iV )  =  ( G ( iV )  +  P ( iV ) )  *  V_Dim ( iV )

      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_FS_G_Kernel


  module procedure Compute_ES_G_Kernel

    !-- Compute_EigenspeedSet_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( EF_P )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        EF_P ( iV )  =  V_Dim ( iV )  +  sqrt ( M_UU_Dim ( iV ) )  *  SS ( iV ) 
        EF_M ( iV )  =  V_Dim ( iV )  -  sqrt ( M_UU_Dim ( iV ) )  *  SS ( iV ) 
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        EF_P ( iV )  =  V_Dim ( iV )  +  sqrt ( M_UU_Dim ( iV ) )  *  SS ( iV ) 
        EF_M ( iV )  =  V_Dim ( iV )  -  sqrt ( M_UU_Dim ( iV ) )  *  SS ( iV ) 
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_ES_G_Kernel


  module procedure Compute_S_UD_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( S_UD_22 )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )  +  P ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )  +  P ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )  +  P ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )  +  P ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_S_UD_Kernel


end submodule Fluid_P__Kernel
