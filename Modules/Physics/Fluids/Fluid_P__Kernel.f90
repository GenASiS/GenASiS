#include "Preprocessor"

submodule ( Fluid_P__Form ) Fluid_P__Kernel

  use Basics

  implicit none

contains


  module procedure Compute_D_S_G_DS_G_Kernel
 	 
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
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV

        if ( N ( iV )  <=  N_Min  .or.  E ( iV )  <=  E_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
          E   ( iV )  =  E_Min
        end if

        D ( iV )  =  N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

        DS ( iV )  =  SB ( iV )  *  N ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

        if ( N ( iV )  <=  N_Min  .or.  E ( iV )  <=  E_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
          E   ( iV )  =  E_Min
        end if

        D ( iV )  =  N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

        DS ( iV )  =  SB ( iV )  *  N ( iV )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_D_S_G_DS_G_Kernel 	 	 


  module procedure Compute_N_V_E_G_A_Kernel

    !-- Compute_DensityC_Velocity_EnergyC_Galileo_All

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
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV

        if ( D ( iV )  <=  N_Min  .or.  G ( iV )  <=  E_Min ) then
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          G   ( iV )  =  E_Min
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
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

        if ( D ( iV )  <=  N_Min  .or.  G ( iV )  <=  E_Min ) then
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
          G   ( iV )  =  E_Min
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

  end procedure Compute_N_V_E_G_A_Kernel


  module procedure Compute_N_V_E_G_S_Kernel

    !-- Compute_DensityC_Velocity_EnergyC_Galileo_Single

    if ( D ( iV )  <=  N_Min  .or.  G ( iV )  <=  E_Min ) then
      D   ( iV )  =  N_Min
      S_1 ( iV )  =  0.0_KDR
      S_2 ( iV )  =  0.0_KDR
      S_3 ( iV )  =  0.0_KDR
      G   ( iV )  =  E_Min
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

  end procedure Compute_N_V_E_G_S_Kernel


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


end submodule Fluid_P__Kernel
