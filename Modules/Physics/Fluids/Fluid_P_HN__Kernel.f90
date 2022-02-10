#include "Preprocessor"

submodule ( Fluid_P_HN__Form ) Fluid_P_HN__Kernel

  use Basics
  
  implicit none
  
contains

  
  module procedure Apply_EOS_PrologueKernel

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
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues

        M ( iV )  =  M_Ref
        
        if ( N ( iV )  <=  N_Min ) then
          N ( iV )   =  N_Min
          YE ( iV )  =  Y_Safe
        end if
        if ( E ( iV )  <=  E_Min ) & 
          E ( iV )  =  E_Min
        if ( T ( iV )  <=  T_Min ) & 
          T ( iV )  =  T_Min
        if ( YE ( iV )  <=  Y_Min ) &
          YE ( iV )  =  Y_Min

        E ( iV )  =  ( E ( iV )  /  ( M ( iV )  *  N ( iV ) )  -  OR_Shift ) &
                     /  SpecificEnergy_CGS
        N ( iV )  =  M ( iV )  *  N ( iV )  /  MassDensity_CGS
        T ( iV )  =  T ( iV )  /  MeV
        P ( iV )  =  P ( iV )  /  Pressure_CGS

      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues

        M ( iV )   =  M_Ref

        if ( N ( iV )  <=  N_Min ) then
          N ( iV )   =  N_Min
          YE ( iV )  =  Y_Safe
        end if
        if ( E ( iV )  <=  E_Min ) & 
          E ( iV )  =  E_Min
        if ( T ( iV )  <=  T_Min ) & 
          T ( iV )  =  T_Min
        if ( YE ( iV )  <=  Y_Min ) &
          YE ( iV )  =  Y_Min

        E ( iV )  =  ( E ( iV )  /  ( M ( iV )  *  N ( iV ) )  -  OR_Shift ) &
                     /  SpecificEnergy_CGS
        N ( iV )  =  M ( iV )  *  N ( iV )  /  MassDensity_CGS
        T ( iV )  =  T ( iV )  /  MeV
        P ( iV )  =  P ( iV )  /  Pressure_CGS

      end do
      !$OMP end parallel do
    
    end if
    
  end procedure Apply_EOS_PrologueKernel
  
  
  module procedure Compute_D_S_G_DE_G_Kernel
 	 
    !-- Compute_DensityB_Momentum_EnergyB_DensityElectronB_Galileo_Kernel

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
          YE  ( iV )  =  Y_Safe
        end if

        D ( iV )  =  N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

        DE ( iV )  =  YE ( iV )  *  N ( iV ) 	 	 
       
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
          YE  ( iV )  =  Y_Safe
        end if

        D ( iV )  =  N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

        DE ( iV )  =  YE ( iV )  *  N ( iV ) 	 	 
       
      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_D_S_G_DE_G_Kernel 	 	 


  module procedure Compute_N_V_E_YE_G_Kernel

    !-- Compute_DensityC_Velocity_EnergyC_ElectronFraction_Galileo

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
          DE  ( iV )  =  Y_Safe * N_Min
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

        YE ( iV )   =  DE ( iV )  /  N ( iV )

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
          DE  ( iV )  =  Y_Safe * N_Min
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

        YE ( iV )   =  DE ( iV )  /  N ( iV )

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Compute_N_V_E_YE_G_Kernel


  module procedure Apply_EOS_EpilogueKernel

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
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues

!        if ( N ( iV ) == 0.0_KDR ) cycle 

        P ( iV )      =  P ( iV ) * Pressure_CGS
        T ( iV )      =  T ( iV ) * MeV
        N ( iV )      =  N ( iV ) / M ( iV ) * MassDensity_CGS
        E ( iV )      =  ( E ( iV ) * SpecificEnergy_CGS  +  OR_Shift ) &
                           * M ( iV ) * N ( iV )
        SS ( iV )     =  sqrt ( SS ( iV ) ) * Speed_CGS
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
        
        !Error_A ( iV )  = Error_A ( iV ) + Error ( iV ) * 1.0_KDR
        
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues

!        if ( N ( iV ) == 0.0_KDR ) cycle 

        P ( iV )      =  P ( iV ) * Pressure_CGS
        T ( iV )      =  T ( iV ) * MeV
        N ( iV )      =  N ( iV ) / M ( iV ) * MassDensity_CGS
        E ( iV )      =  ( E ( iV ) * SpecificEnergy_CGS  +  OR_Shift ) &
                           * M ( iV ) * N ( iV )
        SS ( iV )     =  sqrt ( SS ( iV ) ) * Speed_CGS
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E ( iV ) * MeV

        !Error_A ( iV )  = Error_A ( iV ) + Error ( iV ) * 1.0_KDR
        
      end do
      !$OMP end parallel do
        
    end if
    
  end procedure Apply_EOS_EpilogueKernel


end submodule Fluid_P_HN__Kernel
