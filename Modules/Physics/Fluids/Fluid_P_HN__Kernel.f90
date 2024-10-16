#include "Preprocessor"

submodule ( Fluid_P_HN__Form ) Fluid_P_HN__Kernel

  use Basics
  
  implicit none
  
  real ( KDR ) :: &
      OR_Shift, &
      MassDensity_CGS, &
      SpecificEnergy_CGS, &
      Pressure_CGS, &
      Speed_CGS, &
      MeV

#ifdef ENABLE_OMP_OFFLOAD
  !$OMP declare target to &
  !$OMP   ( OR_Shift, MassDensity_CGS, SpecificEnergy_CGS, Pressure_CGS, &
  !$OMP     Speed_CGS, MeV )
#endif
    
contains


  module procedure InitializeModuleVariablesKernel
  
    !-- Historical Oak Ridge Shift, accounting for nuclear binding energy
    OR_Shift = 8.9_KDR * UNIT % MEGA_ELECTRON_VOLT &
                 / CONSTANT % ATOMIC_MASS_UNIT
      
    MassDensity_CGS     =  UNIT % MASS_DENSITY_CGS
    SpecificEnergy_CGS  =  UNIT % ERG  /  UNIT % GRAM
    Pressure_CGS        =  UNIT % BARYE
    Speed_CGS           =  UNIT % CENTIMETER  /  UNIT % SECOND
    MeV                 =  UNIT % MEGA_ELECTRON_VOLT
    
  end procedure InitializeModuleVariablesKernel

  
  module procedure Apply_EOS_Prologue_A_Kernel

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
    
  end procedure Apply_EOS_Prologue_A_Kernel
  
  
  module procedure Apply_EOS_Prologue_S_Kernel
  

    !$OMP OMP_DECLARE_TARGET

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

  end procedure Apply_EOS_Prologue_S_Kernel
  
  
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


  module procedure Compute_N_V_E_YE_G_A_Kernel

    !-- Compute_DensityC_Velocity_EnergyC_ElectronFraction_All_Galileo

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

  end procedure Compute_N_V_E_YE_G_A_Kernel


  module procedure Compute_N_V_E_YE_G_S_Kernel
    
    !$OMP OMP_DECLARE_TARGET

    !-- Compute_DensityC_Velocity_EnergyC_ElectronFraction_Single_Galileo

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

  end procedure Compute_N_V_E_YE_G_S_Kernel


  module procedure Apply_EOS_Epilogue_A_Kernel

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
        Mu_N  ( iV )  =  Mu_N  ( iV ) * MeV
        Mu_P  ( iV )  =  Mu_P  ( iV ) * MeV
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E  ( iV ) * MeV
        
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
        Mu_N  ( iV )  =  Mu_N  ( iV ) * MeV
        Mu_P  ( iV )  =  Mu_P  ( iV ) * MeV
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E  ( iV ) * MeV

        !Error_A ( iV )  = Error_A ( iV ) + Error ( iV ) * 1.0_KDR
        
      end do
      !$OMP end parallel do
        
    end if
    
  end procedure Apply_EOS_Epilogue_A_Kernel


  module procedure Apply_EOS_Epilogue_S_Kernel
  
    !$OMP OMP_DECLARE_TARGET

!        if ( N ( iV ) == 0.0_KDR ) cycle 

        P ( iV )      =  P ( iV ) * Pressure_CGS
        T ( iV )      =  T ( iV ) * MeV
        N ( iV )      =  N ( iV ) / M ( iV ) * MassDensity_CGS
        E ( iV )      =  ( E ( iV ) * SpecificEnergy_CGS  +  OR_Shift ) &
                           * M ( iV ) * N ( iV )
        SS ( iV )     =  sqrt ( SS ( iV ) ) * Speed_CGS
        Mu_N  ( iV )  =  Mu_N  ( iV ) * MeV
        Mu_P  ( iV )  =  Mu_P  ( iV ) * MeV
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E  ( iV ) * MeV

        !Error_A ( iV )  = Error_A ( iV ) + Error ( iV ) * 1.0_KDR
        
  end procedure Apply_EOS_Epilogue_S_Kernel
  
  
  module procedure ComputeFromBalanced_S_Kernel
  
    !$OMP OMP_DECLARE_TARGET
    call Compute_N_V_E_YE_G_S_Kernel &
             ( D, S_1, S_2, S_3, G, DE, M, M_UU_11, M_UU_22, M_UU_33, &
               N_Min, E_Min, Y_Min, Y_Safe, iV, N, V_1, V_2, V_3, E, YE )
    call Apply_EOS_Prologue_S_Kernel &
           ( M, N, P, T, E, YE, M_Ref, N_Min, E_Min, T_Min, Y_Min, &
             Y_Safe, iV )
    call ComputeFromEnergy_S_Kernel &
           ( FV, EOS, T_L_N, T_L_T, T_Ye, E_Shift, ia_F_I, ia_F_O, &
             ia_E, iSolve, iV = iV )
    call Apply_EOS_Epilogue_S_Kernel &
           ( N, P, T, SS, E, Mu_N, Mu_P, Mu_NP, Mu_E, M, iV )

  
  end procedure ComputeFromBalanced_S_Kernel


end submodule Fluid_P_HN__Kernel
