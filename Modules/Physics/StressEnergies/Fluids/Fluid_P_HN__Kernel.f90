#include "Preprocessor"

submodule ( Fluid_P_HN__Form ) Fluid_P_HN__Kernel

  use Basics
  
  implicit none
  
contains

  
  module procedure Compute_DE_G_Kernel
 	 
    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( DE )

    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nValues
        DE ( iV )  =  YE ( iV )  *  N ( iV )
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = 1, nValues
        DE ( iV )  =  YE ( iV )  *  N ( iV )
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_DE_G_Kernel


  module procedure Compute_YE_G_Kernel
 	 
    !-- Compute_ProtonFraction_Galilean_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
           
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
                         
    nValues = size ( YE )
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nValues
        
        if ( DE ( iV ) < 0.0_KDR ) &
          DE ( iV )  =  0.0_KDR
        
        if ( N ( iV ) > 0.0_KDR ) then
          YE ( iV ) = DE ( iV ) / N ( iV )
        else
          YE ( iV ) = 0.0_KDR
        end if
        
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = 1, nValues
        
        if ( DE ( iV ) < 0.0_KDR ) &
          DE ( iV )  =  0.0_KDR
        
        if ( N ( iV ) > 0.0_KDR ) then
          YE ( iV ) = DE ( iV ) / N ( iV )
        else
          YE ( iV ) = 0.0_KDR
        end if
        
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_YE_G_Kernel
  
  
  module procedure Apply_EOS_HN_T_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp, &
      keyerr, &
      Rank
    real ( KDR ) :: &
      rfeps, &
      Rho_Temp, &
      T_Temp, &
!      cs2, dedt, dpderho, dpdrhoe, munu
      cs2, dedt, dpderho, dpdrhoe, mu_n, mu_p
    logical ( KDL ) :: &
      UseDevice
           
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( P )

    !-- Compute P, E, Gamma, SB from N, T, YE

    rfeps = 1.0e-9_KDR
    keytemp = 1_KDI

    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP& private ( iV, keyerr, Rho_Temp, T_Temp, &
      !$OMP&           cs2, dedt, dpderho, dpdrhoe, mu_n, mu_p ) &
      !$OMP& firstprivate ( rfeps, keytemp, nValues )
      do iV = 1, nValues
        
        if ( N ( iV ) == 0.0_KDR ) cycle 
        
        Rho_Temp = M ( iV ) * N ( iV ) / MassDensity_CGS
        T_Temp   = T ( iV ) / MeV
        E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                   / SpecificEnergy_CGS
        P ( iV ) = P ( iV ) / Pressure_CGS
        ! call nuc_eos_short &
        !        ( N_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
        !          cs2, dedt, dpderho, dpdrhoe, munu, &
        !          keytemp, keyerr, rfeps )
        call nuc_eos_full &
               ( Rho_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
                 cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), X_A ( iV ), &
                 X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
                 mu_n, mu_p, Mu_NP ( iV ), keytemp, keyerr, rfeps )

        !if ( keyerr /= 0 ) then
        !  Rank = PROGRAM_HEADER % Communicator % Rank
        !  call Show ( 'EOS error', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !  call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !  call Show ( 'Apply_EOS_HN_T_Kernel', 'subroutine', &
        !              CONSOLE % WARNING, DisplayRankOption = Rank )
        !  call Show ( Rank, 'Rank', DisplayRankOption = Rank )
        !  call Show ( iV, 'iV', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !end if
  !      call nuc_eos_one ( Rho_Temp, T_Temp, YE ( iV ), Gamma ( iV ), 19 )
        P ( iV )      =  P ( iV ) * Pressure_CGS
        E ( iV )      =  E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
        E ( iV )      =  E ( iV ) * M ( iV ) * N ( iV )
        CS ( iV )     =  sqrt ( cs2 ) * Speed_CGS
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
        
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP& firstprivate ( rfeps, keytemp, nValues )
      do iV = 1, nValues
        
        if ( N ( iV ) == 0.0_KDR ) cycle 
        
        Rho_Temp = M ( iV ) * N ( iV ) / MassDensity_CGS
        T_Temp   = T ( iV ) / MeV
        E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                   / SpecificEnergy_CGS
        P ( iV ) = P ( iV ) / Pressure_CGS
        ! call nuc_eos_short &
        !        ( N_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
        !          cs2, dedt, dpderho, dpdrhoe, munu, &
        !          keytemp, keyerr, rfeps )
        call nuc_eos_full &
               ( Rho_Temp, T_Temp, YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
                 cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), X_A ( iV ), &
                 X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), Mu_E ( iV ), &
                 mu_n, mu_p, Mu_NP ( iV ), keytemp, keyerr, rfeps )

        !if ( keyerr /= 0 ) then
        !  Rank = PROGRAM_HEADER % Communicator % Rank
        !  call Show ( 'EOS error', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !  call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !  call Show ( 'Apply_EOS_HN_T_Kernel', 'subroutine', &
        !              CONSOLE % WARNING, DisplayRankOption = Rank )
        !  call Show ( Rank, 'Rank', DisplayRankOption = Rank )
        !  call Show ( iV, 'iV', CONSOLE % WARNING, &
        !              DisplayRankOption = Rank )
        !end if
  !      call nuc_eos_one ( Rho_Temp, T_Temp, YE ( iV ), Gamma ( iV ), 19 )
        P ( iV )      =  P ( iV ) * Pressure_CGS
        E ( iV )      =  E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
        E ( iV )      =  E ( iV ) * M ( iV ) * N ( iV )
        CS ( iV )     =  sqrt ( cs2 ) * Speed_CGS
        Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
        Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
        
      end do
      !$OMP  end parallel do
    
    end if
      
  end procedure Apply_EOS_HN_T_Kernel
  

  module procedure Apply_EOS_HN_SB_E_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp_e, &
      keytemp_s, &
      keyerr, &
      Rank
    real ( KDR ) :: &
      rfeps, &
      Rho_Temp, &
!      cs2, dedt, dpderho, dpdrhoe, munu
      cs2, dedt, dpderho, dpdrhoe, mu_n, mu_p

    nValues = size ( P )

    !-- Compute P, T, Gamma, SB from N, E, YE

    rfeps = 1.0e-9_KDR
    keytemp_e = 0_KDI
    keytemp_s = 2_KDI

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      Rho_Temp   = M ( iV ) * N ( iV ) / MassDensity_CGS
      E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                 / SpecificEnergy_CGS
      P ( iV ) = P ( iV ) / Pressure_CGS
      T ( iV ) = T ( iV ) / MeV
      if ( Shock ( iV ) > 0.0_KDR ) then
        ! call nuc_eos_short &
        !        ( N_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
        !          cs2, dedt, dpderho, dpdrhoe, munu, &
        !          keytemp_e, keyerr, rfeps )
        call nuc_eos_full &
               ( Rho_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), &
                 SB ( iV ), cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), &
                 X_A ( iV ), X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), &
                 Mu_E ( iV ), mu_n, mu_p, Mu_NP ( iV ), &
                 keytemp_e, keyerr, rfeps )
      else !-- not Shock
        ! call nuc_eos_short &
        !        ( N_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
        !         cs2, dedt, dpderho, dpdrhoe, munu, &
        !         keytemp_s, keyerr, rfeps )
        call nuc_eos_full &
               ( Rho_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), &
                 SB ( iV ), cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), &
                 X_A ( iV ), X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), &
                 Mu_E ( iV ), mu_n, mu_p, Mu_NP ( iV ), &
                 keytemp_s, keyerr, rfeps )
      end if !-- Shock
      if ( keyerr /= 0 ) then
        Rank = PROGRAM_HEADER % Communicator % Rank
        call Show ( 'EOS error', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Apply_EOS_HN_SB_E_Kernel', 'subroutine', &
                    CONSOLE % WARNING, DisplayRankOption = Rank )
        call Show ( Rank, 'Rank', DisplayRankOption = Rank )
        call Show ( iV, 'iV', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
      end if
!       call nuc_eos_one ( N_Temp, T ( iV ), YE ( iV ), Gamma ( iV ), 19 )
      E ( iV )      =  E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
      E ( iV )      =  E ( iV ) * M ( iV ) * N ( iV )
      P ( iV )      =  P ( iV ) * Pressure_CGS
      CS ( iV )     =  sqrt ( cs2 ) * Speed_CGS
      T ( iV )      =  T ( iV ) * MeV
      Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
      Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
    end do
    !$OMP end parallel do
    
  end procedure Apply_EOS_HN_SB_E_Kernel


  module procedure Apply_EOS_HN_E_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp_e, &
      keyerr, &
      Rank
    real ( KDR ) :: &
      rfeps, &
      Rho_Temp, &
!      cs2, dedt, dpderho, dpdrhoe, munu
      cs2, dedt, dpderho, dpdrhoe, mu_n, mu_p

    nValues = size ( P )

    !-- Compute P, T, Gamma, SB from N, E, YE

    rfeps = 1.0e-9_KDR
    keytemp_e = 0_KDI

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      Rho_Temp   = M ( iV ) * N ( iV ) / MassDensity_CGS
      E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                 / SpecificEnergy_CGS
      P ( iV ) = P ( iV ) / Pressure_CGS
      T ( iV ) = T ( iV ) / MeV
      ! call nuc_eos_short &
      !        ( N_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
      !          cs2, dedt, dpderho, dpdrhoe, munu, &
      !          keytemp_e, keyerr, rfeps )
      call nuc_eos_full &
             ( Rho_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), &
               SB ( iV ), cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), &
               X_A ( iV ), X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), &
               Mu_E ( iV ), mu_n, mu_p, Mu_NP ( iV ), &
               keytemp_e, keyerr, rfeps )
      if ( keyerr /= 0 ) then
        Rank = PROGRAM_HEADER % Communicator % Rank
        call Show ( 'EOS error', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Apply_EOS_HN_SB_E_Kernel', 'subroutine', &
                    CONSOLE % WARNING, DisplayRankOption = Rank )
        call Show ( Rank, 'Rank', DisplayRankOption = Rank )
        call Show ( iV, 'iV', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
      end if
!       call nuc_eos_one ( N_Temp, T ( iV ), YE ( iV ), Gamma ( iV ), 19 )
      E ( iV )      =  E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
      E ( iV )      =  E ( iV ) * M ( iV ) * N ( iV )
      P ( iV )      =  P ( iV ) * Pressure_CGS
      CS ( iV )     =  sqrt ( cs2 ) * Speed_CGS
!      CS ( iV )     =  sqrt ( 2.0_KDR * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
      T ( iV )      =  T ( iV ) * MeV
      Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
      Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
    end do
    !$OMP end parallel do
    
  end procedure Apply_EOS_HN_E_Kernel


  module procedure Apply_EOS_HN_SB_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues, &
      keytemp_s, &
      keyerr, &
      Rank
    real ( KDR ) :: &
      rfeps, &
      Rho_Temp, &
!      cs2, dedt, dpderho, dpdrhoe, munu
      cs2, dedt, dpderho, dpdrhoe, mu_n, mu_p

    nValues = size ( P )

    !-- Compute P, T, Gamma, SB from N, E, YE

    rfeps = 1.0e-9_KDR
    keytemp_s = 2_KDI

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      if ( N ( iV ) == 0.0_KDR ) cycle 
      Rho_Temp   = M ( iV ) * N ( iV ) / MassDensity_CGS
      E ( iV ) = ( E ( iV ) / ( M ( iV ) * N ( iV ) )  -  OR_Shift ) &
                 / SpecificEnergy_CGS
      P ( iV ) = P ( iV ) / Pressure_CGS
      T ( iV ) = T ( iV ) / MeV
      ! call nuc_eos_short &
      !        ( N_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), SB ( iV ), &
      !         cs2, dedt, dpderho, dpdrhoe, munu, &
      !         keytemp_s, keyerr, rfeps )
      call nuc_eos_full &
             ( Rho_Temp, T ( iV ), YE ( iV ), E ( iV ), P ( iV ), &
               SB ( iV ), cs2, dedt, dpderho, dpdrhoe, X_He ( iV ), &
               X_A ( iV ), X_N ( iV ), X_P ( iV ), A ( iV ), Z ( iV ), &
               Mu_E ( iV ), mu_n, mu_p, Mu_NP ( iV ), &
               keytemp_s, keyerr, rfeps )
      if ( keyerr /= 0 ) then
        Rank = PROGRAM_HEADER % Communicator % Rank
        call Show ( 'EOS error', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Fluid_P_HN__Form', 'module', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
        call Show ( 'Apply_EOS_HN_SB_E_Kernel', 'subroutine', &
                    CONSOLE % WARNING, DisplayRankOption = Rank )
        call Show ( Rank, 'Rank', DisplayRankOption = Rank )
        call Show ( iV, 'iV', CONSOLE % WARNING, &
                    DisplayRankOption = Rank )
      end if
!       call nuc_eos_one ( N_Temp, T ( iV ), YE ( iV ), Gamma ( iV ), 19 )
      E ( iV )      =  E ( iV ) * SpecificEnergy_CGS  +  OR_Shift
      E ( iV )      =  E ( iV ) * M ( iV ) * N ( iV )
      P ( iV )      =  P ( iV ) * Pressure_CGS
      CS ( iV )     =  sqrt ( cs2 ) * Speed_CGS
!      CS ( iV )     =  sqrt ( 2.0_KDR * P ( iV ) / ( M ( iV ) * N ( iV ) ) )
      T ( iV )      =  T ( iV ) * MeV
      Mu_NP ( iV )  =  Mu_NP ( iV ) * MeV
      Mu_E  ( iV )  =  Mu_E ( iV ) * MeV
    end do
    !$OMP end parallel do
    
  end procedure Apply_EOS_HN_SB_Kernel


  module procedure ComputeRawFluxesKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( F_DE )
    
    if ( UseDevice ) then
       
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV )
      do iV = 1, nValues
        F_DE ( iV )  =  DE ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV )
      do iV = 1, nValues
        F_DE ( iV )  =  DE ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure ComputeRawFluxesKernel


  module procedure ComputeCenterStatesKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      AM_VL, &
      AM_AC, &
      AM_AC_Inv, &
      AP_VR, &
      AP_AC, &
      AP_AC_Inv, &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( AC_I )
    SqrtTiny = sqrt ( tiny ( 0.0_KDR ) )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
      !$OMP& private (  AM_VL, AM_AC, AM_AC_Inv, AP_VR, AP_AC, AP_AC_Inv ) &
      !$OMP& firstprivate ( SqrtTiny )
      do iV = 1, nValues

        AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        DE_ICL ( iV )  =  DE_IL ( iV ) * AM_VL * AM_AC_Inv
        DE_ICR ( iV )  =  DE_IR ( iV ) * AP_VR * AP_AC_Inv

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else 

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
      !$OMP& private (  AM_VL, AM_AC, AM_AC_Inv, AP_VR, AP_AC, AP_AC_Inv ) &
      !$OMP& firstprivate ( SqrtTiny )
      do iV = 1, nValues

        AM_VL     =  AM_I ( iV )  +  V_Dim_IL ( iV )
        AM_AC     =  AM_I ( iV )  +  AC_I ( iV )
  !      AM_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AM_AC ), SqrtTiny ), AM_AC )
        AM_AC_Inv =  1.0_KDR &
                     / max ( abs ( AM_AC ), SqrtTiny )

        AP_VR     =  AP_I ( iV )  -  V_Dim_IR ( iV )
        AP_AC     =  AP_I ( iV )  -  AC_I ( iV )
  !      AP_AC_Inv =  1.0_KDR &
  !                   / sign ( max ( abs ( AP_AC ), SqrtTiny ), AP_AC )
        AP_AC_Inv =  1.0_KDR &
                     / max ( abs ( AP_AC ), SqrtTiny )

        DE_ICL ( iV )  =  DE_IL ( iV ) * AM_VL * AM_AC_Inv
        DE_ICR ( iV )  =  DE_IR ( iV ) * AP_VR * AP_AC_Inv

      end do !-- iV
      !$OMP  end parallel do
    
    end if


  end procedure ComputeCenterStatesKernel


end submodule Fluid_P_HN__Kernel
