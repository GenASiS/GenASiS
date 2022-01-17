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
        N ( iV )  =  max ( N ( iV ), N_Min )
        E ( iV )  =  max ( E ( iV ), E_Min )
        T ( iV )  =  max ( T ( iV ), T_Min )
        Y ( iV )  =  max ( Y ( iV ), Y_Min )

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

        M ( iV )  =  M_Ref
        N ( iV )  =  max ( N ( iV ), N_Min )
        E ( iV )  =  max ( E ( iV ), E_Min )
        T ( iV )  =  max ( T ( iV ), T_Min )
        Y ( iV )  =  max ( Y ( iV ), Y_Min )

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

        D ( iV )  =  N ( iV ) 	 	 
       
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

        D ( iV )  =  N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

        G ( iV )  =  E ( iV )  +  0.5_KDR * (    S_1 ( iV ) * V_1 ( iV )  &
                                              +  S_2 ( iV ) * V_2 ( iV )  &
                                              +  S_3 ( iV ) * V_3 ( iV ) )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_D_S_G_DE_G_Kernel 	 	 


end submodule Fluid_P_HN__Kernel
