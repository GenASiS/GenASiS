#include "Preprocessor"

submodule ( Fluid_D__Form ) Fluid_D__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_M_Kernel 

    !-- Compute_BaryonMass_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( M )
    
    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        M ( iV )  =  M_Ref
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        M ( iV )  =  M_Ref
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure Compute_M_Kernel


  module procedure Compute_D_S_G_Kernel
 	 
    !-- Compute_ConservedDensity_Momentum_Galileo_Kernel

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

        if ( N ( iV )  <  N_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
        end if

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else 

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

        if ( N ( iV )  <  N_Min ) then
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
        end if

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_11 ( iV )  *  V_1 ( iV )
        S_2 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_22 ( iV )  *  V_2 ( iV )
        S_3 ( iV )  =  M ( iV )  *  N ( iV )  *  M_DD_33 ( iV )  *  V_3 ( iV )

      end do !-- iV
      !$OMP end parallel do

    end if

  end procedure Compute_D_S_G_Kernel 	 	 


  module procedure Compute_N_V_G_Kernel

    !-- Compute_ComovingBaryonDensity_Velocity_Galileo

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

        if ( D ( iV )  <  N_Min ) then
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
        end if

call Show ( M ( iV ), '>>> M' )
call Show ( D ( iV ), '>>> D' )
        N ( iV )    =  D ( iV )
        V_1 ( iV )  =  M_UU_11 ( iV )  &
                       *  S_1 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_2 ( iV )  =  M_UU_22 ( iV )  &
                       *  S_2 ( iV )  /  ( M ( iV )  *  D ( iV ) )
        V_3 ( iV )  =  M_UU_33 ( iV )  &
                       *  S_3 ( iV )  /  ( M ( iV )  *  D ( iV ) )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
      
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

        if ( D ( iV )  <  N_Min ) then
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

      end do !-- iV
      !$OMP end parallel do
    
    end if

  end procedure Compute_N_V_G_Kernel


  module procedure ComputeFluxes_G_Kernel

    !-- ComputeFluxes_Galileo_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( F_D )
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        F_D   ( iV )  =  D   ( iV )  *  V_Dim ( iV ) 
        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV ) 
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV ) 
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        F_D   ( iV )  =  D   ( iV )  *  V_Dim ( iV ) 
        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV ) 
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV ) 
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure ComputeFluxes_G_Kernel


  module procedure ComputeEigenspeeds_G_Kernel

    !-- ComputeEigenspeeds_Galileo_Kernel

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
        EF_P ( iV )  =  V_Dim ( iV )
        EF_M ( iV )  =  V_Dim ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        EF_P ( iV )  =  V_Dim ( iV )
        EF_M ( iV )  =  V_Dim ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeEigenspeeds_G_Kernel


end submodule Fluid_D__Kernel
