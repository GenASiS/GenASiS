#include "Preprocessor"

submodule ( Fluid_D__Form ) Fluid_D__Kernel

  use Basics
  
  implicit none
  
contains


  module procedure ComputeBaryonMassKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nValues = size ( M )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        M ( iV ) = 1.0_KDR
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        M ( iV ) = 1.0_KDR
      end do !-- iV
      !$OMP end parallel do
    end if

  end procedure ComputeBaryonMassKernel


  module procedure ComputeDensityMomentumKernel
 	 
 	 	integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( D )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV )  <  0.0_KDR ) &
          N ( iV )  =  0.0_KDR
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV ) * N ( iV ) * V_1 ( iV )
        S_2 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_22 ( iV ) * V_2 ( iV )
        S_3 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_33 ( iV ) * V_3 ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV )  <  0.0_KDR ) &
          N ( iV )  =  0.0_KDR
      end do !-- iV
      !$OMP end parallel do

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV ) * N ( iV ) * V_1 ( iV )	 	 
        S_2 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_22 ( iV ) * V_2 ( iV )
        S_3 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_33 ( iV ) * V_3 ( iV )

      end do !-- iV
      !$OMP end parallel do
      
    end if

  end procedure ComputeDensityMomentumKernel


  module procedure ComputeDensityVelocityKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( N )
    
    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( D ( iV )  >  0.0_KDR ) then
          N ( iV )    =  D ( iV )
          V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
        else
          N   ( iV ) = 0.0_KDR
          V_1 ( iV ) = 0.0_KDR
          V_2 ( iV ) = 0.0_KDR
          V_3 ( iV ) = 0.0_KDR
          D   ( iV ) = 0.0_KDR
          S_1 ( iV ) = 0.0_KDR
          S_2 ( iV ) = 0.0_KDR
          S_3 ( iV ) = 0.0_KDR
        end if
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( D ( iV )  >  0.0_KDR ) then
          N ( iV )    =  D ( iV )
          V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
        else
          N   ( iV ) = 0.0_KDR
          V_1 ( iV ) = 0.0_KDR
          V_2 ( iV ) = 0.0_KDR
          V_3 ( iV ) = 0.0_KDR
          D   ( iV ) = 0.0_KDR
          S_1 ( iV ) = 0.0_KDR
          S_2 ( iV ) = 0.0_KDR
          S_3 ( iV ) = 0.0_KDR
        end if
      end do !-- iV
      !$OMP end parallel do
      
    end if

  end procedure ComputeDensityVelocityKernel


  module procedure ComputeEigenspeedsKernel_D

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( FEP_1 )
    
    if ( UseDevice ) then

      !$OMP  OMP_TARGET_DIRECTIVE parallel do 
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        FEP_1 ( iV ) = V_1 ( iV ) 
        FEP_2 ( iV ) = V_2 ( iV ) 
        FEP_3 ( iV ) = V_3 ( iV ) 
        FEM_1 ( iV ) = V_1 ( iV ) 
        FEM_2 ( iV ) = V_2 ( iV ) 
        FEM_3 ( iV ) = V_3 ( iV ) 
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else
    
      !$OMP  parallel do 
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        FEP_1 ( iV ) = V_1 ( iV ) 
        FEP_2 ( iV ) = V_2 ( iV ) 
        FEP_3 ( iV ) = V_3 ( iV ) 
        FEM_1 ( iV ) = V_1 ( iV ) 
        FEM_2 ( iV ) = V_2 ( iV ) 
        FEM_3 ( iV ) = V_3 ( iV ) 
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure ComputeEigenspeedsKernel_D


  module procedure ComputeRawFluxesKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( F_D )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        F_D   ( iV ) = D   ( iV ) * V_Dim ( iV ) 
        F_S_1 ( iV ) = S_1 ( iV ) * V_Dim ( iV ) 
        F_S_2 ( iV ) = S_2 ( iV ) * V_Dim ( iV ) 
        F_S_3 ( iV ) = S_3 ( iV ) * V_Dim ( iV ) 
      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        F_D   ( iV ) = D   ( iV ) * V_Dim ( iV ) 
        F_S_1 ( iV ) = S_1 ( iV ) * V_Dim ( iV ) 
        F_S_2 ( iV ) = S_2 ( iV ) * V_Dim ( iV ) 
        F_S_3 ( iV ) = S_3 ( iV ) * V_Dim ( iV ) 
      end do !-- iV
      !$OMP end parallel do
    end if      

  end procedure ComputeRawFluxesKernel


end submodule Fluid_D__Kernel
