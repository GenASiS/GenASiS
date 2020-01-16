#include "Preprocessor"

submodule ( Fluid_D__Form ) Fluid_D__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_M_Kernel 

    !-- Compute_BaryonMass_Kernel

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
        M ( iV )  =  BaryonMassReference
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        M ( iV )  =  BaryonMassReference
      end do !-- iV
      !$OMP  end parallel do
    end if

  end procedure Compute_M_Kernel


  module procedure Compute_D_S_G_Kernel
 	 
    !-- Compute_ConservedDensity_Momentum_Galilean_Kernel

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
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do

      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV ) * N ( iV ) * V_1 ( iV )	 	 
        S_2 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_22 ( iV ) * V_2 ( iV )
        S_3 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_33 ( iV ) * V_3 ( iV )

      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else 

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( N ( iV )  <  0.0_KDR ) &
          N ( iV )  =  0.0_KDR
      end do !-- iV
      !$OMP  end parallel do

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues

        D ( iV ) = N ( iV ) 	 	 
       
        S_1 ( iV )  =  M ( iV ) * N ( iV ) * V_1 ( iV )	 	 
        S_2 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_22 ( iV ) * V_2 ( iV )
        S_3 ( iV )  =  M ( iV ) * N ( iV ) * M_DD_33 ( iV ) * V_3 ( iV )

      end do !-- iV
      !$OMP  end parallel do

    end if

  end procedure Compute_D_S_G_Kernel 	 	 


  module procedure Compute_N_V_G_Kernel

    !-- Compute_ComovingBaryonDensity_Velocity_Galilean

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
        if ( D ( iV )  >  N_Min ) then
          N ( iV )    =  D ( iV )
          V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
        else
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV )
      do iV = 1, nValues
        if ( D ( iV )  >  N_Min ) then
          N ( iV )    =  D ( iV )
          V_1 ( iV )  =  S_1 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_2 ( iV )  =  M_UU_22 ( iV ) * S_2 ( iV ) / ( M ( iV ) * D ( iV ) )
          V_3 ( iV )  =  M_UU_33 ( iV ) * S_3 ( iV ) / ( M ( iV ) * D ( iV ) )
        else
          N   ( iV )  =  N_Min
          V_1 ( iV )  =  0.0_KDR
          V_2 ( iV )  =  0.0_KDR
          V_3 ( iV )  =  0.0_KDR
          D   ( iV )  =  N_Min
          S_1 ( iV )  =  0.0_KDR
          S_2 ( iV )  =  0.0_KDR
          S_3 ( iV )  =  0.0_KDR
        end if
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_N_V_G_Kernel


  module procedure Compute_FE_D_G_Kernel

    !-- Compute_FastEigenspeeds_Dust_Galilean_Kernel

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
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        FEP_1 ( iV )  =  V_1 ( iV ) 
        FEP_2 ( iV )  =  V_2 ( iV ) 
        FEP_3 ( iV )  =  V_3 ( iV ) 
        FEM_1 ( iV )  =  V_1 ( iV ) 
        FEM_2 ( iV )  =  V_2 ( iV ) 
        FEM_3 ( iV )  =  V_3 ( iV ) 
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else
    
      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        FEP_1 ( iV )  =  V_1 ( iV ) 
        FEP_2 ( iV )  =  V_2 ( iV ) 
        FEP_3 ( iV )  =  V_3 ( iV ) 
        FEM_1 ( iV )  =  V_1 ( iV ) 
        FEM_2 ( iV )  =  V_2 ( iV ) 
        FEM_3 ( iV )  =  V_3 ( iV ) 
      end do !-- iV
      !$OMP  end parallel do

    end if

  end procedure Compute_FE_D_G_Kernel


  module procedure ComputeRawFluxes_G_Kernel

    !-- Galilean

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
        F_D   ( iV )  =  D   ( iV )  *  V_Dim ( iV ) 
        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV ) 
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV ) 
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
      do iV = 1, nValues
        F_D   ( iV )  =  D   ( iV )  *  V_Dim ( iV ) 
        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV ) 
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV ) 
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV ) 
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure ComputeRawFluxes_G_Kernel


  module procedure ComputeRawFluxes_N_S_Kernel

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      One_FourPi_G
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    One_FourPi_G  =  1.0_KDR  /  ( 4.0_KDR  *  CONSTANT % PI  &
                                   *  CONSTANT % GRAVITATIONAL )

    nValues = size ( F_S_1 )

    select case ( iDimension )
    case ( 1 )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  GradPhi_1 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  GradPhi_1 ( iV ) 
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  GradPhi_1 ( iV ) 
        end do !-- iV
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  GradPhi_1 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  GradPhi_1 ( iV ) 
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  GradPhi_1 ( iV ) 
        end do !-- iV
        !$OMP  end parallel do
      end if
    case ( 2 )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) 
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV )
        end do !-- iV
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) 
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  M_UU_22 ( iV )  *  GradPhi_2 ( iV )
        end do !-- iV
        !$OMP  end parallel do
      end if
    case ( 3 )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) 
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV )
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
        end do !-- iV
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE ) private ( iV ) 
        do iV = 1, nValues
          F_S_1 ( iV )  &
            =  F_S_1 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_1 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) 
          F_S_2 ( iV )  &
            =  F_S_2 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_2 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV )
          F_S_3 ( iV )  &
            =  F_S_3 ( iV )  &
               +  One_FourPi_G  &
                  *  GradPhi_3 ( iV )  *  M_UU_33 ( iV )  *  GradPhi_3 ( iV ) &
               -  One_FourPi_G  *  0.5_KDR  &
                  *  (    GradPhi_1 ( iV ) * GradPhi_1 ( iV ) &
                       +  GradPhi_2 ( iV ) * M_UU_22 ( iV ) * GradPhi_2 ( iV ) &
                       +  GradPhi_3 ( iV ) * M_UU_33 ( iV ) * GradPhi_3 ( iV ) )
        end do !-- iV
        !$OMP  end parallel do
      end if
    end select !-- iDimension

  end procedure ComputeRawFluxes_N_S_Kernel


end submodule Fluid_D__Kernel
