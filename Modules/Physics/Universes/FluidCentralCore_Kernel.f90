#include "Preprocessor"

submodule ( FluidCentralCore_Form ) FluidCentralCore_Kernel

  use Basics
  
  implicit none

contains


  module procedure ComputeTimeStepKernel_G_CSL

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      TimeStepInverse
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    nV = size ( dXL_1 )
    
    TimeStepInverse = - huge ( 0.0_KDR )

    select case ( nDimensions )
    case ( 1 )
      
      !TimeStepInverse &
      !  = maxval ( sqrt ( abs ( GradPhi_1 ) / ( dXL_1 + dXR_1 ) ), &
      !             mask = IsProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV 
          if ( IsProperCell ( iV ) ) &
            cycle
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt ( abs ( GradPhi_1 ( iV ) ) &
                       / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV 
          if ( IsProperCell ( iV ) ) &
            cycle
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt ( abs ( GradPhi_1 ( iV ) ) &
                       / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 2 )
      
      !TimeStepInverse &
      !  = maxval ( sqrt (    abs ( GradPhi_1 ) &
      !                       / ( dXL_1 + dXR_1 ) &
      !                    +  abs ( M_UU_22 * GradPhi_2 ) &
      !                       / ( dXL_2 + dXR_2 ) ), &
      !             mask = IsProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt ( abs ( GradPhi_1 ( iV ) ) &
                             / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                           + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                             / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt (   abs ( GradPhi_1 ( iV ) ) &
                                   / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                           + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                                   / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 3 )
      
      ! TimeStepInverse &
      !   = maxval ( sqrt (   abs ( GradPhi_1 ) / dX_1 &
      !                     + abs ( M_UU_22 * GradPhi_2 ) / dX_2 &
      !                     + abs ( M_UU_33 * GradPhi_3 ) / dX_3 ) ), &
      !              mask = IsProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt (   abs ( GradPhi_1 ( iV ) ) &
                                   / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                           + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                                   / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) &
                           + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
                                   / ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          
          TimeStepInverse &
            = max ( TimeStepInverse, &
                    sqrt (   abs ( GradPhi_1 ( iV ) ) &
                                   / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                           + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                                   / ( dXL_2 ( iV ) + dXR_2 ( iV ) ) &
                           + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
                                   / ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
        end do
        !$OMP  end parallel do
      end if
    
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end procedure ComputeTimeStepKernel_G_CSL
  
  
end submodule FluidCentralCore_Kernel
