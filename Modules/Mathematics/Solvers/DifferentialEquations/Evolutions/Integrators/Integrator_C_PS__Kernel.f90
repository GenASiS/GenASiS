#include "Preprocessor"

submodule ( Integrator_C_PS__Form ) Integrator_C_PS__Kernel

  use Basics
  
  implicit none

contains

  
  module procedure ComputeTimeStepKernel_CSL

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
      
    nV = size ( FEP_1 )
    
    TimeStepInverse = - huge ( 0.0_KDR )

    select case ( nDimensions )
    case ( 1 )

      !TimeStepInverse &
      !  = maxval ( max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ), &
      !             mask = IsProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 2 )

      !TimeStepInverse &
      !  = maxval (   max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ) &
      !             + max ( FEP_2, -FEM_2 ) / ( Crsn_2 * ( dXL_2 + dXR_2 ) ), &
      !             mask = IsProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / ( Crsn_2 ( iV ) * ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / ( Crsn_2 ( iV ) * ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 3 )
      ! TimeStepInverse &
      !   = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
      !              + max ( FEP_2, -FEM_2 ) / dX_2 &
      !              + max ( FEP_3, -FEM_3 ) / dX_3, &
      !              mask = IsProperCell )
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / ( Crsn_2 ( iV ) * ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) &
                      + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                        / ( Crsn_3 ( iV ) * ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( IsProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / ( dXL_1 ( iV ) + dXR_1 ( iV ) ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / ( Crsn_2 ( iV ) * ( dXL_2 ( iV ) + dXR_2 ( iV ) ) ) &
                      + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                        / ( Crsn_3 ( iV ) * ( dXL_3 ( iV ) + dXR_3 ( iV ) ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end procedure ComputeTimeStepKernel_CSL
  
end submodule Integrator_C_PS__Kernel
