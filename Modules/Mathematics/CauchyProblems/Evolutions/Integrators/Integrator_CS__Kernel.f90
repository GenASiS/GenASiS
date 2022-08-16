#include "Preprocessor"

submodule ( Integrator_CS__Form ) Integrator_CS__Kernel

  use Basics
  
  implicit none

contains


  module procedure Compute_dT_CS_CGS_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      dT_Inverse
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( FEP_1 )
    
    dT_Inverse  =  - huge ( 0.0_KDR )

    select case ( nDimensions )
    case ( 1 )

      !dT_Inverse &
      !  = maxval ( max ( FEP_1, -FEM_1 ) / ( dX_1 ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse ) MAP_DT_INVERSE 
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                       max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                       /  dX_1 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                       max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                       /  dX_1 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 2 )

      !dT_Inverse &
      !  = maxval (   max ( FEP_1, -FEM_1 ) / ( dX_1 ) &
      !             + max ( FEP_2, -FEM_2 ) / ( Crsn_2 * dX_2 ) ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse ) MAP_DT_INVERSE
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                          max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                          /  dX_1 ( iV ) &
                       +  max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                          /  dX_2 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                          max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                          /  dX_1 ( iV ) &
                       +  max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                          /  dX_2 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 3 )

      ! dT_Inverse &
      !   = maxval (   max ( FEP_1, -FEM_1 ) / ( dX_1 ) &
      !              + max ( FEP_2, -FEM_2 ) / ( Crsn2 * dX_2 ) &
      !              + max ( FEP_3, -FEM_3 ) / ( Crsn3 * dX_3 ), &
      !              mask = ProperCell )
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse ) MAP_DT_INVERSE
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                          max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                           /  dX_1 ( iV ) &
                       +  max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                          /  dX_2 ( iV ) &
                       +  max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                          /  dX_3 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                          max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                          /  dX_1 ( iV ) &
                       +  max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                          /  dX_2 ( iV ) &
                       +  max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                          /  dX_3 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    end select !-- nDimensions

    dT_Inverse  =  max ( tiny ( 0.0_KDR ), dT_Inverse )
    dT          =  min ( dT, 1.0_KDR  /  dT_Inverse )

  end procedure Compute_dT_CS_CGS_Kernel

  
end submodule Integrator_CS__Kernel
