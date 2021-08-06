#include "Preprocessor"

submodule ( Universe_F_CC__Form ) Universe_F_CC__Kernel

  use Basics
  
  implicit none

contains


  module procedure Compute_dT_G_CGS_Kernel

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
      
    nV  =  size ( dX_1 )
    
    dT_Inverse  =  - huge ( 0.0_KDR )

    select case ( nDimensions )
    case ( 1 )

      !dT_Inverse &
      !  = maxval ( sqrt ( abs ( GradPhi_1 ) / ( dX_1 ) ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                       sqrt ( abs ( M_UU_11 ( iV )  *  GradPhi_1 ( iV ) ) &
                              / dX_1 ( iV ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              =  max ( dT_Inverse, &
                       sqrt ( abs ( M_UU_11 ( iV )  *  GradPhi_1 ( iV ) ) &
                              / dX_1 ( iV ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 2 )

      !dT_Inverse &
      !  = maxval ( sqrt (    abs ( GradPhi_1 ) &
      !                       / ( dX_1 ) &
      !                    +  abs ( M_UU_22 * GradPhi_2 ) &
      !                       / ( Crsn2 * dX_2 ) ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              = max ( dT_Inverse, &
                      sqrt (    abs ( M_UU_11 ( iV )  *  GradPhi_1 ( iV ) ) &
                                /  dX_1 ( iV ) &
                             +  abs ( M_UU_22 ( iV )  *  GradPhi_2 ( iV ) ) &
                                /  dX_2 ( iV ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              = max ( dT_Inverse, &
                      sqrt (    abs ( M_UU_11 ( iV )  *  GradPhi_1 ( iV ) ) &
                                /  dX_1 ( iV ) &
                             +  abs ( M_UU_22 ( iV )  *  GradPhi_2 ( iV ) ) &
                                /  dX_2 ( iV ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 3 )

      !dT_Inverse &
      !  = maxval ( sqrt (    abs ( M_UU_11 * GradPhi_1 ) &
      !                       / ( dX_1 ) &
      !                    +  abs ( M_UU_22 * GradPhi_2 ) &
      !                       / ( Crsn2 * dX_2 ) ), &
      !                    +  abs ( M_UU_11 * GradPhi_2 ) &
      !                       / ( Crsn2 * dX_2 ) ), &
      !             mask = ProperCell )      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              = max ( dT_Inverse, &
                      sqrt (   abs ( M_UU_11 ( iV ) * GradPhi_1 ( iV ) ) &
                               /  dX_1 ( iV ) &
                             + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                               /  dX_2 ( iV ) &
                             + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
                               /  dX_3 ( iV ) ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : dT_Inverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            dT_Inverse &
              = max ( dT_Inverse, &
                      sqrt (   abs ( M_UU_11 ( iV ) * GradPhi_1 ( iV ) ) &
                               /  dX_1 ( iV ) &
                             + abs ( M_UU_22 ( iV ) * GradPhi_2 ( iV ) ) &
                               /  dX_2 ( iV ) &
                             + abs ( M_UU_33 ( iV ) * GradPhi_3 ( iV ) ) &
                               /  dX_3 ( iV ) ) )
        end do
        !$OMP  end parallel do
      end if
      
    end select !-- nDimensions

    dT_Inverse  =  max ( tiny ( 0.0_KDR ), dT_Inverse )
    dT          =  min ( dT, 1.0_KDR  /  dT_Inverse )

  end procedure Compute_dT_G_CGS_Kernel

  
end submodule Universe_F_CC__Kernel
