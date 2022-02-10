#include "Preprocessor"

submodule ( DivergencePart_F_P_V__Form ) DivergencePart_F_P_V__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure Compute_FS_G_Kernel

    !-- Compute_FluxSet_Velocity_Galileo_Kernel

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

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )

        F_G ( iV )  =  G ( iV )  *  V_Dim ( iV )

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )

        F_G ( iV )  =  G ( iV )  *  V_Dim ( iV )

      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure Compute_FS_G_Kernel


  module procedure Compute_S_UD_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( S_UD_22 )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_S_UD_Kernel


end submodule DivergencePart_F_P_V__Kernel
