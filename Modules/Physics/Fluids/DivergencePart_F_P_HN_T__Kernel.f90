#include "Preprocessor"

submodule ( DivergencePart_F_P_HN_T__Form ) DivergencePart_F_P_HN_T__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure Compute_FS_G_Kernel

    !-- Compute_FluxSet_Velocity_Galileo_Kernel

    integer :: &
      Delta_1, Delta_2, Delta_3
    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( F_D )
    
    Delta_1  =  int (    real ( 1 + iDim  -  abs ( 1 - iDim ) )  &
                      /  real ( 1 + iDim  +  abs ( 1 - iDim ) ) )

    Delta_2  =  int (    real ( 2 + iDim  -  abs ( 2 - iDim ) )  &
                      /  real ( 2 + iDim  +  abs ( 2 - iDim ) ) )

    Delta_3  =  int (    real ( 3 + iDim  -  abs ( 3 - iDim ) )  &
                      /  real ( 3 + iDim  +  abs ( 3 - iDim ) ) )

    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3 )
      do iV = 1, nV

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )  +  Delta_1  *  P ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )  +  Delta_2  *  P ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )  +  Delta_3  *  P ( iV )

        F_G ( iV )  =  ( G ( iV )  +  P ( iV ) )  *  V_Dim ( iV )

        F_DE ( iV )  =  DE ( iV )  *  V_Dim ( iV ) 

      end do !-- iV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( Delta_1, Delta_2, Delta_3 )
      do iV = 1, nV

        F_D ( iV )  =  D ( iV )  *  V_Dim ( iV ) 

        F_S_1 ( iV )  =  S_1 ( iV )  *  V_Dim ( iV )  +  Delta_1  *  P ( iV )
        F_S_2 ( iV )  =  S_2 ( iV )  *  V_Dim ( iV )  +  Delta_2  *  P ( iV )
        F_S_3 ( iV )  =  S_3 ( iV )  *  V_Dim ( iV )  +  Delta_3  *  P ( iV )

        F_G ( iV )  =  ( G ( iV )  +  P ( iV ) )  *  V_Dim ( iV )

        F_DE ( iV )  =  DE ( iV )  *  V_Dim ( iV ) 

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
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )  +  P ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )  +  P ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        S_UD_22 ( iV )  =  V_2 ( iV )  *  S_2 ( iV )  +  P ( iV )
        S_UD_33 ( iV )  =  V_3 ( iV )  *  S_3 ( iV )  +  P ( iV )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_S_UD_Kernel


end submodule DivergencePart_F_P_HN_T__Kernel
