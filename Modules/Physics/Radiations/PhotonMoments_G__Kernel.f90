#include "Preprocessor"

submodule ( PhotonMoments_G__Form ) PhotonMoments_G__Kernel

  use Basics 
  
  implicit none

contains


  module procedure Compute_SP_A_Kernel

    !-- Compute_SpectralParameters_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      a
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J )

    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( a )
      do iV = 1, nV
        T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( a )
      do iV = 1, nV
        T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_SP_A_Kernel


  module procedure Compute_SP_S_Kernel

    !-- Compute_SpectralParameters_Single_Kernel

    real ( KDR ) :: &
      a

    !OMP_DECLARE_TARGET
      
    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    T_R ( iV )  =  ( J ( iV )  /  a ) ** ( 0.25_KDR )

  end procedure Compute_SP_S_Kernel


  module procedure Compute_Eq_A_Kernel

    !-- Compute_Equilibrium_All_Kernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny, &
      a
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( J_Eq )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny, a )
      do iV = 1, nV
        J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
        J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                         /  max ( SqrtTiny, J_Eq ( iV ) )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny, a )
      do iV = 1, nV
        J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
        J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                         /  max ( SqrtTiny, J_Eq ( iV ) )
      end do
      !$OMP end parallel do
    end if

  end procedure Compute_Eq_A_Kernel


  module procedure Compute_Eq_S_Kernel

    !-- Compute_Equilibrium_Single_Kernel

    real ( KDR ) :: &
      SqrtTiny, &
      a
    
    !OMP_DECLARE_TARGET
      
    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    a  =  4.0_KDR  *  CONSTANT % STEFAN_BOLTZMANN

    J_Eq  ( iV )  =  a  *  T ( iV ) ** 4
    J_RD  ( iV )  =  abs ( J ( iV )  -  J_Eq ( iV ) )  &
                     /  max ( SqrtTiny, J_Eq ( iV ) )

  end procedure Compute_Eq_S_Kernel


end submodule PhotonMoments_G__Kernel
