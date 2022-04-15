#include "Preprocessor"

submodule ( Gravitation_N_UA__Form ) Gravitation_N_UA__Kernel

  use Basics
  implicit none
  
contains


  module procedure SolveKernel

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nV
    logical ( KDL ) :: &
      UseDevice
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( Phi )

    select case ( nD )
    case ( 1 )
    
      if ( UseDevice ) then
        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nV
          Phi ( iV )        =  A  *  X ( iV )
          GradPhi_1 ( iV )  =  A
        end do
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST )
        do iV = 1, nV
          Phi ( iV )        =  A  *  X ( iV )
          GradPhi_1 ( iV )  =  A
        end do
        !$OMP end parallel do
      end if

    case ( 2 )
      
      if ( UseDevice ) then
        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nV
          Phi ( iV )        =  A  *  Y ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  A
        end do
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST )
        do iV = 1, nV
          Phi ( iV )        =  A  *  Y ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  A
        end do
        !$OMP end parallel do
      end if

    case ( 3 )
      
      if ( UseDevice ) then
        !$OMP OMP_TARGET_DIRECTIVE parallel do &
        !$OMP schedule ( OMP_SCHEDULE_TARGET )
        do iV = 1, nV
          Phi ( iV )        =  A  *  Z ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  0.0_KDR
          GradPhi_3 ( iV )  =  A
        end do
        !$OMP end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP parallel do &
        !$OMP schedule ( OMP_SCHEDULE_HOST )        
        do iV = 1, nV
          Phi ( iV )        =  A  *  Z ( iV )
          GradPhi_1 ( iV )  =  0.0_KDR
          GradPhi_2 ( iV )  =  0.0_KDR
          GradPhi_3 ( iV )  =  A
        end do
        !$OMP end parallel do
      end if

    end select !-- nD

  end procedure SolveKernel


end submodule Gravitation_N_UA__Kernel
