#include "Preprocessor"

submodule ( Step_RK_C_ASC__Template ) Step_RK_C_ASC__Kernel

  use Basics
  
  implicit none

contains


  module procedure RecordDivergenceKernel

    integer ( KDI ) :: &
      iV, &
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( SDV )
    
    if ( UseDevice ) then
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) 
      do iV = 1, nValues
        SDV ( iV )  =  SDV ( iV )  +  Weight_RK  *  IDV ( iV )  /  TimeStep
      end do !-- iV
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) 
      do iV = 1, nValues
        SDV ( iV )  =  SDV ( iV )  +  Weight_RK  *  IDV ( iV )  /  TimeStep
      end do !-- iV
      !$OMP  end parallel do
    
    end if

  end procedure RecordDivergenceKernel


end submodule Step_RK_C_ASC__Kernel
