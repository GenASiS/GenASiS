#include "Preprocessor"

submodule ( DiffusionFactor_RM__Form ) DiffusionFactor_RM__Kernel
  
  use Basics
  
  implicit none
  
contains


  module procedure ComputeKernel

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice      
          
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    nV  =  size ( DF )

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) &
      !$OMP shared ( SqrtTiny )
      do iV  =  1,  nV
        DF ( iV )  &
          =  SF ( iV )  &
             /  max ( sqrt ( M_DD ( iV ) )  *  dX ( iV )  *  TO ( iV ),  &
                      SqrtTiny )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) &
      !$OMP shared ( SqrtTiny )
      do iV  =  1,  nV
        DF ( iV )  &
          =  SF ( iV )  &
             /  max ( sqrt ( M_DD ( iV ) )  *  dX ( iV )  *  TO ( iV ),  &
                      SqrtTiny )
      end do
      !$OMP end parallel do
    end if

  end procedure ComputeKernel


end submodule DiffusionFactor_RM__Kernel
