#include "Preprocessor"

module RelativeDifference_Command

  use Specifiers
  
  public :: &
    RelativeDifference

  interface RelativeDifference
    module procedure RelativeDifference_1D
    module procedure RelativeDifference_2D
  end interface RelativeDifference


contains


  subroutine RelativeDifference_1D ( A, B, C, UseDeviceOption )
  
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
                      
    integer ( KDI ) :: &
      iV
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    nV  =  size ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV  =  1,  nV ( 1 )
        C ( iV )  =  abs ( B ( iV ) -  A ( iV ) ) &
                     / max ( abs ( B ( iV ) ), SqrtTiny )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV ( 1 )
        C ( iV )  =  abs ( B ( iV ) -  A ( iV ) ) &
                     / max ( abs ( B ( iV ) ), SqrtTiny )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine RelativeDifference_1D


  subroutine RelativeDifference_2D ( A, B, C, UseDeviceOption )
  
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
                      
    integer ( KDI ) :: &
      iV, jV
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    real ( KDR ) :: &
      SqrtTiny
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    SqrtTiny  =  sqrt ( tiny ( 0.0_KDR ) )

    nV  =  shape ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          C ( iV, jV )  =  abs ( B ( iV, jV ) -  A ( iV, jV ) ) &
                           / max ( abs ( B ( iV, jV ) ), SqrtTiny )
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          C ( iV, jV )  =  abs ( B ( iV, jV ) -  A ( iV, jV ) ) &
                           / max ( abs ( B ( iV, jV ) ), SqrtTiny )
        end do
      end do
      !$OMP end parallel do
    end if
    
  end subroutine RelativeDifference_2D


end module RelativeDifference_Command
