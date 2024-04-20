#include "Preprocessor"

module RelativeDifference_Command

  use Specifiers
  
  public :: &
    RelativeDifference

  interface RelativeDifference
    module procedure RelativeDifferenceSelf_1D
    module procedure RelativeDifferenceSelf_2D
    module procedure RelativeDifferenceScale_1D
  end interface RelativeDifference


contains


  subroutine RelativeDifferenceSelf_1D ( A, B, C, UseDeviceOption )
  
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
                     / max ( abs ( B ( iV ) ), abs ( A ( iV ) ), SqrtTiny )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV ( 1 )
        C ( iV )  =  abs ( B ( iV ) -  A ( iV ) ) &
                     / max ( abs ( B ( iV ) ), abs ( A ( iV ) ), SqrtTiny )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine RelativeDifferenceSelf_1D


  subroutine RelativeDifferenceSelf_2D ( A, B, C, UseDeviceOption )
  
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
                           / max ( abs ( B ( iV, jV ) ), abs ( A ( iV, jV ) ), &
                                   SqrtTiny )
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          C ( iV, jV )  =  abs ( B ( iV, jV ) -  A ( iV, jV ) ) &
                           / max ( abs ( B ( iV, jV ) ), abs ( A ( iV, jV ) ), &
                                   SqrtTiny )
        end do
      end do
      !$OMP end parallel do
    end if
    
  end subroutine RelativeDifferenceSelf_2D


  subroutine RelativeDifferenceScale_1D ( A, B, C, D, UseDeviceOption )
  
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      C
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      D
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
        D ( iV )  =  abs ( B ( iV ) -  A ( iV ) ) &
                     / max ( abs ( C ( iV ) ), SqrtTiny )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV ( 1 )
        D ( iV )  =  abs ( B ( iV ) -  A ( iV ) ) &
                     / max ( abs ( C ( iV ) ), SqrtTiny )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine RelativeDifferenceScale_1D


end module RelativeDifference_Command
