!-- Clear_Command provides overloaded "clear" routines with intent(out)
!   arguments of intrinsic data types, so that the compiler will use fast
!   operations to set arrays to zero (or to .false. in the case of logical 
!   arrays)

#include "Preprocessor"

module Clear_Command

  use Specifiers
  use Devices

  implicit none
  private

  public :: &
    Clear

  interface Clear
    module procedure ClearInteger_1D
    module procedure ClearInteger_2D
    module procedure ClearInteger_3D
    ! module procedure ClearInteger_4D
    ! module procedure ClearBigInteger_1D
    ! module procedure ClearBigInteger_2D
    ! module procedure ClearBigInteger_3D
    ! module procedure ClearBigInteger_4D
    module procedure ClearReal_1D
    module procedure ClearReal_2D
    module procedure ClearReal_3D
    ! module procedure ClearReal_4D
    ! module procedure ClearComplex_1D
    ! module procedure ClearComplex_2D
    module procedure ClearComplex_3D
    ! module procedure ClearComplex_4D
    module procedure ClearLogical_1D
    ! module procedure ClearLogical_2D
    ! module procedure ClearLogical_3D
    ! module procedure ClearLogical_4D
    ! module procedure ClearTinyLogical_1D
    ! module procedure ClearTinyLogical_2D
  end interface Clear

contains


  subroutine ClearInteger_1D ( A )

    integer ( KDI ), dimension ( : ), intent ( out ) :: &
      A

    !$OMP parallel workshare
    A = 0_KDI
    !$OMP end parallel workshare

  end subroutine ClearInteger_1D

  
  subroutine ClearInteger_2D ( A )

    integer ( KDI ), dimension ( :, : ), intent ( out ) :: &
      A

    !$OMP parallel workshare
    A = 0_KDI
    !$OMP end parallel workshare

  end subroutine ClearInteger_2D


  subroutine ClearInteger_3D ( A )

    integer ( KDI ), dimension ( :, :, : ), intent ( out ) :: &
      A

    !$OMP parallel workshare
    A = 0_KDI
    !$OMP end parallel workshare

  end subroutine ClearInteger_3D


  ! subroutine ClearInteger_4D ( A )

  !   integer ( KDI ), dimension ( :, :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = 0_KDI
  !   !$OMP end parallel workshare

  ! end subroutine ClearInteger_4D


  ! subroutine ClearBigInteger_1D ( A )

  !   integer ( KBI ), dimension ( : ), intent ( out ) :: &
  !     A
    
  !   !$OMP parallel workshare
  !   A = 0_KBI
  !   !$OMP end parallel workshare

  ! end subroutine ClearBigInteger_1D

  
  ! subroutine ClearBigInteger_2D ( A )

  !   integer ( KBI ), dimension ( :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = 0_KBI
  !   !$OMP end parallel workshare

  ! end subroutine ClearBigInteger_2D


  ! subroutine ClearBigInteger_3D ( A )

  !   integer ( KBI ), dimension ( :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = 0_KBI
  !   !$OMP end parallel workshare

  ! end subroutine ClearBigInteger_3D


  ! subroutine ClearBigInteger_4D ( A )

  !   integer ( KBI ), dimension ( :, :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = 0_KBI
  !   !$OMP end parallel workshare

  ! end subroutine ClearBigInteger_4D


  subroutine ClearReal_1D ( A, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( out ) :: &
      A
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV = size ( A )
    
    if ( UseDevice ) then
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP& schedule ( OMP_SCHEDULE )
      do iV = 1, nV
        A ( iV ) = 0.0_KDR
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do private ( iV ) schedule ( OMP_SCHEDULE )
      do iV = 1, nV
        A ( iV ) = 0.0_KDR
      end do
      !$OMP end parallel do
    end if

  end subroutine ClearReal_1D
  
  
  subroutine ClearReal_2D ( A, UseDeviceOption )

    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      A
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV
      
    nV = size ( A, dim = 2 )

    do iV = 1, nV
      call Clear ( A ( :, iV ), UseDeviceOption )
    end do
  
  end subroutine ClearReal_2D
  
  
  subroutine ClearReal_3D ( A )

    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      A

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( A )

    !$OMP parallel do private ( iV, jV, kV ) schedule ( OMP_SCHEDULE )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          A ( iV, jV, kV ) = 0.0_KDR
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ClearReal_3D
  
  
  ! subroutine ClearReal_4D ( A )

  !   real ( KDR ), dimension ( :, :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = 0.0_KDR
  !   !$OMP end parallel workshare

  ! end subroutine ClearReal_4D
  
  
  ! subroutine ClearComplex_1D ( A )

  !   complex ( KDC ), dimension ( : ), intent ( out ) :: &
  !      A

  !   !$OMP parallel workshare
  !    A = 0.0_KDC
  !   !$OMP end parallel workshare

  ! end subroutine ClearComplex_1D


  ! subroutine ClearComplex_2D ( A )

  !   complex ( KDC ), dimension ( :, : ), intent ( out ) :: &
  !      A

  !   !$OMP parallel workshare
  !    A = 0.0_KDC
  !   !$OMP end parallel workshare

  ! end subroutine ClearComplex_2D


  subroutine ClearComplex_3D ( A )

    complex ( KDC ), dimension ( :, :, : ), intent ( out ) :: &
       A

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( A )

    !$OMP parallel do private ( iV, jV, kV ) schedule ( OMP_SCHEDULE )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          A ( iV, jV, kV ) = 0.0_KDC
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine ClearComplex_3D


  ! subroutine ClearComplex_4D ( A )

  !   complex ( KDC ), dimension ( :, :, :, : ), intent ( out ) :: &
  !      A

  !   !$OMP parallel workshare
  !    A = 0.0_KDC
  !   !$OMP end parallel workshare

  !  end subroutine ClearComplex_4D


  subroutine ClearLogical_1D ( A )

    logical ( KDL ), dimension ( : ), intent ( out ) :: &
      A

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( A )

    !$OMP parallel do private ( iV ) schedule ( OMP_SCHEDULE )
    do iV = 1, nV
      A ( iV ) = .false.
    end do
    !$OMP end parallel do

  end subroutine ClearLogical_1D
  
  
  ! subroutine ClearLogical_2D ( A )

  !   logical ( KDL ), dimension ( :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = .false.
  !   !$OMP end parallel workshare

  ! end subroutine ClearLogical_2D
  
  
  ! subroutine ClearLogical_3D ( A )

  !   logical ( KDL ), dimension ( :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = .false.
  !   !$OMP end parallel workshare

  ! end subroutine ClearLogical_3D
  
  
  ! subroutine ClearLogical_4D ( A )

  !   logical ( KDL ), dimension ( :, :, :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = .false.
  !   !$OMP end parallel workshare

  ! end subroutine ClearLogical_4D
  
  
  ! subroutine ClearTinyLogical_1D ( A )

  !   logical ( KTL ), dimension ( : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = .false.
  !   !$OMP end parallel workshare

  ! end subroutine ClearTinyLogical_1D
  
  
  ! subroutine ClearTinyLogical_2D ( A )

  !   logical ( KTL ), dimension ( :, : ), intent ( out ) :: &
  !     A

  !   !$OMP parallel workshare
  !   A = .false.
  !   !$OMP end parallel workshare

  ! end subroutine ClearTinyLogical_2D
  
  
end module Clear_Command
