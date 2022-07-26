!-- Copy_Command provides overloaded "copy" routines with intent(in) and
!   intent(out) arguments of intrinsic data types, so that the compiler will
!   use fast copy methods. 

#include "Preprocessor"

module Copy_Command

  use Specifiers
  use Devices

  implicit none
  private

  public :: &
    Copy, &
    CopyCollapse

  interface Copy
    module procedure CopyInteger_1D
    ! module procedure CopyInteger_2D
    ! module procedure CopyInteger_3D
    ! module procedure CopyInteger_1D_Section
    ! module procedure CopyInteger_2D_Section
    ! module procedure CopyInteger_3D_Section
    module procedure CopyBigInteger_1D
    ! module procedure CopyBigInteger_2D
    ! module procedure CopyBigInteger_3D
    ! module procedure CopyBigInteger_1D_Section
    ! module procedure CopyBigInteger_2D_Section
    ! module procedure CopyBigInteger_3D_Section
    module procedure CopyReal_1D
    module procedure CopyReal_2D
    module procedure CopyReal_3D
!     module procedure CopyReal_1D_Section
!     module procedure CopyReal_2D_Section
!     module procedure CopyReal_3D_Section
! !-- NOTE: This may be needed if KDR is not double precision
! !    module procedure Copy_C_Double_Real_1D
    module procedure CopyReal_3D_1D
    module procedure CopyReal_1D_3D
    module procedure CopyComplex_1D
!     module procedure CopyComplex_2D
!     module procedure CopyComplex_3D
!     module procedure CopyComplex_1D_Section
!     module procedure CopyComplex_2D_Section
!     module procedure CopyComplex_3D_Section
!     module procedure Copy_Character_C_String
  end interface Copy

  interface CopyCollapse
    module procedure CopyCollapseReal_3D
    module procedure CopyCollapseReal_3D_Offset
  end interface CopyCollapse

contains


  subroutine CopyInteger_1D ( A, B )

    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      A
    integer ( KDI ), dimension ( : ), intent ( out ) :: &
      B

    !$OMP parallel workshare
    B = A
    !$OMP end parallel workshare

  end subroutine CopyInteger_1D
  
  
  ! subroutine CopyInteger_2D ( A, B )

  !   integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyInteger_2D
  
  
  ! subroutine CopyInteger_3D ( A, B )

  !   integer ( KDI ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( :, :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyInteger_3D
  
  
  ! subroutine CopyInteger_1D_Section ( A, oSource, oTarget, nValues, B )
    
  !   integer ( KDI ), dimension ( : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KDI ), dimension ( : ), intent ( inout ) :: &
  !     B

  !   call Copy ( A ( oSource + 1 : oSource + nValues ), &
  !               B ( oTarget + 1 : oTarget + nValues ) )

  ! end subroutine CopyInteger_1D_Section
  

  ! subroutine CopyInteger_2D_Section ( A, oSource, oTarget, nValues, B )

  !   integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KDI ), dimension ( :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2  !-- iValue_2

  !   do iV_2 = 1, nValues ( 2 )
  !     call Copy (  &
  !            A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                oSource ( 2 ) + iV_2 ), &
  !            B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                oTarget ( 2 ) + iV_2 ) )
  !   end do

  ! end subroutine CopyInteger_2D_Section
  

  ! subroutine CopyInteger_3D_Section ( A, oSource, oTarget, nValues, B )

  !   integer ( KDI ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KDI ), dimension ( :, :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2, iV_3  !-- iValue_2, iValue_3

  !   do iV_3 = 1, nValues ( 3 )
  !     do iV_2 = 1, nValues ( 2 )
  !       call Copy (  &
  !              A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2, oSource ( 3 ) + iV_3 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2, oTarget ( 3 ) + iV_3 ) )
  !     end do
  !   end do

  ! end subroutine CopyInteger_3D_Section
  

  subroutine CopyBigInteger_1D ( A, B )

    integer ( KBI ), dimension ( : ), intent ( in ) :: &
      A
    integer ( KBI ), dimension ( : ), intent ( out ) :: &
      B

    !$OMP parallel workshare
    B = A
    !$OMP end parallel workshare

  end subroutine CopyBigInteger_1D
  
  
  ! subroutine CopyBigInteger_2D ( A, B )

  !   integer ( KBI ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KBI ), dimension ( :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyBigInteger_2D
  
  
  ! subroutine CopyBigInteger_3D ( A, B )

  !   integer ( KBI ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KBI ), dimension ( :, :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyBigInteger_3D
  
  
  ! subroutine CopyBigInteger_1D_Section ( A, oSource, oTarget, nValues, B )
    
  !   integer ( KBI ), dimension ( : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KBI ), dimension ( : ), intent ( inout ) :: &
  !     B

  !   call Copy ( A ( oSource + 1 : oSource + nValues ), &
  !               B ( oTarget + 1 : oTarget + nValues ) )

  ! end subroutine CopyBigInteger_1D_Section
  

  ! subroutine CopyBigInteger_2D_Section ( A, oSource, oTarget, nValues, B )

  !   integer ( KBI ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KBI ), dimension ( :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2  !-- iValue_2

  !   do iV_2 = 1, nValues ( 2 )
  !     call Copy (  &
  !            A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                oSource ( 2 ) + iV_2 ), &
  !            B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                oTarget ( 2 ) + iV_2 ) )
  !   end do

  ! end subroutine CopyBigInteger_2D_Section
  

  ! subroutine CopyBigInteger_3D_Section ( A, oSource, oTarget, nValues, B )

  !   integer ( KBI ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   integer ( KBI ), dimension ( :, :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2, iV_3  !-- iValue_2, iValue_3

  !   do iV_3 = 1, nValues ( 3 )
  !     do iV_2 = 1, nValues ( 2 )
  !       call Copy (  &
  !              A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2, oSource ( 3 ) + iV_3 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2, oTarget ( 3 ) + iV_3 ) )
  !     end do
  !   end do

  ! end subroutine CopyBigInteger_3D_Section
  

  subroutine CopyReal_1D ( A, B, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV
    logical ( KDL ) :: &
      UseDevice

    nV = size ( A )
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption 
    
    if ( UseDevice ) then 
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nV
        B ( iV ) = A ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else
      !$OMP parallel do private ( iV ) schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nV
        B ( iV ) = A ( iV )
      end do
      !$OMP end parallel do
    end if

  end subroutine CopyReal_1D
  
  
  subroutine CopyReal_2D ( A, B, UseDeviceOption )

    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( :, : ), intent ( out ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( A, dim = 2 )

    do iV = 1, nV
      call Copy ( A ( :, iV ), B ( :, iV ), UseDeviceOption )
    end do
  
  end subroutine CopyReal_2D
  
  
  subroutine CopyReal_3D ( A, B, UseDeviceOption )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    logical ( KDL ) :: &
      UseDevice

    nV = shape ( A )
    
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption 

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP private ( iV, jV, kV ) schedule ( OMP_SCHEDULE_TARGET )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            B ( iV, jV, kV ) = A ( iV, jV, kV )
          end do
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 ) & 
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV = 1, nV ( 3 )
        do jV = 1, nV ( 2 )
          do iV = 1, nV ( 1 )
            B ( iV, jV, kV ) = A ( iV, jV, kV )
          end do
        end do
      end do
      !$OMP end parallel do
    end if

  end subroutine CopyReal_3D
  
  
  ! subroutine CopyReal_1D_Section ( A, oSource, oTarget, nValues, B )
    
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     B

  !   call Copy &
  !          ( A ( oSource + 1 : oSource + nValues ), &
  !            B ( oTarget + 1 : oTarget + nValues ) )
  
  ! end subroutine CopyReal_1D_Section
  

  ! subroutine CopyReal_2D_Section ( A, oSource, oTarget, nValues, B )

  !   real ( KDR ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2  !-- iValue_2

  !   do iV_2 = 1, nValues ( 2 )
  !     call Copy &
  !            ( A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2 ) )
  !   end do

  ! end subroutine CopyReal_2D_Section


  ! subroutine CopyReal_3D_Section ( A, oSource, oTarget, nValues, B )

  !   real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2, &  !-- iValue_2
  !     iV_3     !-- iValue_3

  !   do iV_3 = 1, nValues ( 3 )
  !     do iV_2 = 1, nValues ( 2 )
  !     call Copy &
  !            ( A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2, oSource ( 3 ) + iV_3 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2, oTarget ( 3 ) + iV_3 ) )
  !     end do
  !   end do

  ! end subroutine CopyReal_3D_Section


  ! ! subroutine Copy_C_Double_Real_1D ( A, B )
    
  ! !   real ( c_double ), dimension ( : ), intent ( in ) :: &
  ! !     A
  ! !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  ! !     B
      
  ! !   B = A
  
  ! ! end subroutine Copy_C_Double_Real_1D

  
  subroutine CopyReal_3D_1D &
               ( A, nSource, oSource, oTarget, B, UseDeviceOption )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      A
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nSource, &
      oSource
    integer ( KDI ), intent ( in ) :: &
      oTarget
    !-- This argument is last in the spirit of intent ( out ), but should
    !   remain intent ( inout ) so existing values outside the section
    !   are not corrupted
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    integer ( KDI ) :: &
      iT, &
      iV, jV, kV, &
      iS, jS, kS
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
      
    if ( UseDevice ) then
      
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iS, jS, kS, iT )
      do kV = oSource ( 3 ) + 1, oSource ( 3 ) + nSource ( 3 )
        do jV = oSource ( 2 ) + 1, oSource ( 2 ) + nSource ( 2 )
          do iV = oSource ( 1 ) + 1, oSource ( 1 ) + nSource ( 1 )
            
            iS = iV - oSource ( 1 )
            jS = jV - oSource ( 2 )
            kS = kV - oSource ( 3 )
            
            iT = iS + ( jS - 1 ) * nSource ( 1 ) &
                    + ( kS - 1 ) * nSource ( 1 ) * nSource ( 2 )
            
            B ( oTarget + iT ) = A ( iV, jV, kV )
            
          end do
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do

    else

      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iS, jS, kS, iT )
      do kV = oSource ( 3 ) + 1, oSource ( 3 ) + nSource ( 3 )
        do jV = oSource ( 2 ) + 1, oSource ( 2 ) + nSource ( 2 )
          do iV = oSource ( 1 ) + 1, oSource ( 1 ) + nSource ( 1 )
            
            iS = iV - oSource ( 1 )
            jS = jV - oSource ( 2 )
            kS = kV - oSource ( 3 )
            
            iT = iS + ( jS - 1 ) * nSource ( 1 ) &
                    + ( kS - 1 ) * nSource ( 1 ) * nSource ( 2 )
            
            B ( oTarget + iT ) = A ( iV, jV, kV )
            
          end do
        end do
      end do
      !$OMP end parallel do
    
    end if

  end subroutine CopyReal_3D_1D


  subroutine CopyReal_1D_3D &
               ( A, nTarget, oTarget, oSource, B, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      A
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      nTarget, &
      oTarget
    integer ( KDI ), intent ( in ) :: &
      oSource
    !-- This argument is last in the spirit of intent ( out ), but should
    !   remain intent ( inout ) so existing values outside the section
    !   are not corrupted
    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
    
    integer ( KDI ) :: &
      iS, &
      iV, jV, kV, &
      iT, jT, kT
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption
    
    if ( UseDevice ) then
    
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET ) private ( iT, jT, kT, iS )
      do kV = oTarget ( 3 ) + 1, oTarget ( 3 ) + nTarget ( 3 )
        do jV = oTarget ( 2 ) + 1, oTarget ( 2 ) + nTarget ( 2 )
          do iV = oTarget ( 1 ) + 1, oTarget ( 1 ) + nTarget ( 1 )
          
            iT = iV - oTarget ( 1 )
            jT = jV - oTarget ( 2 )
            kT = kV - oTarget ( 3 )
            
            iS = iT + ( jT - 1 ) * nTarget ( 1 ) &
                    + ( kT - 1 ) * nTarget ( 1 ) * nTarget ( 2 )
            
            B ( iV, jV, kV ) = A ( oSource + iS )
          
          end do
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    
    else
      
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST ) private ( iT, jT, kT, iS )
      do kV = oTarget ( 3 ) + 1, oTarget ( 3 ) + nTarget ( 3 )
        do jV = oTarget ( 2 ) + 1, oTarget ( 2 ) + nTarget ( 2 )
          do iV = oTarget ( 1 ) + 1, oTarget ( 1 ) + nTarget ( 1 )
          
            iT = iV - oTarget ( 1 )
            jT = jV - oTarget ( 2 )
            kT = kV - oTarget ( 3 )
            
            iS = iT + ( jT - 1 ) * nTarget ( 1 ) &
                    + ( kT - 1 ) * nTarget ( 1 ) * nTarget ( 2 )
            
            B ( iV, jV, kV ) = A ( oSource + iS )
          
          end do
        end do
      end do
      !$OMP end parallel do
    
    end if

  end subroutine CopyReal_1D_3D


  subroutine CopyComplex_1D ( A, B )

    complex ( KDC ), dimension ( : ), intent ( in ) :: &
      A
    complex ( KDC ), dimension ( : ), intent ( out ) :: &
      B

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( A )

    !$OMP parallel do private ( iV ) schedule ( OMP_SCHEDULE_HOST )
    do iV = 1, nV
      B ( iV ) = A ( iV )
    end do
    !$OMP end parallel do

  end subroutine CopyComplex_1D
  
  
  ! subroutine CopyComplex_2D ( A, B )

  !   complex ( KDC ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   complex ( KDC ), dimension ( :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyComplex_2D
  
  
  ! subroutine CopyComplex_3D ( A, B )

  !   complex ( KDC ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   complex ( KDC ), dimension ( :, :, : ), intent ( out ) :: &
  !     B

  !   !$OMP parallel workshare
  !   B = A
  !   !$OMP end parallel workshare

  ! end subroutine CopyComplex_3D
  
  
  ! subroutine CopyComplex_1D_Section ( A, oSource, oTarget, nValues, B )
    
  !   complex ( KDC ), dimension ( : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   complex ( KDC ), dimension ( : ), intent ( inout ) :: &
  !     B

  !   call Copy &
  !          ( A ( oSource + 1 : oSource + nValues ), &
  !            B ( oTarget + 1 : oTarget + nValues ) )
  
  ! end subroutine CopyComplex_1D_Section
  

  ! subroutine CopyComplex_2D_Section ( A, oSource, oTarget, nValues, B )

  !   complex ( KDC ), dimension ( :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 2 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   complex ( KDC ), dimension ( :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2  !-- iValue_2

  !   do iV_2 = 1, nValues ( 2 )
  !     call Copy &
  !            ( A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2 ) )
  !   end do

  ! end subroutine CopyComplex_2D_Section


  ! subroutine CopyComplex_3D_Section ( A, oSource, oTarget, nValues, B )

  !   complex ( KDC ), dimension ( :, :, : ), intent ( in ) :: &
  !     A
  !   integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues
  !   !-- This argument is last in the spirit of intent ( out ), but should
  !   !   remain intent ( inout ) so existing values outside the section
  !   !   are not corrupted
  !   complex ( KDC ), dimension ( :, :, : ), intent ( inout ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iV_2, &  !-- iValue_2
  !     iV_3     !-- iValue_3

  !   do iV_3 = 1, nValues ( 3 )
  !     do iV_2 = 1, nValues ( 2 )
  !     call Copy &
  !            ( A ( oSource ( 1 ) + 1 : oSource ( 1 ) + nValues ( 1 ), &
  !                  oSource ( 2 ) + iV_2, oSource ( 3 ) + iV_3 ), &
  !              B ( oTarget ( 1 ) + 1 : oTarget ( 1 ) + nValues ( 1 ), &
  !                  oTarget ( 2 ) + iV_2, oTarget ( 3 ) + iV_3 ) )
  !     end do
  !   end do

  ! end subroutine CopyComplex_3D_Section


  ! subroutine Copy_Character_C_String ( A, B )
    
  !   character ( * ), intent ( in ) :: &
  !     A
  !   character ( c_char ), dimension ( : ), intent ( out ) :: &
  !     B
      
  !   integer ( KDI ) :: &
  !     iC
      
  !   do iC = 1, len_trim ( A ) 
  !     B ( iC ) = A ( iC : iC )
  !   end do 
  !   B ( iC ) = c_null_char
  
  ! end subroutine Copy_Character_C_String
  

  subroutine CopyCollapseReal_3D ( A, B )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      B

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( B )

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 ) &
    !$OMP schedule ( OMP_SCHEDULE_HOST )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          B ( iV, jV, kV ) = A ( iV, jV, kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine CopyCollapseReal_3D


  subroutine CopyCollapseReal_3D_Offset ( A, B, oV )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      A
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      B
    integer ( KDI ), dimension ( 3 ), intent ( in ) :: &
      oV

    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV

    nV = shape ( B )

    !$OMP parallel do private ( iV, jV, kV ) collapse ( 3 ) &
    !$OMP schedule ( OMP_SCHEDULE_HOST )
    do kV = 1, nV ( 3 )
      do jV = 1, nV ( 2 )
        do iV = 1, nV ( 1 )
          B ( iV, jV, kV ) = A ( oV ( 1 ) + iV, oV ( 2 ) + jV, oV ( 3 ) + kV )
        end do
      end do
    end do
    !$OMP end parallel do

  end subroutine CopyCollapseReal_3D_Offset


end module Copy_Command
