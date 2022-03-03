!-- Add provides an overloaded interface to add matrices, in order to 
!   expose elemental variables to the compiler and include threading.

#include "Preprocessor"

module Add_Command

  use Specifiers
  
  public :: &
    Add
    
  interface Add
    module procedure AddReal_1D
    module procedure AddReal_1D_InPlace
    module procedure AddReal_2D
    module procedure AddReal_2D_InPlace
    module procedure AddReal_3D
    module procedure AddReal_3D_InPlace
    ! module procedure AddSectionReal_1D
    ! module procedure AddSectionReal_1D_InPlace
    ! module procedure AddReal_2D_InPlace
  end interface

contains  


  subroutine AddReal_1D ( A, B, C, UseDeviceOption )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      A, &
      B
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      C
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
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV  =  1,  nV
        C ( iV )  =  A ( iV )  +  B ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV
        C ( iV )  =  A ( iV )  +  B ( iV )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine AddReal_1D


  subroutine AddReal_1D_InPlace ( A, B, UseDeviceOption )
  
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
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
      !$OMP OMP_TARGET_DIRECTIVE parallel do &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do iV  =  1,  nV
        A ( iV )  =  A ( iV )  +  B ( iV )
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do iV  =  1,  nV
        A ( iV )  =  A ( iV )  +  B ( iV )
      end do
      !$OMP end parallel do
    end if
    
  end subroutine AddReal_1D_InPlace


  subroutine AddReal_2D ( A, B, C, UseDeviceOption )
  
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
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV  =  shape ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          C ( iV, jV )  =  A ( iV, jV )  +  B ( iV, jV )
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          C ( iV, jV )  =  A ( iV, jV )  +  B ( iV, jV )
        end do
      end do
      !$OMP end parallel do
    end if
    
  end subroutine AddReal_2D


  subroutine AddReal_2D_InPlace ( A, B, UseDeviceOption )
  
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
                      
    integer ( KDI ) :: &
      iV, jV
    integer ( KDI ), dimension ( 2 ) :: &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV  =  shape ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          A ( iV, jV )  =  A ( iV, jV )  +  B ( iV, jV )
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          A ( iV, jV )  =  A ( iV, jV )  +  B ( iV, jV )
        end do
      end do
      !$OMP end parallel do
    end if
    
  end subroutine AddReal_2D_InPlace


  subroutine AddReal_3D ( A, B, C, UseDeviceOption )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      A, &
      B
    real ( KDR ), dimension ( :, :, : ), intent ( out ) :: &
      C
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
    
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV  =  shape ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do kV  =  1,  nV ( 3 )
        do jV  =  1,  nV ( 2 )
          do iV  =  1,  nV ( 1 )
            C ( iV, jV, kV )  =  A ( iV, jV, kV )  +  B ( iV, jV, kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV  =  1,  nV ( 3 )
        do jV  =  1,  nV ( 2 )
          do iV  =  1,  nV ( 1 )
            C ( iV, jV, kV )  =  A ( iV, jV, kV )  +  B ( iV, jV, kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    end if
    
  end subroutine AddReal_3D


  subroutine AddReal_3D_InPlace ( A, B, UseDeviceOption )
  
    real ( KDR ), dimension ( :, :, : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      B
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
                      
    integer ( KDI ) :: &
      iV, jV, kV
    integer ( KDI ), dimension ( 3 ) :: &
      nV
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nV  =  shape ( A )

    if ( UseDevice ) then
      !$OMP OMP_TARGET_DIRECTIVE parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_TARGET )
      do kV  =  1,  nV ( 3 )
        do jV  =  1,  nV ( 2 )
          do iV  =  1,  nV ( 1 )
            A ( iV, jV, kV )  =  A ( iV, jV, kV )  +  B ( iV, jV, kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 3 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do kV  =  1,  nV ( 3 )
        do jV  =  1,  nV ( 2 )
          do iV  =  1,  nV ( 1 )
            A ( iV, jV, kV )  =  A ( iV, jV, kV )  +  B ( iV, jV, kV )
          end do !-- iV
        end do !-- jV
      end do !-- kV
      !$OMP end parallel do
    end if

  end subroutine AddReal_3D_InPlace


  ! subroutine AddSectionReal_1D &
  !              ( A, B, oSource_A, oSource_B, oTarget, nValues, C )
  
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     A, &
  !     B
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource_A, &
  !     oSource_B, &
  !     oTarget, &
  !     nValues
  !   real ( KDR ), dimension ( : ), intent ( out ) :: &
  !     C

  !   call Add &
  !          ( A ( oSource_A + 1 : oSource_A + nValues ), &
  !            B ( oSource_B + 1 : oSource_B + nValues ), &
  !            C ( oTarget + 1 : oTarget + nValues ) )                      
  
  ! end subroutine AddSectionReal_1D
                            

  ! subroutine AddSectionReal_1D_InPlace ( A, B, oSource, oTarget, nValues )
  
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     A
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     B
  !   integer ( KDI ), intent ( in ) :: &
  !     oSource, &
  !     oTarget, &
  !     nValues

  !   call Add &
  !          ( A ( oSource + 1 : oSource + nValues ), &
  !            B ( oTarget + 1 : oTarget + nValues ) )                      
  
  ! end subroutine AddSectionReal_1D_InPlace
                            

  ! subroutine AddReal_2D_InPlace ( A, B )
  
  !   real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
  !     A
  !   real ( KDR ), dimension ( :, : ), intent ( in ) :: &
  !     B
                      
  !   !$OMP parallel workshare
  !   A = A + B
  !   !$OMP end parallel workshare
  
  ! end subroutine AddReal_2D_InPlace
                            

end module Add_Command
