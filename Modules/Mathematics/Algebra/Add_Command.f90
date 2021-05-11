!-- Add provides an overloaded interface to add matrices, in order to 
!   expose elemental variables to the compiler and include threading.

module Add_Command

  use Basics
  
  public :: &
    Add
    
  interface Add
    module procedure AddReal_1D
    module procedure AddReal_1D_InPlace
    ! module procedure AddSectionReal_1D
    ! module procedure AddSectionReal_1D_InPlace
    ! module procedure AddReal_2D_InPlace
  end interface

contains  


  subroutine AddReal_1D ( A, B, C )
  
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      A, &
      B
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      C
                      
    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( A )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      C ( iV )  =  A ( iV ) + B ( iV )
    end do
    !$OMP end parallel do
  
  end subroutine AddReal_1D
                            

  subroutine AddReal_1D_InPlace ( A, B )
  
    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      B
                      
    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( A )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      A ( iV )  =  A ( iV ) + B ( iV )
    end do
    !$OMP end parallel do
  
  end subroutine AddReal_1D_InPlace
                            

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
