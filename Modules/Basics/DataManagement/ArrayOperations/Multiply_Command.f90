!-- Multiply provides an overloaded interface to add matrices, in order to 
!   expose elemental variables to the compiler and include threading.

#include "Preprocessor"

module Multiply_Command

  use Specifiers
  
  public :: &
    Multiply
    
  interface Multiply
    module procedure MultiplyReal_0D_2D_InPlace
  end interface


contains  


  subroutine MultiplyReal_0D_2D_InPlace ( A, B, UseDeviceOption )
  
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      A
    real ( KDR ), intent ( in ) :: &
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
          A ( iV, jV )  =  B  *  A ( iV, jV )
        end do
      end do
      !$OMP end OMP_TARGET_DIRECTIVE parallel do
    else 
      !$OMP parallel do collapse ( 2 ) &
      !$OMP schedule ( OMP_SCHEDULE_HOST )
      do jV  =  1,  nV ( 2 )
        do iV  =  1,  nV ( 1 )
          A ( iV, jV )  =  B  *  A ( iV, jV )
        end do
      end do
      !$OMP end parallel do
    end if
    
  end subroutine MultiplyReal_0D_2D_InPlace


end module Multiply_Command
