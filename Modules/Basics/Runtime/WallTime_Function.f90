!-- Walltime return the current time in seconds since an arbitrary time
!   in the pass. Currently this wraps MPI_WTIME(). 

module WallTime_Function
  
  use MPI
  use Specifiers

  implicit none
  private
  
  public :: &
    WallTime
  
  interface WallTime
    module procedure WallTime_MPI
  end interface WallTime
  
contains


  function WallTime_MPI ( ) result ( WT )
    
    type ( MeasuredValueForm ) :: &
      WT
    
    call WT % Initialize ( 's', MPI_WTIME ( ) )
  
  end function WallTime_MPI


end module WallTime_Function
