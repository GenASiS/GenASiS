!-- Sort_Commands implements an overloaded "Sort" subroutine to sort arrays
!   of intrinsict types

module Sort_Command

  use Specifiers
  
  implicit none 
  private
  
  public :: & 
    Sort
    
  interface Sort
    module procedure Sort_I
    module procedure Sort_R
  end interface Sort 
  
  private :: &
    PartitionIndex_I, &
    PartitionIndex_R
  
contains

  !-- Quick sort, based on http://www.fortran.com/qsort_c.f95
  !   by Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03, 
  !   based in turn on algorithm from Cormen et al., 
  !   Introduction to Algorithms, 1997 printing 
  
  recursive pure subroutine Sort_I ( A ) 

    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      A
    
    integer ( KDI ) :: &
      iPartition
  
    if ( size ( A ) > 1 ) then
      call PartitionIndex_I ( A, iPartition )
      call Sort ( A ( : iPartition - 1 ) )
      call Sort ( A ( iPartition : ) ) 
    end if
    
  end subroutine Sort_I
  
  
  recursive pure subroutine Sort_R ( A ) 

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      A
    
    integer ( KDI ) :: &
      iPartition
  
    if ( size ( A ) > 1 ) then
      call PartitionIndex_R ( A, iPartition )
      call Sort ( A ( : iPartition - 1 ) )
      call Sort ( A ( iPartition : ) ) 
    end if
    
  end subroutine Sort_R


  pure subroutine PartitionIndex_I ( A, iPartition )
    
    integer ( KDI ), dimension ( : ), intent ( inout )  :: &
      A
    integer ( KDI ), intent ( out ) :: &
      iPartition

    integer ( KDI ) :: &
      iLow, iHigh
    integer ( KDI ) :: &
      Swap, &
      Pivot
      
    Pivot = A ( 1 ) 
    iLow  = 0
    iHigh = size ( A ) + 1
    
    do

      iHigh = iHigh - 1
      do while ( A ( iHigh ) > Pivot ) 
        iHigh = iHigh - 1
      end do

      iLow = iLow + 1
      do while ( A ( iLow ) < Pivot ) 
        iLow = iLow + 1
      end do

      if ( iLow < iHigh ) then
        Swap        = A ( iLow ) 
        A ( iLow )  = A ( iHigh ) 
        A ( iHigh ) = Swap
      else if ( iLow == iHigh ) then
        iPartition = iLow + 1
        return
      else
        iPartition = iLow
        return
      end if

    end do
      
  end subroutine PartitionIndex_I
  
  
  pure subroutine PartitionIndex_R ( A, iPartition )
    
    real ( KDR ), dimension ( : ), intent ( inout )  :: &
      A
    integer ( KDI ), intent ( out ) :: &
      iPartition

    integer ( KDI ) :: &
      iLow, iHigh
    real ( KDR ) :: &
      Swap, &
      Pivot
      
    Pivot = A ( 1 ) 
    iLow  = 0
    iHigh = size ( A ) + 1
    
    do

      iHigh = iHigh - 1
      do while ( A ( iHigh ) > Pivot ) 
        iHigh = iHigh - 1
      end do

      iLow = iLow + 1
      do while ( A ( iLow ) < Pivot ) 
        iLow = iLow + 1
      end do

      if ( iLow < iHigh ) then
        Swap        = A ( iLow ) 
        A ( iLow )  = A ( iHigh ) 
        A ( iHigh ) = Swap
      else if ( iLow == iHigh ) then
        iPartition = iLow + 1
        return
      else
        iPartition = iLow
        return
      end if

    end do

      
  end subroutine PartitionIndex_R
  

end module Sort_Command 
