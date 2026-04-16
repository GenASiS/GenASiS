!-- Search_Command implements an overloaded "search" subroutine to find the 
!   index of a sorted array (of intrinsic types) corresponding to an input 
!   value

#include "Preprocessor"

module Search_Command

  !-- based on Numerical Recipes, Fortran (1992), Fortran 90 (1996)  
  !   routine "locate", but with modification
  
  !-- Given a sorted array A and a value Value, returns iValue such that 
  !   A ( iValue ) <= Value < A ( iValue + 1 )

  use Specifiers

  implicit none
  private

  public :: &
    Search
    
  interface Search
    module procedure SearchInteger
    module procedure SearchReal
  end interface Search

contains


  pure subroutine SearchInteger ( A, Value, iValue ) 
    
    integer ( KDI ), dimension ( : ), intent ( in )  :: &
      A
    integer ( KDI ), intent ( in )  :: &
      Value
    integer ( KDI ), intent ( out )  :: &
      iValue
    
    integer ( KDI )  :: &
      nValues, &
      iLow, &
      iMiddle, &
      iHigh
    logical ( KDL )  :: &
      Ascending
      
    !OMP_DECLARE_TARGET
      
    nValues = size ( A ) 
    Ascending = ( A ( nValues )  >= A ( 1 )  ) 
    iLow = 0
    iHigh = nValues + 1
    do while ( iHigh - iLow > 1 ) 
      iMiddle =  ( iHigh + iLow )  / 2
      if ( Ascending .and. ( Value >= A ( iMiddle )  )  ) then
        iLow = iMiddle
      else
        iHigh = iMiddle
      end if
    end do
    if ( Value == A ( 1 ) ) then
      iValue = 1
    else if ( Value == A ( nValues ) ) then
      iValue = nValues  !-- modification from Numerical Recipes
    else
      iValue = iLow
    end if
  
  end subroutine SearchInteger


  pure subroutine SearchReal ( A, Value, iValue ) 
  
    real ( KDR ), dimension ( : ), intent ( in )  :: &
      A
    real ( KDR ), intent ( in )  :: &
      Value
    integer ( KDI ), intent ( out )  :: &
      iValue
    
    integer ( KDI )  :: &
      nValues, &
      iLow, &
      iMiddle, &
      iHigh
    logical ( KDL )  :: &
      Ascending
      
    !OMP_DECLARE_TARGET
      
    nValues = size ( A ) 
    Ascending = ( A ( nValues )  >= A ( 1 )  ) 
    iLow = 0
    iHigh = nValues + 1
    do while ( iHigh - iLow > 1 ) 
      iMiddle =  ( iHigh + iLow )  / 2
      if ( Ascending .and. ( Value >= A ( iMiddle )  )  ) then
        iLow = iMiddle
      else
        iHigh = iMiddle
      end if
    end do
    if ( Value == A ( 1 ) ) then
      iValue = 1
    else if ( Value == A ( nValues ) ) then
      iValue = nValues  !-- modification from Numerical Recipes
    else
      iValue = iLow
    end if
  
  end subroutine SearchReal


!--  pure subroutine Search ( A, Value, iValue ) 
!--
!--    class ( * ), dimension ( : ), intent ( in )  :: &
!--      A
!--    class ( * ), intent ( in )  :: &
!--      Value
!--    integer ( KDI ), intent ( out )  :: &
!--      iValue
!--    
!--    integer ( KDI )  :: &
!--      nValues, &
!--      iLow, &
!--      iMiddle, &
!--      iHigh
!--    logical ( KDL )  :: &
!--      Ascending
!--      
!--    !OMP_DECLARE_TARGET
!--      
!--    select type ( A )
!--    
!--    type is ( integer ( KDI ) )
!--      
!--      select type ( Value )
!--      type is ( integer ( KDI ) )
!--        nValues = size ( A ) 
!--        Ascending = ( A ( nValues )  >= A ( 1 )  ) 
!--        iLow = 0
!--        iHigh = nValues + 1
!--        do while ( iHigh - iLow > 1 ) 
!--          iMiddle =  ( iHigh + iLow )  / 2
!--          if ( Ascending .and. ( Value >= A ( iMiddle )  )  ) then
!--            iLow = iMiddle
!--          else
!--            iHigh = iMiddle
!--          end if
!--        end do
!--        if ( Value == A ( 1 ) ) then
!--          iValue = 1
!--        else if ( Value == A ( nValues ) ) then
!--          iValue = nValues  !-- modification from Numerical Recipes
!--        else
!--          iValue = iLow
!--        end if
!--      end select
!--    
!--    type is ( real ( KDR ) )
!--      
!--      select type ( Value )
!--      type is ( real ( KDR ) )
!--        nValues = size ( A ) 
!--        Ascending = ( A ( nValues )  >= A ( 1 )  ) 
!--        iLow = 0
!--        iHigh = nValues + 1
!--        do while ( iHigh - iLow > 1 ) 
!--          iMiddle =  ( iHigh + iLow )  / 2
!--          if ( Ascending .and. ( Value >= A ( iMiddle )  )  ) then
!--            iLow = iMiddle
!--          else
!--            iHigh = iMiddle
!--          end if
!--        end do
!--        if ( Value == A ( 1 ) ) then
!--          iValue = 1
!--        else if ( Value == A ( nValues ) ) then
!--          iValue = nValues  !-- modification from Numerical Recipes
!--        else
!--          iValue = iLow
!--        end if
!--      end select
!--    
!--    end select
!--    
!--  end subroutine Search


end module Search_Command
