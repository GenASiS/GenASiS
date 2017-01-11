!-- Split_Command implements subroutine Split to split a string by 
!   delimiter string

module Split_Command

  use KIND_DEFAULT_Singleton
  
  implicit none
  private
  
  public :: &
    Split
    
contains


  pure subroutine Split ( String, Delimiter, Piece )
  
    character ( * ), intent ( in ) :: &
      String, &
      Delimiter
    character ( * ), intent ( out ), dimension ( : ), allocatable :: &
      Piece
    
    integer ( KDI ) :: &
      iNP, &  !-- iNextPosition
      iP, &   !-- iPiece
      oC, &   !-- oCharacter
      lD, &   !-- lenDelimiter
      nPieces
    
    nPieces = 1
    oC = 0
    do
      iNP = index ( trim ( String ( oC + 1 : ) ), Delimiter )
      if ( iNP > 0 ) then
        nPieces = nPieces + 1
      else
        exit  
      end if
      oC = oC + iNP   
    end do
    
    if ( allocated ( Piece ) ) deallocate ( Piece ) 
    allocate ( Piece ( nPieces ) )
    
    lD = len ( Delimiter )
    
    oC = 0
    do iP = 1, nPieces - 1
      iNP = index ( trim ( String ( oC + 1 : ) ), Delimiter )
      Piece ( iP ) = adjustl ( String ( oC + 1 : oC + iNP - 1 ) )
      oC = oC + iNP + lD - 1
    end do
    Piece ( nPieces ) = adjustl ( String ( oC + 1 : ) ) 
  
  end subroutine Split


end module Split_Command
