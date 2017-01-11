!-- Join_Command implements subroutine Join to join an array of strings 
!   by glue string

module Join_Command
  
  use KIND_DEFAULT_Singleton
  
  implicit none
  private
  
  public :: &
    Join
    
contains


  pure subroutine Join ( Piece, Glue, String )
  
    character ( * ), dimension ( : ), intent ( in ) :: &
      Piece
    character ( * ), intent ( in ) :: &
      Glue
    character ( * ), intent ( out ) :: &
      String
    
    integer ( KDI ) :: &
      iP, &  !-- iPiece
      nPieces
      
    nPieces = size ( Piece )
    String = Piece ( 1 )
    do iP = 2, nPieces
      String &
        = trim ( String ) // Glue // trim ( adjustl ( Piece ( iP ) ) )
    end do
    
  end subroutine Join


end module Join_Command
