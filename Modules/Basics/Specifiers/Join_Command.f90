!-- Join_Command implements subroutine Join to join an array of strings 
!   by glue string

module Join_Command
  
  use KIND_DEFAULT_Singleton
  use KIND_BIG_Singleton
  
  implicit none
  private
  
  public :: &
    Join, &
    Join_KBCH

  interface Join
    module procedure Join_KDCH
!    module procedure Join_KBCH
  end interface !-- Join
  
contains


  pure subroutine Join_KDCH ( Piece, Glue, String )
  
    character ( *, KDCH ), dimension ( : ), intent ( in ) :: &
      Piece
    character ( *, KDCH ), intent ( in ) :: &
      Glue
    character ( *, KDCH ), intent ( out ) :: &
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
    
  end subroutine Join_KDCH


  pure subroutine Join_KBCH ( Piece, Glue, String )
  
    character ( *, KBCH ), dimension ( : ), intent ( in ) :: &
      Piece
    character ( *, KBCH ), intent ( in ) :: &
      Glue
    character ( *, KBCH ), intent ( out ) :: &
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
    
  end subroutine Join_KBCH


end module Join_Command
