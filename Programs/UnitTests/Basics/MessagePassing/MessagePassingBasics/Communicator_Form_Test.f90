program Communicator_Form_Test

  use Specifiers
  use Display
  use Communicator_Form

  implicit none

  integer ( KDI ) :: &
    iR  !-- iRank
  type ( CommunicatorForm ), allocatable :: &
    C, &
    SC_S, &   !-- Single
    SC_FH, &  !-- FirstHalf
    SC_SH     !-- SecondHalf

  allocate ( C )
  call C % Initialize ( )
  call CONSOLE % Initialize ( C % Rank )
  call CONSOLE % SetVerbosity ( 'INFO_2' )
!  call CONSOLE % SetDisplayRank ( 0 )
  call CONSOLE % SetDisplayRank ( C % Size - 1 )
  
  allocate ( SC_S )
  call SC_S % Initialize ( C, [ C % Rank ], NameOption = 'Single' )
  
  allocate ( SC_FH )
  call SC_FH % Initialize &
         ( C, [ ( iR, iR = 0, C % Size / 2 - 1 ) ], &
           NameOption = 'FirstHalf' )

  allocate ( SC_SH )
  call SC_SH % Initialize &
         ( C, [ ( iR, iR = C % Size / 2, C % Size - 1 ) ], &
           NameOption = 'SecondHalf' )

  deallocate ( SC_SH )
  deallocate ( SC_FH )
  deallocate ( SC_S )
  deallocate ( C )

end program Communicator_Form_Test
