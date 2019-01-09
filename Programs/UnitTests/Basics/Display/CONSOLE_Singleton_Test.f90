program CONSOLE_Singleton_Test

  use Specifiers
  use CONSOLE_Singleton

  implicit none

  include 'mpif.h'

  integer ( KDI ) :: &
    Rank, &
    Error

  interface
    subroutine ShowMessage &
                 ( Character, IgnorabilityOption, &
                   DisplayRankOption, nLeadingLinesOption, &
                   nTrailingLinesOption )
      use Specifiers
      character ( * ), intent ( in ) :: &
        Character
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption, &
        DisplayRankOption, &
        nLeadingLinesOption, &
        nTrailingLinesOption
    end subroutine ShowMessage
  end interface

  call MPI_INIT ( Error )
  call MPI_COMM_RANK ( MPI_COMM_WORLD, Rank, Error )

  call CONSOLE % Initialize ( ProcessRank = Rank )
  call CONSOLE % SetVerbosity ( 'INFO_7' )
  call CONSOLE % SetDisplayRank ( 0 )

  call CONSOLE % Mute ( )
  call ShowMessage ( 'This should not appear' )

  call CONSOLE % Unmute ( )
  call ShowMessage ( 'This should appear' )

  call MPI_FINALIZE ( Error )

end program CONSOLE_Singleton_Test
