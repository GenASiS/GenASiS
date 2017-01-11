!-- DelayFileAccess provides a mechanism to delay  operaration 
!   (e.g. open, read, write) based on process rank, suitable to prevent 
!   e.g. metadata operation congestion on parallel file system.

module DelayFileAccess_Command

  use VariableManagement

  implicit none
  private
  
  public :: &
    DelayFileAccess

    integer ( KDI ), private, parameter :: &
      N_PROCESSES_FILE_ACCESS = 128, &
      FILE_ACCESS_DELAY_UNIT  = 10000

contains


  subroutine DelayFileAccess ( ProcessRank )

    integer ( KDI ), intent ( in ) :: &
      ProcessRank

    integer ( KDI ) :: &
      iDelay, &
      nFileAccessDelayIterations
    real ( KDR ) :: &
      Number, &
      Dummy
      
    nFileAccessDelayIterations &
      =  ( ProcessRank / N_PROCESSES_FILE_ACCESS ) * FILE_ACCESS_DELAY_UNIT
    
    do iDelay = 1, nFileAccessDelayIterations
      !-- wait my turn
      call random_number ( Number )
      Dummy = Number * Number
    end do

  end subroutine DelayFileAccess


end module DelayFileAccess_Command
