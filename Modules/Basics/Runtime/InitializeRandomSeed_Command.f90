!-- InitializeRandomSeed initialize Fortran random number generator via 
!   Unix's pseudorandom generator /dev/urandom.

module InitializeRandomSeed_Command
  
  use VariableManagement
  use Display
  use MessagePassing
  use FileSystem

  implicit none
  private
  
  public :: &
    InitializeRandomSeed
    
    character ( LDF ), private, parameter :: &
      RANDOM_FILE = '/dev/urandom'
                                 
contains


  subroutine InitializeRandomSeed ( C )
      
    type ( CommunicatorForm ), intent ( in ) :: &
      C 

    integer ( KDI ) :: &
      FileUnit, &
      Status, &
      SeedSize
    integer ( KDI ), dimension ( : ), allocatable :: &
      Seed
    real ( KDR ) :: &
      TestRandom

    call random_seed ( size = SeedSize )
    allocate ( Seed ( SeedSize ) )

    call DelayFileAccess ( C % Rank )

    open ( newunit = FileUnit, file = RANDOM_FILE, access = 'stream', &
           form = 'UNFORMATTED', iostat = Status )
    read ( FileUnit ) Seed
    close ( FileUnit )

    call random_seed ( put = Seed )
    call random_number ( TestRandom )
    
    call Show ( Seed, 'Random seed', CONSOLE % INFO_4 )
    call Show ( TestRandom, 'TestRandom', CONSOLE % INFO_4 )

  end subroutine InitializeRandomSeed


end module InitializeRandomSeed_Command
