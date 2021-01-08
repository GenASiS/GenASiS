!-- InitializeRandomSeed initialize Fortran random number generator via 
!   Unix's pseudorandom generator /dev/urandom.

module InitializeRandomSeed_Command
  
  use Specifiers
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
      iS, &
      FileUnit, &
      Status, &
      SeedSize
    integer ( KDI ), dimension ( : ), allocatable :: &
      Seed
    real ( KDR ) :: &
      TestRandom
      
    call Show ( 'Initializing random seed', CONSOLE % INFO_1 )

    call random_seed ( size = SeedSize )
    allocate ( Seed ( SeedSize ) )

    call DelayFileAccess ( C % Rank )

    open ( newunit = FileUnit, file = RANDOM_FILE, access = 'stream', &
           form = 'UNFORMATTED', status = 'old', iostat = Status )
    if ( Status /= 0 ) & !-- Workaround for IBM XL compiler
      open ( newunit = FileUnit, file = RANDOM_FILE, access = 'sequential', &
           form = 'UNFORMATTED', status = 'old', iostat = Status )
    
    if ( Status == 0 ) then
      read ( FileUnit ) Seed
      close ( FileUnit )
    else
      do iS = 1, SeedSize
        call system_clock ( Seed ( iS ), count_rate = TestRandom )
      end do
      Seed = Seed * ( C % Rank + 1 )
    end if
    
    call random_seed ( put = Seed )
    call random_number ( TestRandom )
    
    call Show ( Seed, 'Random seed', CONSOLE % INFO_4 )
    call Show ( TestRandom, 'TestRandom', CONSOLE % INFO_1 )

  end subroutine InitializeRandomSeed


end module InitializeRandomSeed_Command
