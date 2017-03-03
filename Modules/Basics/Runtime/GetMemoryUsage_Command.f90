!-- GetMemoryUsage gets per process and across processes memory usage via
!   reading and parsing of memory info file.

module GetMemoryUsage_Command
  
  use MPI
  use VariableManagement
  use Display
  use MessagePassing

  implicit none
  private
  
  public :: &
    GetMemoryUsage
    
    integer ( KDI ), private, save :: &
      FileUnit = 99
    character ( 7 ), private, parameter :: &
      FORMAT_BUFFER = '(a1023)'  !-- must correspond to 
                                 !   LDB = LEN_DEFAULT % BUFFER
    character ( LDF ), private, parameter :: &
      MEMORY_INFO_FILE = '/proc/self/status'
                                 
contains


  subroutine GetMemoryUsage &
               ( HWM, RSS, Ignorability, C_Option, Max_HWM_Option, &
                 Min_HWM_Option, Mean_HWM_Option, Max_RSS_Option, &
                 Min_RSS_Option, Mean_RSS_Option )
    
    type ( MeasuredValueForm ), intent ( out ) :: &
      HWM, &
      RSS
    integer ( KDI ), intent ( in ) :: &
      Ignorability
    type ( CommunicatorForm ), intent ( in ), optional :: &
      C_Option
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      Max_HWM_Option, &
      Min_HWM_Option, &
      Mean_HWM_Option, &
      Max_RSS_Option, &
      Min_RSS_Option, &
      Mean_RSS_Option
      
    integer ( KDI ) :: &
      Status
    real ( KDR ), dimension ( : ), pointer :: &
      Incoming
    logical ( KDL ) :: &
      MemoryInfoFileExists
    character ( LDB ) :: &
      Buffer
    type ( CollectiveOperation_R_Form ) :: &
      CO

    call HWM % Initialize ( 'kB', 0.0_KDR )
    call RSS % Initialize ( 'kB', 0.0_KDR )
    if ( present ( Max_HWM_Option ) ) &
      call Max_HWM_Option % Initialize ( 'kB', 0.0_KDR )    
    if ( present ( Min_HWM_Option ) ) &
      call Min_HWM_Option % Initialize ( 'kB', 0.0_KDR )
    if ( present ( Mean_HWM_Option ) ) &
      call Mean_HWM_Option % Initialize ( 'kB', 0.0_KDR )    
    if ( present ( Max_RSS_Option ) ) &
      call Max_RSS_Option % Initialize ( 'kB', 0.0_KDR )
    if ( present ( Min_RSS_Option ) ) &
      call Min_RSS_Option % Initialize ( 'kB', 0.0_KDR )    
    if ( present ( Mean_RSS_Option ) ) &
      call Mean_RSS_Option % Initialize ( 'kB', 0.0_KDR )
    
    inquire ( file = MEMORY_INFO_FILE, exist = MemoryInfoFileExists )
    if ( .not. MemoryInfoFileExists ) then
      call Show &
             ( 'Memory info file does not exist on this system', &
               Ignorability )
      return
    end if
     
    open ( FileUnit, file = MEMORY_INFO_FILE, iostat = Status )
    
    if ( Status /= 0 ) then
      call Show ( 'Failure in opening memory info file', Ignorability )
      return
    end if
     
    do while ( .true. )
      read ( FileUnit, '(a)', iostat = Status ) Buffer
      if ( Buffer ( : 6 ) == 'VmHWM:' ) &
        read ( Buffer ( 8 : ), * ) HWM % Number
      if ( Buffer ( : 6 ) == 'VmRSS:' ) &
        read ( buffer ( 8 : ), * ) RSS % Number 
      if ( Status < 0 ) exit 
    end do
     
    close ( FileUnit )
    
    !-- across processes value :
    
    if ( present ( C_Option ) &
         .and. ( present ( Max_HWM_Option ) .or. present ( Min_HWM_Option ) &
         .or. present ( Mean_HWM_Option ) .or. present ( Max_RSS_Option ) &
         .or. present ( Min_RSS_Option ) .or. present (  Mean_RSS_Option ) ) ) &
    then
      call CO % Initialize &
             ( C_Option, nOutgoing = [ 1 ], nIncoming = [ 1 ], &
               RootOption = CONSOLE % DisplayRank )
    end if
      
    if ( present ( Max_HWM_Option ) ) then
      CO % Outgoing % Value = [ HWM % Number ]
      call CO % Reduce ( REDUCTION % MAX )
      Incoming => CO % Incoming % Value
      Max_HWM_Option % Number = Incoming ( 1 )
    end if
    
    if ( present ( Min_HWM_Option ) ) then
      CO % Outgoing % Value = [ HWM % Number ]
      call CO % Reduce ( REDUCTION % MIN )
      Incoming => CO % Incoming % Value
      Min_HWM_Option % Number = Incoming ( 1 )
    end if
    
    if ( present ( Mean_HWM_Option ) .and. present ( C_Option ) ) then
      CO % Outgoing % Value = [ HWM % Number ]
      call CO % Reduce ( REDUCTION % SUM )
      Incoming => CO % Incoming % Value
      Mean_HWM_Option % Number = Incoming ( 1 ) / C_Option % Size
    end if
    
    if ( present ( Max_RSS_Option ) ) then
      CO % Outgoing % Value = [ RSS % Number ]
      call CO % Reduce ( REDUCTION % MAX )
      Incoming => CO % Incoming % Value
      Max_RSS_Option % Number = Incoming ( 1 )
    end if
    
    if ( present ( Min_RSS_Option ) ) then
      CO % Outgoing % Value = [ RSS % Number ]
      call CO % Reduce ( REDUCTION % MIN )
      Incoming => CO % Incoming % Value
      Min_RSS_Option % Number = Incoming ( 1 )
    end if
    
    if ( present ( Mean_RSS_Option ) .and. present ( C_Option ) ) then
      CO % Outgoing % Value = [ RSS % Number ]
      call CO % Reduce ( REDUCTION % SUM )
      Incoming => CO % Incoming % Value
      Mean_RSS_Option % Number = Incoming ( 1 ) / C_Option % Size
    end if

    nullify ( Incoming )

  end subroutine GetMemoryUsage


end module GetMemoryUsage_Command
