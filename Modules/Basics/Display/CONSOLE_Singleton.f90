!-- ConsoleSingleton extends ConsoleHeaderForm and specifies the process rank
!   that displays to the "console" (i.e. STDOUT) with the accompanying methods
!   to modify the display properties.

module CONSOLE_Singleton

  use VariableManagement
  use ConsoleHeader_Form

  implicit none
  private

  type, extends ( ConsoleHeaderForm ), public :: ConsoleSingleton
    integer ( KDI ) :: &
      Verbosity   = 5, &  !-- INFO_3, see ConsoleHeader_Form.f90
      DisplayRank = 0, &
      ProcessRank = 0
    logical ( KDL ) :: &
      Muted = .false.
    procedure ( AbortInterface ), nopass, pointer :: &
      Abort => null ( )
  contains
    procedure, public, nopass :: &
      Initialize
    procedure, public, nopass :: &
      SetVerbosity
    procedure, public, nopass :: &
      SetDisplayRank
    procedure, public, nopass :: &
      Mute
    procedure, public, nopass :: &
      Unmute
    procedure, private, nopass :: &
      Abort_C
    final :: &
      Finalize
  end type ConsoleSingleton
  
  type ( ConsoleSingleton ), public, protected :: &
    CONSOLE

  interface

    subroutine ShowInteger &
                 ( Integer, Description, IgnorabilityOption, &
                   DisplayRankOption, nLeadingLinesOption, &
                   nTrailingLinesOption )
      use VariableManagement
      integer ( KDI ), intent ( in ) :: &
        Integer
      character ( * ), intent ( in ) :: &
        Description
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption, &
        DisplayRankOption, &
        nLeadingLinesOption, &
        nTrailingLinesOption
    end subroutine ShowInteger

    subroutine ShowCharacter &
                 ( Character, Description, IgnorabilityOption, &
                   DisplayRankOption, nLeadingLinesOption, &
                   nTrailingLinesOption )
      use VariableManagement
      character ( * ), intent ( in ) :: &
        Character, &
        Description
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption, &
        DisplayRankOption, &
        nLeadingLinesOption, &
        nTrailingLinesOption
    end subroutine ShowCharacter

    subroutine ShowMessage &
                 ( Character, IgnorabilityOption, &
                   DisplayRankOption, nLeadingLinesOption, &
                   nTrailingLinesOption )
      use VariableManagement
      character ( * ), intent ( in ) :: &
        Character
      integer ( KDI ), intent ( in ), optional :: &
        IgnorabilityOption, &
        DisplayRankOption, &
        nLeadingLinesOption, &
        nTrailingLinesOption
    end subroutine ShowMessage
    
    subroutine AbortInterface ( )
      use Specifiers
    end subroutine AbortInterface
    
  end interface
  

contains


  subroutine Initialize ( ProcessRank, AbortOption )

    integer ( KDI ), intent ( in ) :: &
      ProcessRank
    procedure ( AbortInterface ), intent ( in ), pointer, optional :: &
      AbortOption

    CONSOLE % Verbosity = CONSOLE % INFO_3
    CONSOLE % ProcessRank = ProcessRank
    !CONSOLE % Abort => Abort_C
    if ( present ( AbortOption ) ) CONSOLE % Abort => AbortOption

  end subroutine Initialize

 
  subroutine SetVerbosity ( Verbosity )

    character ( * ), intent ( in ) :: &
      Verbosity

    call ShowMessage ( 'Modifying CONSOLE', CONSOLE % INFO_1 )

    select case ( trim ( Verbosity ) )
    case ( 'ERROR' ) 
      CONSOLE % Verbosity = CONSOLE % ERROR
    case ( 'WARNING' ) 
      CONSOLE % Verbosity = CONSOLE % WARNING
    case ( 'INFO_1' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_1
    case ( 'INFO_2' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_2
    case ( 'INFO_3' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_3
    case ( 'INFO_4' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_4
    case ( 'INFO_5' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_5
    case ( 'INFO_6' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_6
    case ( 'INFO_7' ) 
      CONSOLE % Verbosity = CONSOLE % INFO_7
    case default
      call ShowMessage &
             ( 'Unknown display verbosity. Reverting to default.', &
               CONSOLE % WARNING )
    end select

    call ShowCharacter &
           ( CONSOLE % LABEL ( CONSOLE % Verbosity ), 'Verbosity', &
             CONSOLE % INFO_1 )

  end subroutine SetVerbosity
  
  
  subroutine SetDisplayRank ( DisplayRank )

    integer ( KDI ), intent ( in ) :: &
      DisplayRank

    integer ( KDI ) :: &
      DisplayRankOld

    DisplayRankOld = CONSOLE % DisplayRank

    CONSOLE % DisplayRank = DisplayRank

    call ShowMessage ( 'Modifying CONSOLE', CONSOLE % INFO_1 )
    call ShowInteger ( CONSOLE % DisplayRank, 'DisplayRank', CONSOLE % INFO_1 )

  end subroutine SetDisplayRank
  
  
  subroutine Mute ( )

    call ShowMessage ( 'Muting CONSOLE', CONSOLE % INFO_1 )

    CONSOLE % Muted = .true.

  end subroutine Mute

  
  subroutine Unmute ( )

    CONSOLE % Muted = .false.
    call ShowMessage ( 'Unmuting CONSOLE', CONSOLE % INFO_1 )

  end subroutine Unmute
  

  elemental subroutine Finalize ( C )

    type ( ConsoleSingleton ), intent ( inout ) :: &
      C

    nullify ( C % Abort )

  end subroutine Finalize

  
  subroutine Abort_C ( )
    
    stop 1
  
  end subroutine Abort_C


end module CONSOLE_Singleton
