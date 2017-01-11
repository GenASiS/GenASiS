!-- TableStreamForm is a class to read file in (ASCII) table form.

module TableStream_Form
  
  use VariableManagement
  use Display
  use DelayFileAccess_Command

  implicit none
  private

    integer ( KDI ), parameter, private :: &
      HANDLE_UNINITIALIZED = - huge ( 1_KDI )

  type, public :: TableStreamForm
    integer ( KDI ) :: &
      Handle      = HANDLE_UNINITIALIZED, &
      ProcessRank =  0, &
      nRows, &
      nColumns
    character ( LDF ) :: &
      Filename = '', &
      Path     = ''
  contains
    procedure, public, pass :: &
      Initialize
    procedure, private, pass :: &
      PrepareRead
    procedure, public, pass :: &
      Read
    final :: &
      Finalize
  end type TableStreamForm
     
    private :: &
      GetNumberOfRows, &
      GetNumberOfColumns, &
      IsIgnored
  
    character ( 1 ), private, parameter :: &
      DELIMITER_COMMENT_1 = '!', &
      DELIMITER_COMMENT_2 = '#', &
      DELIMITER_UNIT      = '@'
    character ( 7 ), private, parameter :: &
      FORMAT_BUFFER = '(a1023)'  !-- must correspond to 
                                 !   LDB = LEN_DEFAULT % BUFFER

contains


  subroutine Initialize ( TS, Filename, ProcessRank, PathOption )

    class ( TableStreamForm ), intent ( inout ) :: &
      TS
    character ( * ), intent ( in ) :: &
      Filename
    integer ( KDI ), intent ( in ) :: &
      ProcessRank
    character ( * ), intent ( in ), optional :: &
      PathOption

    integer ( KDI ) :: &
      Status
    character ( LDB ) :: &
      ErrorMessage
      
    TS % ProcessRank = ProcessRank
    
    TS % Filename = trim ( Filename )
    
    TS % Path = ''
    if ( present ( PathOption ) ) TS % Path = trim ( PathOption )

    call Show ( 'Opening a table file', CONSOLE % INFO_2 )
    call Show ( TS % Filename, 'Name', CONSOLE % INFO_2 )

    call DelayFileAccess ( ProcessRank )
    open &
      ( newunit = TS % Handle, &
        file = trim ( TS % Path ) // trim ( TS % Filename ), &
        action = 'read', status = 'old', iostat = Status, &
        iomsg = ErrorMessage )

    if ( Status /= 0 ) then
      TS % Handle = HANDLE_UNINITIALIZED
      call Show &
             ( 'An error occurred when opening table file', &
               CONSOLE % ERROR )
      call Show ( 'The parameters file may not exist', CONSOLE % ERROR )
      call Show ( TS % Filename, 'Name', CONSOLE % INFO_1 )
      call Show ( ErrorMessage, 'Error Message', CONSOLE % INFO_1 )
    end if

  end subroutine Initialize
  
  
  subroutine PrepareRead ( TS, nRowsOption, nColumnsOption, oRowOption )
    
    class ( TableStreamForm ), intent ( inout ) :: &
      TS
    integer ( KDI ), intent ( in ), optional :: &
      nRowsOption, &
      nColumnsOption, &
      oRowOption
    
    call Show ( 'Reading table file', CONSOLE % INFO_5 )
    call Show ( &
           trim ( TS % Path ) // trim ( TS % Filename ), 'Name', &
           CONSOLE % INFO_5 )
    
    if ( present ( nRowsOption ) ) then
      TS % nRows = nRowsOption
    else
      call GetNumberOfRows ( TS, oRowOption )
    end if
    
    if ( present ( nColumnsOption ) )then
      TS % nColumns = nColumnsOption
    else
      call GetNumberOfColumns ( TS, oRowOption )
    end if
    
    call Show ( TS % nRows, 'nRows', CONSOLE % INFO_2 )
    call Show ( TS % nColumns, 'nColumns', CONSOLE % INFO_2 )
    
  end subroutine PrepareRead 
  
  
  subroutine Read ( TS, Value, nRowsOption, nColumnsOption, oRowOption, &
                    SuccessOption )
  
    class ( TableStreamForm ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( :, : ), allocatable, intent ( out ) :: &
      Value
    integer ( KDI ), intent ( in ), optional :: &
      nRowsOption, &
      nColumnsOption, &
      oRowOption
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    integer ( KDI ) :: &
      iRow, &
      Status
    character ( LDB ) :: &
      Buffer
    logical ( KDL ) :: &
      Success
      
    if ( TS % Handle == HANDLE_UNINITIALIZED ) then
      if ( present ( SuccessOption ) ) SuccessOption = .false.
      return
    end if
    
    call TS % PrepareRead ( nRowsOption, nColumnsOption, oRowOption )

    allocate ( Value ( TS % nRows, TS % nColumns ) )
    
    Success = .false. 
    
    rewind TS % Handle

    if ( present ( oRowOption ) ) then
      do iRow = 1, oRowOption
      read &
        ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      end do
    end if

    iRow = 1
    do
      read &
        ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      if ( Status < 0 ) exit
      if ( IsIgnored ( Buffer ) ) cycle
      read ( Buffer, * ) Value ( iRow, 1 : TS % nColumns )
      if ( iRow == TS % nRows ) then
        Success = .true.
        exit
      end if
      iRow = iRow + 1
    end do
    
    if ( present ( SuccessOption ) ) SuccessOption = Success
  
  end subroutine Read
  
  
  impure elemental subroutine Finalize ( TS )
  
    type ( TableStreamForm ), intent ( inout ) :: &
      TS

    if ( TS % Handle == HANDLE_UNINITIALIZED ) return
    
    call Show &
           ( 'Closing a table file', CONSOLE % INFO_5 )
    call Show &
           ( trim ( TS % Path ) // trim ( TS % Filename ), 'Name', &
             CONSOLE % INFO_5 )
             
    call DelayFileAccess ( TS % ProcessRank )
    close ( unit = TS % Handle )
    
  end subroutine Finalize 
  
  
  subroutine GetNumberOfRows ( TS, oRowOption )
  
    type ( TableStreamForm ), intent ( inout ) :: &
      TS
    integer ( KDI ), intent ( in ), optional :: &
      oRowOption

    integer ( KDI ) :: &
      iRow, &
      Status
    character ( LDB ) :: &
      Buffer
    
    TS % nRows = 0
    
    rewind TS % Handle

    if ( present ( oRowOption ) ) then
      do iRow = 1, oRowOption
      read &
        ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      end do
    end if

    do
      read ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      if ( Status < 0 ) exit
      if ( IsIgnored ( Buffer ) ) cycle
      TS % nRows = TS % nRows + 1
    end do
    
  end subroutine GetNumberOfRows
  
  
  subroutine GetNumberOfColumns ( TS, oRowOption )
  
    type ( TableStreamForm ), intent ( inout ) :: &
      TS
    integer ( KDI ), intent ( in ), optional :: &
      oRowOption
    
    integer ( KDI ) :: &
      iRow, &
      Blank, &
      Status
    logical ( KDL ) :: &
      ReadSingleLine
    character ( LDB ) :: &
      Buffer
    
    ReadSingleLine = .false.
    TS %nColumns = 0
    
    rewind TS % Handle

    if ( present ( oRowOption ) ) then
      do iRow = 1, oRowOption
      read &
        ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      end do
    end if

    do
      read ( unit = TS % Handle, fmt = FORMAT_BUFFER, iostat = Status ) Buffer
      if ( Status < 0 ) exit
      if ( IsIgnored ( Buffer ) ) cycle
      ReadSingleLine = .true.
      exit
    end do
    
    if ( ReadSingleLine ) then
      do while ( len_trim ( Buffer ) > 0 )
        Buffer = adjustl ( Buffer )
        Blank = index ( Buffer, ' ' )
        if ( Blank > 0 )then
          TS % nColumns = TS % nColumns + 1
          Buffer = Buffer ( Blank + 1 : )
        end if
      end do
    end if
    
  end subroutine GetNumberOfColumns
  
  
  function IsIgnored ( Buffer ) result ( II )
  
    character ( * ), intent ( inout ) :: &
      Buffer
    logical ( KDL ) :: &
      II
        
    Buffer = adjustl ( Buffer )
    
    II = ( Buffer ( 1 : 1 ) == DELIMITER_COMMENT_1 &
           .or. Buffer ( 1 : 1 ) == DELIMITER_COMMENT_2 &
           .or. len_trim ( Buffer ) == 0 )
  
  end function IsIgnored
  

end module TableStream_Form
