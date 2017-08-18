!-- GridImageStreamTemplate is an abstract class to handle related file access
!   (i.e. stream) for GridImage (parallel and serial). This class' members 
!   will be needed and inherited by any implementation of GridImageStream 
!   in the form of its extension.

module GridImageStream_Template

  use iso_c_binding
  use VariableManagement
  use Display
  use MessagePassing
  use FileSystemBasics

  implicit none
  private 
  
  type, public, abstract :: GridImageStreamTemplate
    integer ( KDI ) :: &
      ACCESS_SET_GRID, &
      ACCESS_CREATE, &
      ACCESS_WRITE, &
      ACCESS_READ, &
      ACCESS_UNOPENED, &
      HANDLE_UNDEFINED, &
      IGNORABILITY, &
      nBlocks, &
      lName, &
      lDescription, &
      Number, &
      AccessMode, &
      MeshBlockHandle, &
      MultiMeshHandle
    character ( LDL ) :: &
      MeshBlockFileSuffix, &
      MultiMeshFileSuffix
    character ( LDF ) :: &
      Name = '', &
      Description = '', &
      WorkingDirectory = ''
    character ( LDF ), dimension ( : ), allocatable :: &
      MeshBlocksPrefix
    logical ( KDL ) :: &
      Parallel, &
      WorkingDirectoryCreated
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Open
    procedure, public, pass :: &
      IsReadable
    procedure, public, pass :: &
      IsWritable
    procedure ( CloseTemplate ), public, pass, deferred :: &
      Close
    !-- The following subroutines are conceptually "private" but need to be 
    !   declared public so that they are call-able and overidable by 
    !   the child type (extensions)
    procedure ( OpenForWritingTemplate ), public, pass, deferred :: &
      OpenForWriting
    procedure ( OpenForReadingTemplate ), public, pass, deferred :: &
      OpenForReading
    procedure, public, pass :: &
      CloseTemplate
    procedure, public, pass :: &
      OpenForWritingTemplate
    procedure, public, pass :: &
      OpenForReadingTemplate
    procedure, public, pass :: &
      FileCreate
    procedure, public, pass :: &
      FileOpenWrite
    procedure, public, pass :: &
      FileOpenRead
    procedure, public, pass :: &
      FileClose
    procedure, public, pass :: &
      GetNumberOfBlocks
  end type GridImageStreamTemplate
  
  interface

    integer ( c_int ) function MakeSystemDirectory ( Path, Mode ) &
                                 bind ( c, name = 'mkdir' ) 
      use iso_c_binding
      implicit none
      character ( kind = c_char ) :: &
        Path ( * )
      integer ( c_int ), value :: &
        Mode
    end function MakeSystemDirectory 

  end interface
  
  integer ( KDI ), private, parameter :: &
    ACCESS_MODE_UNOPENED = -1, &
    ACCESS_MODE_SET_GRID =  1, &
    ACCESS_MODE_CREATE   =  2, &
    ACCESS_MODE_WRITE    =  3, &
    ACCESS_MODE_READ     =  4
  
  integer ( KDI ), public, protected :: &
    HANDLE_UNDEFINED     = -1
  
contains

  
  subroutine Initialize &
               ( GIS, Name, CommunicatorOption, DescriptionOption, &
                 WorkingDirectoryOption )
    
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    character ( * ), intent ( in ) :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      DescriptionOption, &
      WorkingDirectoryOption
    
    GIS % IGNORABILITY = CONSOLE % INFO_2

    call Show ( 'Initializing a GridImageStream', GIS % IGNORABILITY )
    call Show ( Name, 'Name', GIS % IGNORABILITY )

    GIS % ACCESS_SET_GRID  = ACCESS_MODE_SET_GRID
    GIS % ACCESS_CREATE    = ACCESS_MODE_CREATE
    GIS % ACCESS_WRITE     = ACCESS_MODE_WRITE
    GIS % ACCESS_READ      = ACCESS_MODE_READ
    GIS % ACCESS_UNOPENED  = ACCESS_MODE_UNOPENED
    GIS % HANDLE_UNDEFINED = HANDLE_UNDEFINED
    
    GIS % lName = len_trim ( Name )

    GIS % lDescription = GIS % lName
    if ( present ( DescriptionOption ) ) &
      GIS % lDescription = len_trim ( DescriptionOption )

    GIS % Number          = -1
    GIS % AccessMode      = ACCESS_MODE_UNOPENED
    GIS % MeshBlockHandle = HANDLE_UNDEFINED
    GIS % MultiMeshHandle = HANDLE_UNDEFINED

    GIS % Name = trim ( Name )
    
    GIS % Description = GIS % Name
    if ( present ( DescriptionOption ) ) &
      GIS % Description = DescriptionOption
      
    GIS % WorkingDirectory = '../Output/'
    if ( present ( WorkingDirectoryOption ) ) &
      GIS % WorkingDirectory = WorkingDirectoryOption
    call Show ( GIS % WorkingDirectory, 'WorkingDirectory', GIS % IGNORABILITY )
      
    if ( present ( CommunicatorOption ) ) then
      GIS % nBlocks = CommunicatorOption % Size
      GIS % Parallel = .true.
      GIS % Communicator => CommunicatorOption
    else
      GIS % nBlocks = 0
      GIS % Parallel = .false.
      GIS % Communicator => null ( )
    end if
    
    GIS % WorkingDirectoryCreated = .false.
    
  end subroutine Initialize
  
  
  subroutine Open ( GIS, AccessMode, SeriesOption, NumberOption, &
                    BlockNumberOption )
               
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( in ) :: &
      AccessMode
    logical ( KDL ), intent ( in ), optional :: &
      SeriesOption
    integer ( KDI ), intent ( in ), optional :: &
      NumberOption, &
      BlockNumberOption
    
    if ( AccessMode == ACCESS_MODE_SET_GRID ) then
      GIS % AccessMode = AccessMode
      call Show ( 'Setting a GridImage grid', GIS % IGNORABILITY )
      call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )
    end if

    if ( AccessMode == ACCESS_MODE_CREATE &
         .or. AccessMode == ACCESS_MODE_WRITE ) &
      call GIS % OpenForWriting ( AccessMode, SeriesOption, NumberOption )
      
    if ( AccessMode == ACCESS_MODE_READ ) &
      call GIS % OpenForReading ( NumberOption, BlockNumberOption )
    
  end subroutine Open
  
  
  function IsReadable ( GIS, CheckMultiMeshOption ) result ( IR )
  
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    logical ( KDL ), intent ( in ), optional :: &
      CheckMultiMeshOption
    logical ( KDL ) :: &
      IR
    
    logical ( KDL ) :: &
      CheckMultiMesh
      
    CheckMultiMesh = .false.
    if ( present ( CheckMultiMeshOption ) ) &
      CheckMultiMesh = CheckMultiMeshOption
    
    IR = ( ( GIS % AccessMode == GIS % ACCESS_CREATE &
             .or. GIS % AccessMode == GIS % ACCESS_WRITE &
             .or. GIS % AccessMode == GIS % ACCESS_READ ) &
           .and. GIS % MeshBlockHandle /= GIS % HANDLE_UNDEFINED )
    
    if ( CheckMultiMesh ) &
      IR = IR .and. ( GIS % MultiMeshHandle /= GIS % HANDLE_UNDEFINED )
         
  end function IsReadable
  
  
  function IsWritable ( GIS, CheckMultiMeshOption ) result ( IW )
  
    class ( GridImageStreamTemplate ), intent ( in ) :: &
      GIS
    logical ( KDL ), intent ( in ), optional :: &
      CheckMultiMeshOption
    logical ( KDL ) :: &
      IW
    
    logical ( KDL ) :: &
      CheckMultiMesh
      
    CheckMultiMesh = .false.
    if ( present ( CheckMultiMeshOption ) ) &
      CheckMultiMesh = CheckMultiMeshOption
    
    IW = ( ( GIS % AccessMode == GIS % ACCESS_CREATE &
             .or. GIS % AccessMode == GIS % ACCESS_WRITE ) &
           .and. GIS % MeshBlockHandle /= GIS % HANDLE_UNDEFINED )
    
    if ( CheckMultiMesh ) &
      IW = IW .and. ( GIS % MultiMeshHandle /= GIS % HANDLE_UNDEFINED )
         
  end function IsWritable
  
  
  subroutine CloseTemplate ( GIS )
    
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    
    integer ( KDI ) :: &
      Error
    character ( LDN + 1 ) :: &
      FileNumberString

    if ( GIS % AccessMode == ACCESS_MODE_SET_GRID ) then
      call Show ( 'Finished setting a GridImage grid', GIS % IGNORABILITY )
      call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )
      GIS % AccessMode = ACCESS_MODE_UNOPENED
      return
    end if

    FileNumberString = ''
    if ( GIS % Number >= 0 ) &
      write ( FileNumberString, fmt = '(a1,i7.7)' ) '_', GIS % Number
    
    call Show ( 'Closing a GridImage file', GIS % IGNORABILITY )
    call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )
    call Show ( FileNumberString, 'Number', GIS % IGNORABILITY )
    
    if ( GIS % MeshBlockHandle == HANDLE_UNDEFINED ) return
    
    call GIS % FileClose ( GIS % MeshBlockHandle, Error )
    if ( Error == 0 ) GIS % MeshBlockHandle = HANDLE_UNDEFINED
      
    if ( GIS % MultiMeshHandle /= HANDLE_UNDEFINED ) then
      call GIS % FileClose ( GIS % MultiMeshHandle, Error )
      if ( Error == 0 ) GIS % MultiMeshHandle = HANDLE_UNDEFINED
    end if
    
    GIS % AccessMode = ACCESS_MODE_UNOPENED
    
    if ( allocated ( GIS % MeshBlocksPrefix ) ) &
      deallocate ( GIS % MeshBlocksPrefix )
      
  end subroutine CloseTemplate
  
  
  subroutine OpenForWritingTemplate &
               ( GIS, AccessMode, SeriesOption, NumberOption )
    
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( in ) :: &
      AccessMode
    logical ( KDL ), intent ( in ), optional :: &
      SeriesOption
    integer ( KDI ), intent ( in ), optional :: &
      NumberOption
      
    integer ( KDI ) :: &
      iB, &  !-- iBlock
      lPathName, &
      Error
    logical ( KDL ) :: &
      Series
    character ( LDN + 1 ) :: &
      FileNumberString, &
      BlockNumberString
    character ( LDF ) :: &
      MeshBlockDirectory, &
      PathSuffix, &
      PathName
    
    if ( AccessMode == ACCESS_MODE_CREATE ) GIS % Number = GIS % Number + 1
    
    GIS % AccessMode = AccessMode
    
    Series = .true.
    if ( present ( SeriesOption ) ) Series = SeriesOption
    if ( present ( NumberOption ) ) GIS % Number = NumberOption
    
    if ( Series ) then 
      write ( FileNumberString, fmt = '(a1,i7.7)' ) '_', GIS % Number
    else
      FileNumberString = ''
    end if
    
    if ( GIS % Parallel ) then
      write ( BlockNumberString, fmt = '(a1,i7.7)' ) &
        '_', GIS % Communicator % Rank
    else
      BlockNumberString = ''
    end if
    
    if ( GIS % AccessMode == ACCESS_MODE_CREATE ) then
      call Show ( 'Creating a GridImage file for writing', &
                  GIS % IGNORABILITY )
    else
      call Show ( 'Opening a GridImage file for writing', &
                  GIS % IGNORABILITY )
    end if
    call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )
    call Show ( FileNumberString, 'Number', GIS % IGNORABILITY )

    if ( .not. GIS % WorkingDirectoryCreated ) then
      call Show ( 'Creating working directory if it does not exist', &
                  GIS % IGNORABILITY + 2 )
      call Show ( GIS % WorkingDirectory,  'Working Directory', &
                  GIS % IGNORABILITY + 2 )
      Error = MakeSystemDirectory &
                ( trim ( GIS % WorkingDirectory ) // c_null_char, &
                  int( O'777', c_int ) )
      GIS % WorkingDirectoryCreated = .true.
    end if  
    
    if ( GIS % Parallel ) then
    
      MeshBlockDirectory &
        = trim ( adjustl ( GIS % Name ) ) // trim ( FileNumberString ) &
          // '_MeshBlocks/'
      
      if ( GIS % Communicator % Rank == 0 ) then
        
        !-- Create multi-block mesh file and thereby obtain 
        !   GIS % MultiMeshHandle

        PathSuffix &
          = '_MultiMesh' // trim ( FileNumberString ) &
            // trim ( GIS % MultiMeshFileSuffix )
        PathName &
          = trim ( GIS % WorkingDirectory ) // trim ( adjustl ( GIS % Name ) ) & 
            // trim ( PathSuffix )
        lPathName = len_trim ( PathName )
        
        if ( GIS % AccessMode == ACCESS_MODE_CREATE ) then
          call GIS % FileCreate ( GIS % MultiMeshHandle, PathName, Error )
          Error = MakeSystemDirectory &
                    ( trim ( GIS % WorkingDirectory ) &
                      // trim ( MeshBlockDirectory ) // c_null_char, &
                      int( O'777', c_int ) )
        else 
          call GIS % FileOpenWrite ( GIS % MultiMeshHandle, PathName, Error )
        end if
        
      end if
      
      call GIS % Communicator % Synchronize ( )
      
      !-- Create file for a mesh block and thereby obtain 
      !   GIS % MeshBlockHandle
      
      PathSuffix &
        = trim ( FileNumberString ) // '_MeshBlock' &
          // trim ( BlockNumberString ) // trim ( GIS % MeshBlockFileSuffix )
      PathName &
        = trim ( GIS % WorkingDirectory ) // trim ( MeshBlockDirectory ) &
          // trim ( adjustl ( GIS % Name ) ) // trim ( PathSuffix )
      lPathName = len_trim ( PathName )
       
      if ( GIS % AccessMode == ACCESS_MODE_CREATE ) then
        call GIS % FileCreate ( GIS % MeshBlockHandle, PathName, Error )
      else
        call GIS % FileOpenWrite ( GIS % MeshBlockHandle, PathName, Error )
      end if
        
      if ( GIS % Communicator % Rank == 0 ) then
        allocate ( GIS % MeshBlocksPrefix ( GIS % nBlocks ) )
        do iB = 1, GIS % nBlocks
          write ( BlockNumberString, fmt='(a1,i7.7)' ) '_', iB - 1
          PathSuffix &
            = trim ( FileNumberString ) // '_MeshBlock' &
              // trim ( BlockNumberString ) &
              // trim ( GIS % MeshBlockFileSuffix )
          GIS % MeshBlocksPrefix ( iB ) &
            = trim ( adjustl ( GIS % Name ) ) // trim ( FileNumberString ) &
              // '_MeshBlocks/' // trim ( adjustl ( GIS % Name ) ) &
              // trim ( PathSuffix ) // ':'
        end do
      end if
    
    else  !== .not. GIS % Parallel
    
      PathSuffix &
        = trim ( FileNumberString ) // trim ( GIS % MeshBlockFileSuffix )
      PathName  &
        = trim ( GIS % WorkingDirectory ) &
            // trim ( adjustl ( GIS % Name ) ) // trim ( PathSuffix )
      lPathName = len_trim ( PathName )
      
      !-- Create file for a mesh block and thereby obtain 
      !   GIS % MeshBlockHandle
      
      if ( GIS % AccessMode == ACCESS_MODE_CREATE ) then
        call GIS % FileCreate ( GIS % MeshBlockHandle, PathName, Error )
      else
        call GIS % FileOpenWrite ( GIS % MeshBlockHandle, PathName, Error )
      end if
    
    end if
    
  end subroutine OpenForWritingTemplate
  
  
  subroutine OpenForReadingTemplate &
               ( GIS, NumberOption, BlockNumberOption )
    
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( in ), optional :: &
      NumberOption, &
      BlockNumberOption
         
    integer ( KDI ) :: &
      lPathName, &
      BlockNumber, &
      Error
    character ( LDF ) :: &
      MeshBlocksPrefix, &
      PathName, &
      PathSuffix
    character ( LDN + 1 ) :: &
      FileNumberString, &
      BlockNumberString
    
    if ( GIS % AccessMode == ACCESS_MODE_CREATE &
         .or. GIS % AccessMode == ACCESS_MODE_WRITE ) &
      return
      
    if ( GIS % AccessMode == ACCESS_MODE_READ ) call GIS % Close ( )
    
    FileNumberString = ''
    if ( present ( NumberOption ) ) then
      GIS % Number = NumberOption
      write ( FileNumberString, fmt = '(a1,i7.7)' ) '_', GIS % Number
    end if
    
    if ( GIS % Parallel ) then
      BlockNumber = GIS % Communicator % Rank
      if ( present ( BlockNumberOption ) ) BlockNumber = BlockNumberOption
      write ( BlockNumberString, fmt = '(a1,i7.7)' ) '_', BlockNumber
      !-- Setting of GIS % nBlocks for Parallel case is deferred until 
      !   MultiMesh handle is acquired
    else
      BlockNumberString = ''
      GIS % nBlocks = 0
    end if
        
    GIS % AccessMode = ACCESS_MODE_READ
    
    call Show ( 'Opening a GridImage file for reading', GIS % IGNORABILITY )
    call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )
    call Show ( FileNumberString, 'Number', GIS % IGNORABILITY )
    
    if ( GIS % Parallel ) then
      
      PathSuffix &
        = trim ( FileNumberString ) // '_MeshBlock' &
          // trim ( BlockNumberString ) // trim ( GIS % MeshBlockFileSuffix )
      MeshBlocksPrefix &
        = trim ( adjustl ( GIS % Name ) ) // trim ( FileNumberString ) &
          // '_MeshBlocks/'
      PathName  &
        = trim ( GIS % WorkingDirectory ) // trim ( MeshBlocksPrefix ) &
          // trim ( adjustl ( GIS % Name ) ) // trim ( PathSuffix )
      lPathName = len_trim ( PathName )
    
    else  !-- .not. GIS % Parallel
      
      MeshBlocksPrefix = ''
      PathSuffix &
        = trim ( FileNumberString ) // trim ( GIS % MeshBlockFileSuffix )
      PathName  &
        = trim ( GIS % WorkingDirectory ) // trim ( MeshBlocksPrefix ) &
          // trim ( adjustl ( GIS % Name ) ) // trim ( PathSuffix )
      lPathName = len_trim ( PathName )
      
    end if
    
    if ( GIS % Parallel ) &
      call DelayFileAccess ( GIS % Communicator % Rank )
    
    call GIS % FileOpenRead ( GIS % MeshBlockHandle, PathName, Error )
    
    !-- Open MultiMesh file
    
    if ( GIS % Parallel ) then
      
      call DelayFileAccess ( GIS % Communicator % Rank )
    
      PathSuffix &
        = '_MultiMesh' // trim ( FileNumberString ) &
          // trim ( GIS % MultiMeshFileSuffix )
      PathName &
        = trim ( GIS % WorkingDirectory ) // trim ( adjustl ( GIS % Name ) ) &
          // trim ( PathSuffix )
      lPathName = len_trim ( PathName )
      
      call GIS % FileOpenRead ( GIS % MultiMeshHandle, PathName, Error )
      
      if ( index ( GIS % Name, '1D' ) > 0 ) then
        GIS % nBlocks = 1
      else
        call GIS % GetNumberOfBlocks ( )
      end if
      
    end if
      
  end subroutine OpenForReadingTemplate
  
  
  subroutine FileCreate ( GIS, Handle, PathName, Status )
  
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status
      
    open ( newunit = Handle, file = trim ( PathName ), &
           action = 'write', status = 'replace', iostat = Status )    
  
  end subroutine FileCreate
  
  
  subroutine FileOpenWrite ( GIS, Handle, PathName, Status )
  
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status
    
    open ( newunit = Handle, file = trim ( PathName ), &
           action = 'write', status = 'old', iostat = Status )    
  
  end subroutine FileOpenWrite
  
  
  subroutine FileOpenRead ( GIS, Handle, PathName, Status )
  
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status
    
    open &
      ( newunit = Handle, file = trim ( PathName ), &
        action = 'read', status = 'old', iostat = Status )    
  
  end subroutine FileOpenRead
  
  
  subroutine FileClose ( GIS, Handle, Status )
  
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    integer ( KDI ), intent ( out ) :: &
      Status
    
    close ( Handle, iostat = Status )
  
  end subroutine FileClose
  
  
  subroutine GetNumberOfBlocks ( GIS )
    
    class ( GridImageStreamTemplate ), intent ( inout ) :: &
      GIS
      
    if ( GIS % MultiMeshHandle == HANDLE_UNDEFINED ) then
      GIS % nBlocks = -1
      return
    end if
    
    if ( associated ( GIS % Communicator ) ) &
      GIS % nBlocks = GIS % Communicator % Size
    
!     call SiloWrappers_DBGetNumberOfBlcks ( GIS % MultiMeshHandle, nBlocks )
  
  end subroutine GetNumberOfBlocks
  
  
end module GridImageStream_Template
