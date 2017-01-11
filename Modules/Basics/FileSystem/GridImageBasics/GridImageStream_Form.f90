!-- GridImageStreamSilo is a class to handle access Silo file format and its 
!   structures.

module GridImageStream_Form

  use iso_c_binding
  use VariableManagement
  use Display
  use MessagePassing
  use FileSystemBasics
  use GridImageStream_Template
  use DB_TableOfContentsType_Silo_C

  implicit none
  private 
  
  include 'silo_f9x.inc'
  
  type, public, extends ( GridImageStreamTemplate ) :: GridImageStreamForm
    character ( LDL ) :: &
      ContentType
    character ( LDL ), dimension ( : ), allocatable :: &
      ContentList
    character ( LDF ) :: &
      CurrentDirectory
    character ( LDF ), dimension ( : ), allocatable, private :: &
      CreatedDirectory
  contains
    procedure, public, pass :: &
      MakeDirectory
    procedure, public, pass :: &
      ChangeDirectory
    procedure, public, pass :: &
      ListContents
    procedure, public, pass :: &
      AccessSiloPointer
    procedure, public, pass :: &
      WriteVariableGroup
    procedure, public, pass :: &
      ReadVariableGroup
    procedure, public, pass :: &
      Close
    procedure, public, pass :: &
      OpenForWriting
    procedure, public, pass :: &
      OpenForReading
    procedure, public, pass :: &
      FileCreate
    procedure, public, pass :: &
      FileOpenWrite
    procedure, public, pass :: &
      FileOpenRead
    procedure, public, pass :: &
      FileClose
    final :: &
      Finalize
  end type GridImageStreamForm
  
  interface

    type ( c_ptr ) function DB_FortranAccessPointer ( Handle ) &
                              bind ( c, name = 'DBFortranAccessPointer' )
      use iso_c_binding
      implicit none
      integer ( c_int ), value :: &
        Handle
    end function DB_FortranAccessPointer
    
    type ( c_ptr ) function DB_GetVariable ( Handle, VariableName ) &
                              bind ( c, name = 'DBGetVar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Handle
      character ( kind = c_char ) :: &
        VariableName ( * )
    end function DB_GetVariable
    
    subroutine DB_FreeVariable ( Handle ) bind ( c, name = 'free' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        Handle
    end subroutine DB_FreeVariable
    
  end interface
    
  private :: &
    DB_FortranAccessPointer
    
  private :: &
    TrimDoubleSlash
  
  integer ( KDI ), private, parameter :: &
    MAX_DIRECTORIES      = 512
  
  character ( 2 ), private, parameter :: & 
  !-- This must match with the definition in SiloWrappers.c
    SD = '@@', &   ! Section Delimiter
    ND = ';;', &   ! Name Delimiter
    VD = '##'      ! Value Delimiter

contains

  
  recursive subroutine MakeDirectory ( GIS, Name )
    
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    character ( * ), intent ( in ) :: &
      Name
    
    integer ( KDI ) :: &
      iP, &
      iD, &
      Error
    logical ( KDL ) :: &
      DirectoryExist
    character ( LDF ) :: &
      AbsolutePathName
    character ( LDF ), dimension ( : ), allocatable :: &
      Parent
    
    if ( .not. GIS % IsWritable ( ) ) return 
    
    if ( len_trim ( Name ) == 0 ) return
    if ( trim ( Name ) == '/' ) return
    if ( trim ( Name ) == '.' .or. trim ( Name ) == './' ) return
    if ( trim ( Name ) == '..' .or. trim ( Name ) == '../' ) return
    
    if ( index ( Name, '/' ) > 0 ) then
      call Split ( Name, '/', Parent )
      do iP = 1, size ( Parent )
        call GIS % MakeDirectory ( Parent ( iP ) )
      end do
      return  
    end if
    
    AbsolutePathName = trim ( GIS % CurrentDirectory ) // trim ( Name ) // '/'
    call TrimDoubleSlash ( AbsolutePathName ) 
    
    DirectoryExist = .false.
    do iD = 1, MAX_DIRECTORIES
      if ( GIS % CreatedDirectory ( iD ) == AbsolutePathName ) then
        DirectoryExist = .true.
        exit
      end if
      if ( GIS % CreatedDirectory ( iD ) == '' ) exit
    end do
    
    if ( .not. DirectoryExist ) &
      Error = DBMKDIR &
                ( GIS % MeshBlockHandle, trim ( Name ), len_trim ( Name ), &
                  Error )
    Error = DBSETDIR &
              ( GIS % MeshBlockHandle, trim ( Name ), len_trim ( Name ) )
    
    GIS % CurrentDirectory = AbsolutePathName
    GIS % CreatedDirectory ( iD ) = AbsolutePathName
    
    if ( .not. GIS % Parallel ) return
    
    if ( GIS % Communicator % Rank == 0 ) then
      if ( GIS % MultiMeshHandle == HANDLE_UNDEFINED ) return
      if ( .not. DirectoryExist ) &
        Error = DBMKDIR &
                  ( GIS % MultiMeshHandle, trim ( Name ), len_trim ( Name ), &
                    Error )
      Error = DBSETDIR &
                ( GIS % MultiMeshHandle, trim ( Name ), len_trim ( Name ) )
    end if
    
    call Show &
           ( trim ( GIS % CurrentDirectory ), 'CurrentDirectory', &
             CONSOLE % INFO_7 )
    
  end subroutine MakeDirectory
  
  
  subroutine ChangeDirectory ( GIS, PathName, SuccessOption )
    
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    character ( * ), intent ( in ) :: &
      PathName
    logical ( KDL ), intent ( out ), optional :: &
      SuccessOption
    
    integer ( KDI ) :: &
      Error, &
      LastSlash
    character ( LDF ) :: &
      Buffer
    
    if ( len_trim ( PathName ) == 0 ) return
    
    !-- Early return means change directory fails
    if ( present ( SuccessOption ) ) SuccessOption = .false.
    
    if ( .not. GIS % IsReadable ( ) ) return 
    
    if ( trim ( Pathname ) == '.' .or. trim ( Pathname ) == './' )  return

    Error = DBSETDIR &
              ( GIS % MeshBlockHandle, trim ( PathName ), &
                len_trim ( PathName ) )
    if ( Error /= 0 ) return
    
    !-- Set GIS % CurrentDirectory to reflect the new directory we are in
    
    if ( Pathname ( : 1 ) == '/' ) then
      GIS % CurrentDirectory = trim ( Pathname ) // '/'
      call TrimDoubleSlash ( GIS % CurrentDirectory )
    else if ( PathName ( : 2 ) == '..' ) then
      Buffer = GIS % CurrentDirectory
      LastSlash &
        = index ( Buffer ( : len_trim ( Buffer ) - 1 ), '/', back = .true. )
      GIS % CurrentDirectory = GIS % CurrentDirectory ( : LastSlash )
      if ( len_trim ( PathName ) > 3 ) &
        GIS % CurrentDirectory &
          = trim ( GIS % CurrentDirectory ) // trim ( PathName ( 4: ) ) // '/'
      call TrimDoubleSlash ( GIS % CurrentDirectory )
    else
      GIS % CurrentDirectory &
        = trim ( GIS % CurrentDirectory ) // trim ( PathName ) // '/'
      call TrimDoubleSlash ( GIS % CurrentDirectory )
    end if
    
    if ( GIS % Parallel ) then
      if ( GIS % Communicator % Rank == 0 ) then
        if ( GIS % MultiMeshHandle == HANDLE_UNDEFINED ) return
        Error = DBSETDIR &
                  ( GIS % MultiMeshHandle, trim ( PathName ), &
                    len_trim ( PathName ) )
        if ( Error /= 0 ) return
      end if
    end if
    
    call Show &
           ( trim ( GIS % CurrentDirectory ), 'CurrentDirectory', &
             CONSOLE % INFO_7 )
      
    if ( present ( SuccessOption ) ) SuccessOption = .true.
    
  end subroutine ChangeDirectory
  

  subroutine ListContents ( GIS, ContentTypeOption )

    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    character ( * ), intent ( in ), optional :: &
      ContentTypeOption
   
    integer ( KDI ) :: &
      iC, &
      iCh, &
      nContents
    character ( LDL ) :: &
      Name, &
      ContentType
    character ( LDB ), pointer :: &
      Buffer
    type ( c_ptr ) :: &
      DB_File = c_null_ptr, &
      DB_TOC = c_null_ptr
    type ( c_ptr ), dimension ( : ), pointer :: &
      StringArray 
    type ( DB_TableOfContentsType ), pointer :: &
      Contents
     
    ContentType = 'Directory'
    if ( present ( ContentTypeOption ) ) then
      ContentType = ContentTypeOption
      GIS % ContentType = ContentType
    end if
   
    DB_File = AccessSiloPointer ( GIS, GIS % MeshBlockHandle )
    DB_TOC = DB_Get_TOC ( DB_File )
    call c_f_pointer ( DB_TOC, Contents )

    nContents = 0
    if ( trim ( ContentType ) == 'Directory' ) then
      call c_f_pointer &
        ( Contents % Directory, StringArray, [ Contents % nDirectories ] )
      nContents = Contents % nDirectories
    else if ( trim ( ContentType ) == 'Curve' ) then
      call c_f_pointer &
        ( Contents % Curve, StringArray, [ Contents % nCurves ] )
      nContents = Contents % nCurves
    else if ( trim ( ContentType ) == 'PointGridVariable' ) then
      call c_f_pointer &
        ( Contents % PointVariable, StringArray, &
          [ Contents % nPointVariables ] )
      nContents = Contents % nPointVariables
    else if ( trim ( ContentType ) == 'StructuredGridVariable' ) then
      call c_f_pointer &
        ( Contents % QuadVariable, StringArray, [ Contents % nQuadVariables ] )
      nContents = Contents % nQuadVariables
    else if ( trim ( ContentType ) == 'UnstructuredGridVariable' ) then
      call c_f_pointer &
        ( Contents % UnstructuredVariable, StringArray, &
          [ Contents % nUnstructuredVariables ] )
      nContents = Contents % nUnstructuredVariables
    else if ( trim ( ContentType ) == 'Variable' ) then
      call c_f_pointer &
        ( Contents % Variable, StringArray, &
          [ Contents % nVariables ] )
      nContents = Contents % nVariables
    end if
   
    if ( allocated ( GIS % ContentList ) ) deallocate ( GIS % ContentList )  
    allocate ( GIS % ContentList ( nContents ) )
       
    do iC = 1, nContents
      call c_f_pointer ( StringArray ( iC ), Buffer )
      iCh = 1
      do while ( Buffer ( iCh : iCh ) /= c_null_char )
        iCh = iCh + 1
      end do
      GIS % ContentList ( iC ) = Buffer ( 1 : iCh - 1 )
    end do

    call Show &
           ( GIS % ContentList, 'List of ' // trim ( ContentType ), &
             GIS % IGNORABILITY + 2 )
  
  end subroutine ListContents


  function AccessSiloPointer ( GIS, Handle ) result ( SP )
  
    class ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    integer ( KDI ), intent ( in ) :: &
      Handle
    
    type ( c_ptr ) :: &
      SP
      
    SP = DB_FortranAccessPointer ( int ( Handle,  c_int ) )
  
  end function AccessSiloPointer
  
  
  subroutine WriteVariableGroup ( GIS, VG, DirectoryOption )
  
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    class ( VariableGroupForm ), dimension ( : ), intent ( in ) :: &
      VG
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    
    integer ( KDI )  :: &
      iG, &
      iS, &
      iVrbl, &
      Error
    logical ( KDL ) :: &
      Success
    character ( LDF ) :: &
      WorkingDirectory
      
    WorkingDirectory = GIS % CurrentDirectory
    
    if ( present ( DirectoryOption ) ) then
      call GIS % ChangeDirectory ( DirectoryOption, SuccessOption = Success )
      if ( .not. Success ) &
        call GIS % MakeDirectory ( DirectoryOption )
    end if
    
    do iG = 1, size ( VG )
      call GIS % MakeDirectory ( VG ( iG ) % Name )
      do iS = 1, VG ( iG ) % nVariables
        iVrbl = VG ( iG ) % iaSelected ( iS )
        Error = DBWRITE &
                  ( GIS % MeshBlockHandle, &
                    VG ( iG ) % Variable ( iVrbl ), & 
                    VG ( iG ) % lVariable ( iVrbl ), &
                    VG ( iG ) % Value ( :, iVrbl ), &
                    [ VG ( iG ) % nValues ], 1, DB_DOUBLE )
      end do
      call GIS % ChangeDirectory ( '../' )
    end do
    
    call GIS % ChangeDirectory ( WorkingDirectory )
  
  end subroutine WriteVariableGroup
  
  
  subroutine ReadVariableGroup ( GIS, VG, DirectoryOption )
  
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    class ( VariableGroupForm ), dimension ( : ), intent ( out ) :: &
      VG
    character ( * ), intent ( in ), optional :: &
      DirectoryOption
    
    integer ( KDI )  :: &
      iG, &
      iD, &   !-- iDirectory
      iV, &
      iVrbl, &
      nValues, &
      Error
    logical ( KDL ) :: &
      Success
    character ( LDL ), dimension ( : ), allocatable :: &
      Directory
    character ( LDF ) :: &
      WorkingDirectory
    type ( c_ptr ) :: &
      Value, &
      DB_File
    real ( c_double ), dimension ( : ), pointer :: &
      VariableValue
    
    WorkingDirectory = GIS % CurrentDirectory
    
    if ( present ( DirectoryOption ) ) &
      call GIS % ChangeDirectory ( DirectoryOption )
    
    call GIS % ListContents ( )
    allocate ( Directory ( size ( GIS % ContentList ) ) )
    Directory = GIS % ContentList
    
    DB_File = GIS % AccessSiloPointer ( GIS % MeshBlockHandle )
    
    iG = 0
    do iD = 1, size ( Directory )
      
      call GIS % ChangeDirectory ( Directory ( iD ) )
      call GIS % ListContents ( 'Variable' )
      if ( size ( GIS % ContentList ) == 0 ) then
        call GIS % ChangeDirectory ( '../' )
        cycle
      end if
      
      iG = iG + 1
      !-- assume all variables in this directory has the same number of 
      !   elements since they were created from a VariableGroup object
      Error = DBINQLEN &
                  ( GIS % MeshBlockHandle, GIS % ContentList ( 1 ), &
                    len_trim ( GIS % ContentList ( 1 ) ), nValues )
      
      call VG ( iG ) % Initialize &
             ( [ nValues, size ( GIS % ContentList ) ], &
               NameOption = trim ( Directory ( iD ) ), &
               VariableOption = GIS % ContentList )
      
      do iVrbl = 1, VG ( iG ) % nVariables 
        Value = DB_GetVariable &
                  ( DB_File, &
                    trim ( VG ( iG ) % Variable ( iVrbl ) ) // c_null_char )
        call c_f_pointer ( Value, VariableValue, [ nValues ] )
        do iV = 1, nValues
          VG ( iG ) % Value ( iV, iVrbl ) = VariableValue ( iV ) 
        end do
        call DB_FreeVariable ( Value )
      end do 
      
      call GIS % ChangeDirectory ( '../' )
      
    end do
    
    call GIS % ChangeDirectory ( WorkingDirectory )
    
  end subroutine ReadVariableGroup
  
  
  subroutine Close ( GIS )
    
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
      
    call GIS % CloseTemplate ( )
    
    if ( allocated ( GIS % CreatedDirectory ) ) &
      deallocate ( GIS % CreatedDirectory )
    
  end subroutine Close
  
  
  impure elemental subroutine Finalize ( GIS )
  
    type ( GridImageStreamForm ), intent ( inout ) :: &
      GIS

    nullify ( GIS % Communicator ) 

    if ( allocated ( GIS % CreatedDirectory ) ) &
      deallocate ( GIS % CreatedDirectory )
    if ( allocated ( GIS % MeshBlocksPrefix ) ) &
      deallocate ( GIS % MeshBlocksPrefix )
    if ( allocated ( GIS % ContentList ) ) &
      deallocate ( GIS % ContentList )

    if ( GIS % Name == '' ) return

    call Show ( 'Finalizing a GridImageStream', GIS % IGNORABILITY )
    call Show ( GIS % Name, 'Name', GIS % IGNORABILITY )

  end subroutine Finalize
  
  
  subroutine OpenForWriting ( GIS, AccessMode, SeriesOption, NumberOption )
    
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( in ) :: &
      AccessMode
    logical ( KDL ), intent ( in ), optional :: &
      SeriesOption
    integer ( KDI ), intent ( in ), optional :: &
      NumberOption
      
    GIS % MeshBlockFileSuffix = '.silo'
    GIS % MultiMeshFileSuffix = '.silo'
      
    call GIS % OpenForWritingTemplate &
           ( AccessMode, SeriesOption, NumberOption )
    
    GIS % CurrentDirectory = '/'
    
    allocate ( GIS % CreatedDirectory ( MAX_DIRECTORIES ) )
    GIS % CreatedDirectory = ''
    
  end subroutine OpenForWriting
  
  
  subroutine OpenForReading ( GIS, NumberOption, BlockNumberOption )
    
    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( in ), optional :: &
      NumberOption, &
      BlockNumberOption
         
    GIS % MeshBlockFileSuffix = '.silo'
    GIS % MultiMeshFileSuffix = '.silo'
    
    call GIS % OpenForReadingTemplate ( NumberOption, BlockNumberOption )
    
    GIS % CurrentDirectory = '/'
    
  end subroutine OpenForReading


  subroutine FileCreate ( GIS, Handle, PathName, Status )

    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status
    
    Status = DBCREATE &
              ( PathName, len_trim ( PathName ), DB_CLOBBER, DB_LOCAL, &
                GIS % Description, GIS % lDescription, DB_PDB, Handle )

  end subroutine FileCreate


  subroutine FileOpenWrite ( GIS, Handle, PathName, Status )

    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status

    Status = DBOPEN &
              ( PathName, len_trim ( PathName ), DB_PDB, DB_APPEND, Handle )

  end subroutine FileOpenWrite


  subroutine FileOpenRead ( GIS, Handle, PathName, Status )

    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS
    integer ( KDI ), intent ( inout ) :: &
      Handle
    character ( * ), intent ( in ) :: &
      PathName
    integer ( KDI ), intent ( out ) :: &
      Status

    Status = DBOPEN &
              ( PathName, len_trim ( PathName ), DB_PDB, DB_READ, Handle )
  
  end subroutine FileOpenRead

  
  subroutine FileClose ( GIS, Handle, Status )   

    class ( GridImageStreamForm ), intent ( inout ) :: &
      GIS 
    integer ( KDI ), intent ( inout ) :: &
      Handle
    integer ( KDI ), intent ( out ) :: &
      Status

    Status = DBCLOSE ( Handle )

  end subroutine FileClose

  
  subroutine TrimDoubleSlash ( Path )
  
    character ( * ), intent ( inout ) :: &
      Path
      
    integer ( KDI ) :: &
      TrimmedLength
      
    TrimmedLength = len_trim ( Path )
    
    if ( TrimmedLength <= 1 ) return
    
    if ( Path ( TrimmedLength-1:TrimmedLength ) == '//' ) then
      Path = Path ( : TrimmedLength - 1 )
    end if
  
  end subroutine TrimDoubleSlash
  
    
end module GridImageStream_Form
