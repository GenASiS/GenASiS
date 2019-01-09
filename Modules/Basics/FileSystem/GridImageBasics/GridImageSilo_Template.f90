!-- GridImageSiloTemplate is an abstract class whose members will be needed 
!   and inherited by the different types of Silo GridImage in the form of its
!   extension.

module GridImageSilo_Template

  use iso_c_binding
  use Specifiers
  use DataManagement
  use Display
  use GridImageStream_Form
  use DB_QuadMeshType_Silo_C
  use GridImage_Template
        
  implicit none
  private
  
  include 'silo_f9x.inc'

  type, public, abstract, extends ( GridImageTemplate ) :: &
    GridImageSiloTemplate
    integer ( KDI ) :: &
      lDirectory = 0, &
      MultiMeshType, &
      MultiVariableType
    character ( LDF ) :: &
      Directory = ''
    type ( GridImageStreamForm ), pointer :: &
      Stream => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      WriteHeader
    procedure, public, pass :: &
      WriteMultiMesh
    procedure, public, pass :: &
      WriteMultiVariable
    procedure, public, pass :: &
      WriteVectorVariable
    procedure ( WriteInterface ), public, pass, deferred :: &
      Write
    procedure ( SetReadAttributesInterface ), public, pass, deferred :: &
      SetReadAttributes
    procedure, public, pass :: &
      ReadHeader
    procedure ( ReadInterface ), public, pass, deferred :: &
      Read
  end type GridImageSiloTemplate
  
  
  abstract interface 
  
    subroutine SetReadAttributesInterface ( GI, Directory, oValue )
      use Specifiers
      import GridImageSiloTemplate
      class ( GridImageSiloTemplate ), intent ( inout ) :: &
        GI  
      character ( * ), intent ( in ) :: &
        Directory
      integer ( KDI ), intent ( in ) :: &
        oValue
    end subroutine
    
    subroutine WriteInterface ( GI, TimeOption, CycleNumberOption )
      use Specifiers
      import GridImageSiloTemplate
      class ( GridImageSiloTemplate ), intent ( inout ) :: &
        GI
      type ( MeasuredValueForm ), intent ( in ), optional :: &
        TimeOption
      integer ( KDI ), intent ( in ), optional :: &
        CycleNumberOption
    end subroutine WriteInterface
    
    subroutine ReadInterface ( GI, TimeOption, CycleNumberOption )
      use Specifiers
      import GridImageSiloTemplate
      class ( GridImageSiloTemplate ), intent ( inout ) :: &
        GI
      type ( MeasuredValueForm ), intent ( out ), optional :: &
        TimeOption
      integer ( KDI ), intent ( out ), optional :: &
        CycleNumberOption
    end subroutine ReadInterface
  
  end interface 
  
contains

  subroutine Initialize ( GI, S )
  
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    class ( GridImageStreamForm ), intent ( in ), target :: &
      S
    
    call GI % InitializeTemplate ( )
    GI % Stream => S
  
  end subroutine Initialize

  
  subroutine WriteHeader ( GI, TimeOption, CycleNumberOption )
  
    !-- WriteHeader writes a fake mesh with time and cycle number information
    !   at the root directory of Silo file since VisIt cannot read time 
    !   and cycle number of mesh that does not reside in the root directory.
  
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      nSiloOptions, &
      lHeaderName, &
      SiloOptionList, &
      Error
    real ( KDR ), dimension ( 1 ) :: &
      Dummy_X, &
      Dummy_Y
    real ( KDR ), dimension ( 1, 1 ) :: &
      DummyValue
    character ( LDF ) :: &
      WorkingDirectory, &
      HeaderName
    
    if ( .not. GI % Stream % IsWritable ( ) ) return
    
    nSiloOptions = 1
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
    
    Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
    Error = DBADDIOPT ( SiloOptionList, DBOPT_HIDE_FROM_GUI, 1 )
    if ( present ( CycleNumberOption ) ) &
      Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
    if ( present ( TimeOption ) ) &
      Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
    
    Dummy_X    = 0.0_KDR
    Dummy_Y    = 0.0_KDR
    DummyValue = 0.0_KDR

    WorkingDirectory = GI % Stream % CurrentDirectory
    call GI % Stream % ChangeDirectory ( '/' )

    HeaderName = 'Header'
    lHeaderName = len_trim ( HeaderName )

!-- FIXME: Tweak for NAG 5.3.1    
!    Error &
!      = DBPUTQM &
!          ( GI % Stream % MeshBlockHandle, HeaderName, lHeaderName, &
!            'x', 1, 'y', 1, 'z', 1, Dummy_X, Dummy_Y, DB_F77NULL, &
!            [ 1, 1 ], 2, DB_DOUBLE, DB_COLLINEAR, SiloOptionList, Error )
    Error &
      = DBPUTQM &
          ( GI % Stream % MeshBlockHandle, HeaderName, lHeaderName, &
            'x', 1, 'y', 1, 'z', 1, Dummy_X, Dummy_Y, Dummy_Y, &
            [ 1, 1 ], 2, DB_DOUBLE, DB_COLLINEAR, SiloOptionList, Error )
    
    Error = DBFREEOPTLIST ( SiloOptionList )
    
    call WriteMultiMesh &
           ( GI, HeaderName, TimeOption, CycleNumberOption, &
             HideOption = .true. )
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )
    
  end subroutine WriteHeader
  
  
  subroutine WriteMultiMesh &
               ( GI, Name, TimeOption, CycleNumberOption, HideOption )
    
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    !-- Break convention on argument ordering for consistency with similar
    !   subroutines in StructuredGridImageSiloForm and 
    !   UnstructuredGridImageSiloForm
    logical ( KDL ), intent ( in ), optional :: &
      HideOption

    integer ( KDI ) :: &
      iB, &
      nBlocks, &
      nSiloOptions, &
      SiloOptionList, &
      Error
    integer ( KDI ), dimension ( : ), allocatable :: &
      lMultiMeshName, &
      MultiMeshType
    character ( LDF ), dimension ( : ), allocatable :: &
      MultiMeshName
      
    if ( .not. GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) ) &
      return
      
    nBlocks = GI % Stream % nBlocks
    allocate ( MultiMeshName ( nBlocks ) )
    allocate ( lMultiMeshName ( nBlocks ) )
    allocate ( MultiMeshType ( nBlocks ) )
    
    nSiloOptions = 0 
    SiloOptionList = DB_F77NULL
    
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( HideOption ) ) &
      nSiloOptions = nSiloOptions + 1
      
    if ( nSiloOptions > 0 ) then
      Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
      if ( present ( CycleNumberOption ) ) &
        Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
      if ( present ( TimeOption ) ) &
        Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
      if ( present ( HideOption ) ) then
        if ( HideOption ) then
          Error = DBADDIOPT ( SiloOptionList, DBOPT_HIDE_FROM_GUI, 1 )
        else
          Error = DBADDIOPT ( SiloOptionList, DBOPT_HIDE_FROM_GUI, 0 )
        end if
      end if
    end if
    
    do iB = 1, nBlocks
      MultiMeshName ( iB ) & 
        = trim ( GI % Stream % MeshBlocksPrefix ( iB ) ) &    
          // trim ( GI % Stream % CurrentDirectory ) // trim ( Name )
      lMultiMeshName  ( iB ) = len_trim ( MultiMeshName ( iB ) )
      MultiMeshType ( iB ) = GI % MultiMeshType
    end do
    
    Error = DBSET2DSTRLEN ( LDF )
    Error = DBPUTMMESH &
              ( GI % Stream % MultiMeshHandle, trim ( Name ), &
                len_trim ( Name ), nBlocks, MultiMeshName, &
                lMultiMeshName, MultiMeshType, SiloOptionList, Error )

    if ( SiloOptionList /= DB_F77NULL ) &
      Error = DBFREEOPTLIST ( SiloOptionList )
      
  end subroutine WriteMultiMesh
  
  
  subroutine WriteMultiVariable ( GI, Name, TimeOption, CycleNumberOption )
      
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
      
    integer ( KDI ) :: &
      iB, &
      nBlocks, &
      nSiloOptions, &
      SiloOptionList, &
      Error
    integer ( KDI ), dimension ( : ), allocatable :: &
      lMultiVariableName, &
      MultiVariableType
    character ( LDF ), dimension ( : ), allocatable :: &
      MultiVariableName
    
    if ( .not. GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) ) &
      return
      
    nBlocks = GI % Stream % nBlocks
    allocate ( MultiVariableName ( nBlocks ) )
    allocate ( lMultiVariableName ( nBlocks ) )
    allocate ( MultiVariableType ( nBlocks ) )
      
    nSiloOptions = 0 
    SiloOptionList = DB_F77NULL
    
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
      
    if ( nSiloOptions > 0 ) then
      Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
      if ( present ( CycleNumberOption ) ) &
        Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
      if ( present ( TimeOption ) ) &
        Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
    end if
    
    do iB = 1, nBlocks
      MultiVariableName ( iB ) & 
        = trim ( GI % Stream % MeshBlocksPrefix ( iB ) ) &    
          // trim ( GI % Stream % CurrentDirectory ) // trim ( Name )
      lMultiVariableName  ( iB ) = len_trim ( MultiVariableName ( iB ) )
      MultiVariableType ( iB ) = GI % MultiVariableType
    end do
    
    Error = DBSET2DSTRLEN ( LDF )
    Error = DBPUTMVAR &
              ( GI % Stream % MultiMeshHandle, trim ( Name ), &
                len_trim ( Name ), nBlocks, MultiVariableName, &
                lMultiVariableName,  MultiVariableType, SiloOptionList, &
                Error )

    if ( SiloOptionList /= DB_F77NULL ) &
      Error = DBFREEOPTLIST ( SiloOptionList )
    
  end subroutine WriteMultiVariable
  
  
  subroutine WriteVectorVariable ( GI, S )
  
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    class ( StorageForm ), intent ( in ) :: &
      S
    
    integer ( KDI ) :: &
      iD, &
      iVctr, &
      oS, &  !-- oString
      Error
    integer ( KDI ), dimension ( : ), allocatable :: &
      lVectorName, &
      lVectorDefinition
    character ( LDB ), dimension ( : ), allocatable :: &
      VectorName, &
      VectorDefinition, &
      Scratch
    
    
    if ( .not. GI % Stream % IsWritable ( CheckMultiMeshOption = .true. ) ) &
      return
    
    if ( S % nVectors == 0 ) return
    
    allocate ( VectorName ( S % nVectors ) )
    allocate ( VectorDefinition ( S % nVectors ) )
    allocate ( lVectorName ( S % nVectors ) )
    allocate ( lVectorDefinition ( S % nVectors ) )
    
    allocate ( Scratch ( GI % nDimensions ) )
    
    oS = 0
    if ( GI % Stream % CurrentDirectory ( 1 : 1 ) == '/' ) oS = 1
    
    do iVctr = 1, S % nVectors
      !-- FIXME: XL compiler has trouble using implicit do loop
      do iD = 1, GI % nDimensions
        Scratch ( iD ) &
          = ( '<' // trim ( GI % Stream % CurrentDirectory ( oS + 1 : ) ) &
              // trim ( S % Variable ( S % VectorIndices ( iVctr ) &
                                         % Value ( iD ) ) ) // '>' )
      end do
      call Join ( Scratch , ',', VectorDefinition ( iVctr ) ) 
      VectorDefinition ( iVctr ) &
        = '{' // trim ( VectorDefinition ( iVctr ) ) // '}'
      lVectorDefinition ( iVctr ) = len_trim ( VectorDefinition ( iVctr ) )
      VectorName ( iVctr ) &
        = trim ( GI % Stream % CurrentDirectory ( oS + 1 : ) ) &
            // trim ( S % Vector ( iVctr ) )
      lVectorName ( iVctr ) = len_trim ( VectorName ( iVctr ) )
    end do
    Error = DBSET2DSTRLEN ( LDB )
    Error = DBPUTDEFVARS &
              ( GI % Stream % MultiMeshHandle, trim ( S % Name ), &
                S % lName, S % nVectors, VectorName, lVectorName, &
                [ ( DB_VARTYPE_VECTOR, iVctr = 1, S % nVectors ) ], &
                VectorDefinition, lVectorDefinition, &
                [ ( DB_F77NULL, iVctr = 1, S % nVectors ) ], Error )
   
  end subroutine WriteVectorVariable
  
  
  subroutine ReadHeader ( GI, TimeOption, CycleNumberOption )
  
    class ( GridImageSiloTemplate ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      lHeaderName
    character ( LDF ) :: &
      WorkingDirectory
    character ( kind = c_char, len = LDF ) :: &
      HeaderName
    type ( c_ptr ) :: &
      DB_File      = c_null_ptr, &
      DB_QM_Handle = c_null_ptr
    type ( DB_QuadMeshType ), pointer :: &
      DB_QM
      
    WorkingDirectory = GI % Stream % CurrentDirectory
    call GI % Stream % ChangeDirectory ( '/' )

    HeaderName = 'Header' // c_null_char
    lHeaderName = len_trim ( HeaderName )

    DB_File &
      = GI % Stream % AccessSiloPointer ( GI % Stream % MeshBlockHandle )
    DB_QM_Handle = DB_GetQuadMesh ( DB_File, HeaderName ) 
    call c_f_pointer ( DB_QM_Handle, DB_QM )
    
    if ( present ( CycleNumberOption ) ) &
      CycleNumberOption = DB_QM % CycleNumber
      
    if ( present ( TimeOption ) ) &
!      TimeOption = DB_QM % D_Time
!-- FIXME: What about the units?
      call TimeOption % Initialize ( '', DB_QM % TimeInDoublePrecision )
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )
    
    nullify ( DB_QM )
    
  end subroutine ReadHeader
  
  
end module GridImageSilo_Template
