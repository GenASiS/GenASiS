!-- PointGridImageForm is a class for handling point grid 
!   in Silo format, suitable for visualization tool such as VisIt.

module PointGridImage_Form
  
  use iso_c_binding
  use VariableManagement
  use Display
  use GridImageBasics
  
  implicit none 
  private
  
  include 'silo_f9x.inc'
  
  type, public, extends ( GridImageSiloTemplate ) :: PointGridImageForm 
    integer ( KDI ) :: &
      nCells = 0
  contains
    procedure, public, pass :: &
      SetGrid
    procedure, public, pass :: &
      SetReadAttributes
    procedure, public, pass :: &
      Write
    procedure, private, pass :: &
      WriteMesh
    procedure, private, pass :: &
      WriteStorage
    procedure, public, pass :: &
      Read
    procedure, private, pass :: &
      ReadMesh
    procedure, private, pass :: &
      ReadStorage
    final :: &
      Finalize
  end type PointGridImageForm
    
contains


  subroutine SetGrid &
               ( PGI, Directory, Coordinate, nCells, nDimensions, oValue, &
                 CoordinateLabelOption, CoordinateUnitOption )

    class ( PointGridImageForm ), intent ( inout ), target :: &
      PGI
    character ( * ), intent ( in ) :: &
      Directory    
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Coordinate
    integer ( KDI ), intent ( in ) :: &
      nCells, &
      nDimensions, &
      oValue
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption

    PGI % oValue      = oValue
    PGI % nDimensions = nDimensions
    PGI % nTotalCells = nCells
    PGI % nGhostCells = 0
    PGI % lDirectory  = len_trim ( Directory )

    PGI % MultiMeshType = DB_POINTMESH
    PGI % MultiVariableType = DB_POINTVAR

    if ( allocated ( PGI % NodeCoordinate_1 ) ) &
      deallocate ( PGI % NodeCoordinate_1 )
    if ( allocated ( PGI % NodeCoordinate_2 ) ) &
      deallocate ( PGI % NodeCoordinate_2 )
    if ( allocated ( PGI % NodeCoordinate_3 ) ) &
      deallocate ( PGI % NodeCoordinate_3 )
    allocate ( PGI % NodeCoordinate_1 ( nCells ) )
    allocate ( PGI % NodeCoordinate_2 ( nCells ) )
    allocate ( PGI % NodeCoordinate_3 ( nCells ) )
    PGI % NodeCoordinate_1 = Coordinate ( oValue + 1 : oValue + nCells, 1 )
    PGI % NodeCoordinate_2 = Coordinate ( oValue + 1 : oValue + nCells, 2 )
    if ( nDimensions > 2 ) then
      PGI % NodeCoordinate_3 &
        = Coordinate ( oValue + 1 : oValue + nCells, 3 )
    else
      PGI % NodeCoordinate_3 = 0.0_KDR
    end if

    PGI % nCells = nCells
    
    PGI % Directory = Directory
    if ( trim ( PGI % Directory ) == '/' ) PGI % Directory = ''
    
    if ( present ( CoordinateUnitOption ) ) &
      PGI % CoordinateUnit ( 1 : size ( CoordinateUnitOption ) ) &
        = CoordinateUnitOption
    
    if ( present ( CoordinateLabelOption ) ) &
      PGI % CoordinateLabel ( 1 : size ( CoordinateLabelOption ) ) &
        = CoordinateLabelOption

  end subroutine SetGrid


  subroutine SetReadAttributes ( GI, Directory, oValue )
  
    class ( PointGridImageForm ), intent ( inout ) :: &
      GI 
    character ( * ), intent ( in ) :: &
      Directory
    integer ( KDI ), intent ( in ) :: &
      oValue
      
    GI % oValue      = oValue
    
    GI % lDirectory  = len_trim ( Directory )
    
    GI % Directory   = Directory
    
  end subroutine SetReadAttributes
  
  
  subroutine Write ( GI, TimeOption, CycleNumberOption )
  
    class ( PointGridImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( in ) , optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption

    character ( LDF ) :: &
      WorkingDirectory
      
    if ( .not. GI % Stream % IsWritable ( ) ) return  
    
    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % MakeDirectory ( GI % Directory ) 
    
    call GI % WriteHeader ( TimeOption, CycleNumberOption )
    
    call GI % WriteMesh ( TimeOption, CycleNumberOption )
    
    call GI % WriteStorage ( TimeOption, CycleNumberOption ) 
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory ) 
  
  end subroutine Write
  
  
  subroutine Read ( GI, TimeOption, CycleNumberOption ) 
    
    class ( PointGridImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iStrg, &
      nVariables 
    character ( LDL ), dimension ( : ), allocatable :: &
      StorageName, &
      VariableName
    character ( LDF ) :: &
      WorkingDirectory
    
    if ( .not. GI % Stream % IsReadable ( ) ) return    

    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % ChangeDirectory ( GI % Directory )  

    call GI % ReadHeader ( TimeOption, CycleNumberOption )

    call GI % ReadMesh ( )
    
    !-- prepare Storage to read into
    if ( GI % nStorages == 0 ) then 
      call GI % Stream % ListContents ( ContentTypeOption = 'Directory' )
      GI % nStorages = size ( GI % Stream % ContentList )
!-- FIXME: NAG 5.3.1 should support sourced allocation
!      allocate ( StorageName, source = GI % Stream % ContentList )
      allocate ( StorageName ( size ( GI % Stream % ContentList ) ) )
      StorageName = GI % Stream % ContentList
      do iStrg = 1, GI % nStorages
        if ( len_trim ( StorageName ( iStrg ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( StorageName ( iStrg ) )
        call GI % Stream % ListContents &
               ( ContentTypeOption = 'PointGridVariable' )
        if ( allocated ( VariableName ) ) deallocate ( VariableName )
!        allocate ( VariableName, source = GI % Stream % ContentList )
        allocate ( VariableName ( size ( GI % Stream % ContentList ) ) )
        VariableName = GI % Stream % ContentList 
        nVariables = size ( GI % Stream % ContentList )
        if ( nVariables == 0 ) then
          GI % nStorages = 0 
        else
          call GI % Storage ( iStrg ) % Initialize &
                 ( [ GI % nCells, nVariables ], &
                   VariableOption = VariableName, &
                   NameOption = StorageName ( iStrg ) ) 
        end if
        if ( len_trim ( StorageName ( iStrg ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( '..' )
      end do
    end if  
      
    call GI % ReadStorage ( )
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )

  end subroutine Read

  
  impure elemental subroutine Finalize ( PGI )
    
    type ( PointGridImageForm ), intent ( inout ) :: & 
      PGI 
    
    nullify ( PGI % Stream )

    if ( allocated ( PGI % Storage ) ) &
      deallocate ( PGI % Storage )
    
    if ( allocated ( PGI % NodeCoordinate_3 ) ) &
      deallocate ( PGI % NodeCoordinate_3 )
    if ( allocated ( PGI % NodeCoordinate_2 ) ) &
      deallocate ( PGI % NodeCoordinate_2 )
    if ( allocated ( PGI % NodeCoordinate_1 ) ) &
      deallocate ( PGI % NodeCoordinate_1 )
      
  end subroutine Finalize 
  
  
  subroutine WriteMesh ( PGI, TimeOption, CycleNumberOption )
    
    class ( PointGridImageForm ), intent ( inout ) :: &
      PGI
    type ( MeasuredValueForm ), intent ( in ) , optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iD, &     !-- iDimension
      nSiloOptions, &
      SiloOptionList, &
      Error
      
    nSiloOptions = 3
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
    do iD = 1, 3
      if ( len_trim ( PGI % CoordinateUnit ( iD ) % Label ) > 0 ) &
        nSiloOptions = nSiloOptions + 1
    end do
      
    Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
    
    Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_XLABEL, &
                  trim ( PGI % CoordinateLabel ( 1 ) ), &
                  len_trim ( PGI % CoordinateLabel ( 1 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_YLABEL, &
                trim ( PGI % CoordinateLabel ( 2 ) ), &
                len_trim ( PGI % CoordinateLabel ( 2 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_ZLABEL, &
                trim ( PGI % CoordinateLabel ( 3 ) ), &
                len_trim ( PGI % CoordinateLabel ( 3 ) ) )
    
    if ( len_trim ( PGI % CoordinateUnit ( 1 ) % Label ) > 0 ) &
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_XUNITS, &
                  trim ( PGI % CoordinateUnit ( 1 ) % Label ), &
                  len_trim ( PGI % CoordinateUnit ( 1 ) % Label ) )
    if ( len_trim ( PGI % CoordinateUnit ( 2 ) % Label ) > 0 ) &   
      Error = DBADDCOPT & 
                ( SiloOptionList, DBOPT_YUNITS, &
                  trim ( PGI % CoordinateUnit ( 2 ) % Label ), &
                  len_trim ( PGI % CoordinateUnit ( 2 ) % Label ) )
    if ( len_trim ( PGI % CoordinateUnit ( 3 ) % Label ) > 0 ) &   
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_ZUNITS, &
                  trim ( PGI % CoordinateUnit ( 3 ) % Label ), &
                  len_trim ( PGI % CoordinateUnit ( 3 ) % Label ) )

    if ( present ( CycleNumberOption ) ) &
      Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
    if ( present ( TimeOption ) ) &
      Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
    
    Error = DBPUTPM &
              ( PGI % Stream % MeshBlockHandle, 'Mesh', 4, &
                PGI % nDimensions, &
                PGI % NodeCoordinate_1 / PGI % CoordinateUnit ( 1 ) % Number, &
                PGI % NodeCoordinate_2 / PGI % CoordinateUnit ( 2 ) % Number, &
                PGI % NodeCoordinate_3 / PGI % CoordinateUnit ( 3 ) % Number, &
                PGI % nCells, DB_DOUBLE, SiloOptionList, Error )
    
    Error = DBFREEOPTLIST ( SiloOptionList )
    
    call PGI % WriteMultiMesh ( 'Mesh', TimeOption, CycleNumberOption )
  
  end subroutine WriteMesh
  
  
  subroutine ReadMesh ( PGI )
    
    class ( PointGridImageForm ), intent ( inout ) :: &
      PGI
      
    real ( c_double ), dimension ( : ), pointer :: &
      NC_1, &
      NC_2, &
      NC_3   
    type ( c_ptr ) :: &
      DB_File, &
      DB_PM_Handle
    type ( DB_PointMeshType ), pointer :: &
      DB_PM
    
    DB_File &
      = PGI % Stream % AccessSiloPointer ( PGI % Stream % MeshBlockHandle )
    DB_PM_Handle = DB_GetPointMesh ( DB_File, c_char_'Mesh' // c_null_char )
    call c_f_pointer ( DB_PM_Handle, DB_PM )
    
    PGI % nDimensions = DB_PM % nDimensions
    PGI % nTotalCells = DB_PM % nElements
    PGI % nCells      = DB_PM % nElements
    PGI % nGhostCells = 0
    
    if ( DB_PM % nDimensions >= 1 ) then
      call c_f_pointer &
             ( DB_PM % NodeCoordinate ( 1 ), NC_1, [ DB_PM % nElements ] )
      if ( allocated ( PGI % NodeCoordinate_1 ) ) &  
        deallocate ( PGI % NodeCoordinate_1 )
      allocate &
        ( PGI % NodeCoordinate_1 ( DB_PM % nElements ), &
          Source = NC_1 ( 1 : DB_PM % nElements ) )
    else 
      allocate ( PGI % NodeCoordinate_1 ( 0 ) )
    end if
    
    if ( DB_PM % nDimensions >= 2 ) then
      call c_f_pointer &
             ( DB_PM % NodeCoordinate ( 2 ), NC_2, &
               [ DB_PM % nElements ] )
      if ( allocated ( PGI % NodeCoordinate_2 ) ) &  
        deallocate ( PGI % NodeCoordinate_2 )
      allocate &
        ( PGI % NodeCoordinate_2 ( DB_PM % nElements ), &
          Source = NC_2 ( 1 : DB_PM % nElements ) )
    else 
      allocate ( PGI % NodeCoordinate_2 ( 0 ) )
    end if

    if ( DB_PM % nDimensions == 3 ) then
      call c_f_pointer &
             ( DB_PM % NodeCoordinate ( 3 ), NC_3, &
               [ DB_PM % nElements ] )
      if ( allocated ( PGI % NodeCoordinate_3 ) ) &  
        deallocate ( PGI % NodeCoordinate_3 )
      allocate &
        ( PGI % NodeCoordinate_3 ( DB_PM % nElements ), &
          Source = NC_3 ( 1 : DB_PM % nElements ) )
    else 
      allocate ( PGI % NodeCoordinate_3 ( 0 ) )
    end if

    call DB_FreePointMesh ( DB_PM_Handle )
    
    !-- FIXME: Here we make the assumption that the CoordinateUnit for reading 
    !          is the same as the ones used to write. A better way would be to 
    !          read the unit directly from Silo attribute (but that's harder   
    !          to do esp. with bugs on C compatibility for CCE compiler, hence 
    !          more wrappers would had to be written)
    PGI % NodeCoordinate_1 &
      = PGI % NodeCoordinate_1 * PGI % CoordinateUnit ( 1 ) % Number
    PGI % NodeCoordinate_2 &
      = PGI % NodeCoordinate_2 * PGI % CoordinateUnit ( 2 ) % Number
    PGI % NodeCoordinate_3 &
      = PGI % NodeCoordinate_3 * PGI % CoordinateUnit ( 3 ) % Number
    
    DB_PM_Handle = c_null_ptr
    DB_File      = c_null_ptr

    nullify ( NC_3, NC_2, NC_1 )
    nullify ( DB_PM )

  end subroutine ReadMesh
  
  
  subroutine WriteStorage ( PGI, TimeOption, CycleNumberOption ) 

    class ( PointGridImageForm ), intent ( inout ) :: &
      PGI
    type ( MeasuredValueForm ), intent ( in ) , optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
  
    integer ( KDI ) :: &
      iV, &      !-- iValue
      iVrbl, &   !-- iVariable
      iS, &      !-- iSelected
      iStrg, &      !-- iStorage
      nSiloOptions, &
      SiloOptionList, &
      Error
    real ( KDR ), dimension ( : ), pointer :: &
      Value
    character ( LDF ) :: &
      MeshDirectory
  
    SiloOptionList = DB_F77NULL
  
    MeshDirectory = PGI % Stream % CurrentDirectory
  
    do iStrg = 1, PGI % nStorages
    
      associate ( S => PGI % Storage ( iStrg ) )
    
      call Show ( 'Writing a Storage (point)', CONSOLE % INFO_5 )
      call Show ( iStrg, 'iStorage', CONSOLE % INFO_5 )
      call Show ( S % Name, 'Name', CONSOLE % INFO_5 )

      call PGI % Stream % MakeDirectory ( S % Name ) 
    
      do iS = 1, S % nVariables
        
        iVrbl = S % iaSelected ( iS )
        
        call Show ( 'Writing a Variable (point)', CONSOLE % INFO_6 )
        call Show ( iS, 'iSelected', CONSOLE % INFO_6 )
        call Show ( S % Variable ( iVrbl ), 'Name', CONSOLE % INFO_6 )

        Value => S % Value ( PGI % oValue + 1 :, iVrbl )
        
        nSiloOptions = 0
        if ( present ( TimeOption ) ) &
          nSiloOptions = nSiloOptions + 1
        if ( present ( CycleNumberOption ) ) &
          nSiloOptions = nSiloOptions + 1    
        if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) &
          nSiloOptions = nSiloOptions + 1
          
        if ( nSiloOptions > 0 ) &
          Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
          
        if ( present ( TimeOption ) ) &
          Error = DBADDDOPT &
                    ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
        if ( present ( CycleNumberOption ) ) &
          Error = DBADDIOPT &
                    ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
        if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) &
          Error = DBADDCOPT &
                    ( SiloOptionList, DBOPT_UNITS, &
                      trim ( S % Unit ( iVrbl ) % Label ), &
                      len_trim ( S % Unit ( iVrbl ) % Label ) )
        
        call Show &
               ( trim ( S % Variable ( iVrbl ) ), 'Variable', &
                 CONSOLE % INFO_6 )
        call Show &
               ( S % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
        call Show &
               ( trim ( MeshDirectory ) // 'Mesh', 'MeshDirectory', &
                 CONSOLE % INFO_6 )
        call Show &
               ( len_trim ( MeshDirectory ) + 4, 'lDirectory', &
                 CONSOLE % INFO_6 )
        call Show ( nSiloOptions, 'nSiloOptions', CONSOLE % INFO_6 )
        if ( len_trim ( S % Unit ( iVrbl ) % Label ) > 0 ) then
          call Show &
                 ( trim ( S % Unit ( iVrbl ) % Label ), 'Unit', &
                   CONSOLE % INFO_6 )
          call Show &
                 ( len_trim ( S % Unit ( iVrbl ) % Label ), 'lUnit', &
                   CONSOLE % INFO_6 )
        end if

        Error = DBPUTPV1 &
                  ( PGI % Stream % MeshBlockHandle, &
                    trim ( S % Variable ( iVrbl ) ), &
                    S % lVariable ( iVrbl ), &
                    trim ( MeshDirectory ) // 'Mesh', &
                    len_trim ( MeshDirectory ) + 4, &
                    Value / S % Unit ( iVrbl ) % Number, PGI % nCells, &
                    DB_DOUBLE, SiloOptionList, Error )

        if ( SiloOptionList /= DB_F77NULL ) then
          call Show &
                 ( SiloOptionList, 'SiloOptionList before free', &
                   CONSOLE % INFO_6 )        
          Error = DBFREEOPTLIST ( SiloOptionList )
          call Show &
                 ( SiloOptionList, 'SiloOptionList after free', &
                   CONSOLE % INFO_6 )        
          SiloOptionList = DB_F77NULL
          call Show &
                 ( SiloOptionList, 'SiloOptionList after nullification', &
                   CONSOLE % INFO_6 )        
          call Show ( DB_F77NULL, 'DB_F77NULL', CONSOLE % INFO_6 )
        end if
          
        call PGI % WriteMultiVariable &
               ( S % Variable ( iVrbl ), TimeOption, CycleNumberOption )
        
      end do
      
      call PGI % WriteVectorVariable ( S )
      
      call PGI % Stream % ChangeDirectory ( '../' )
      
      end associate 
  
    end do
  
    nullify ( Value )
  
  end subroutine WriteStorage 
  
  
  subroutine ReadStorage ( PGI ) 

    class ( PointGridImageForm ), intent ( inout ) :: &
      PGI
      
    integer ( KDI ) :: &
      iV, &      !-- iValue
      iVrbl, &   !-- iVariable
      iStrg, &      !-- iSelected
      iS, &      !-- iStorage   
      iA, &      !-- iArray   
      oV, &
      nProperCells
    real ( c_double ), dimension ( : ), pointer :: &
      VariableValue
    character ( LDF ) :: &
      MeshDirectory
    character ( kind = c_char, len = LDL ) :: &
      VariableName
    type ( c_ptr ) :: &
      DB_File, &
      DB_MV_Handle
    type ( c_ptr ), dimension ( : ), pointer :: &
      ValueArrays
    type ( DB_MeshVariableType ), pointer :: &
      DB_MV

    MeshDirectory = PGI % Stream % CurrentDirectory

    do iStrg = 1, PGI % nStorages
      
      associate ( S => PGI % Storage ( iStrg ), &
                  nDims => PGI % nDimensions )

      call Show ( 'Reading a Storage', CONSOLE % INFO_5 )
      call Show ( iStrg, 'iStorage', CONSOLE % INFO_5 )
      call Show ( S % Name, 'Name', CONSOLE % INFO_5 )

      if ( len_trim ( S % Name ) > 0 ) &
        call PGI % Stream % ChangeDirectory ( S  % Name )

      DB_File &
        = PGI % Stream % AccessSiloPointer ( PGI % Stream % MeshBlockHandle )

      do iS = 1, S % nVariables

        iVrbl = S % iaSelected ( iS )
        VariableName = trim ( S % Variable ( iVrbl ) ) // c_null_char
        
        DB_MV_Handle = DB_GetPointVariable ( DB_File, VariableName )
        call c_f_pointer ( DB_MV_Handle, DB_MV )
      
        call c_f_pointer ( DB_MV % Value, ValueArrays, [ DB_MV % nValues ] )
        
        do iA = 1, DB_MV % nValues
          !-- FIXME: An assumption is made that the unit used to write
          !          and read are the same. A better way would be to read
          !          the unit directly from Silo file.
          call c_f_pointer &
                 ( ValueArrays ( iA ), VariableValue, [ DB_MV % nElements ] )
          call Copy &
                 ( VariableValue * S % Unit ( iVrbl ) % Number, &
                   S % Value ( PGI % oValue + 1 &
                                : PGI % oValue + DB_MV % nElements, iVrbl ) )
        end do
        
        call DB_FreeMeshVariable ( DB_MV_Handle )
      
      end do
      
      if ( len_trim ( S % Name ) > 0 ) &
        call PGI % Stream % ChangeDirectory ( '../' )

      end associate

    end do

    DB_File      = c_null_ptr
    DB_MV_Handle = c_null_ptr
         
  end subroutine ReadStorage


end module PointGridImage_Form
