!-- UnstructuredGridImageForm is a class for handling unstructured grid 
!   in Silo format, suitable for visualization tool such as VisIt.

module UnstructuredGridImage_Form
  
  use iso_c_binding
  use VariableManagement
  use Display
  use GridImageBasics

  implicit none
  private

  include 'silo_f9x.inc'

  type, public, extends ( GridImageSiloTemplate ) &
    :: UnstructuredGridImageForm
      integer ( KDI ) :: &
        nNodes          = 0, &
        NodeListSize    = 0
      integer ( KDI ), dimension ( : ), allocatable :: &
        NodeList
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
      WriteVariableGroup
    procedure, public, pass :: &
      Read
    procedure, private, pass :: &
      ReadMesh
    procedure, private, pass :: &
      ReadVariableGroup 
    procedure, public, pass :: &
      ClearGrid
    final :: &
      Finalize
  end type UnstructuredGridImageForm


contains

  
  subroutine SetGrid &
               ( UGI, Directory, NodeCoordinate, &
                 NodeConnectivity, nDimensions, nProperCells, nGhostCells, &
                 oValue, CoordinateUnitOption )

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI
    character ( * ), intent ( in ) :: &
      Directory
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      NodeCoordinate
    integer ( KDI ), dimension ( :, : ), intent ( in ) :: &
      NodeConnectivity
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nProperCells, &
      nGhostCells, &
      oValue
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption

    UGI % oValue      = oValue
    UGI % nDimensions = nDimensions
    UGI % nTotalCells = nProperCells + nGhostCells
    UGI % nGhostCells = nGhostCells
    UGI % nNodes      = size ( NodeCoordinate, dim = 1 )
    UGI % lDirectory  = len_trim ( Directory )
    
    UGI % MultiMeshType = DB_UCDMESH
    UGI % MultiVariableType = DB_UCDVAR

    if ( UGI % nTotalCells > 0 ) then
      UGI % NodeListSize = 2 ** nDimensions * UGI % nTotalCells
      allocate ( UGI % NodeList ( UGI % NodeListSize ) )
      UGI % NodeList &
        = reshape &
            ( NodeConnectivity &
                ( 1 : 2**nDimensions, 1 : UGI % nTotalCells ), & 
                  [ UGI % NodeListSize ] )
    else !-- nTotalCells == 0
      UGI % NodeListSize = 2**nDimensions
      allocate ( UGI % NodeList ( UGI % NodeListSize ) )
      UGI % NodeList &
        = reshape &
            ( NodeConnectivity &
                ( 1 : 2 ** nDimensions, 1 ), &
                  [ UGI % NodeListSize ] )
    end if
      
    allocate ( UGI % NodeCoordinate_1 ( UGI % nNodes ) )
    allocate ( UGI % NodeCoordinate_2 ( UGI % nNodes ) )
    allocate ( UGI % NodeCoordinate_3 ( UGI % nNodes ) )
    UGI % NodeCoordinate_1 = NodeCoordinate ( :, 1 )
    UGI % NodeCoordinate_2 = NodeCoordinate ( :, 2 )
    UGI % NodeCoordinate_3 = NodeCoordinate ( :, 3 )

    UGI % Directory = Directory
    if ( trim ( UGI % Directory ) == '/' ) UGI % Directory = ''

    if ( present ( CoordinateUnitOption ) ) &
      UGI % CoordinateUnit ( 1 : size ( CoordinateUnitOption ) ) &
        = CoordinateUnitOption
    
  end subroutine SetGrid
  
  
  subroutine SetReadAttributes ( GI, Directory, oValue )

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
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

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( in ), optional :: &
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
    
    call GI % WriteVariableGroup ( TimeOption, CycleNumberOption )

    call GI % Stream % ChangeDirectory ( WorkingDirectory )
    
  end subroutine Write 
  
  
  subroutine Read ( GI, TimeOption, CycleNumberOption )
  
    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      GI
    type ( MeasuredValueForm ), intent ( out ), optional :: &
      TimeOption
    integer ( KDI ), intent ( out ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iG, &
      nVariables 
    character ( LDL ), dimension ( : ), allocatable :: &
      GroupName, &
      VariableName
    character ( LDF ) :: &
      WorkingDirectory
      
    if ( .not. GI % Stream % IsReadable ( ) ) return
    
    WorkingDirectory = GI % Stream % CurrentDirectory
    if ( GI % lDirectory > 0 ) &
      call GI % Stream % ChangeDirectory ( GI % Directory )
      
    call GI % ReadHeader ( TimeOption, CycleNumberOption )
    
    call GI % ReadMesh ( )
    
    !-- prepare VariableGroup to read into
    if ( GI % nVariableGroups == 0 ) then
      call GI % Stream % ListContents ( ContentTypeOption = 'Directory' )
      GI % nVariableGroups = size ( GI % Stream % ContentList )
!-- FIXME: NAG 5.3.1 should support sourced allocation
!      allocate ( GroupName, source = GI % Stream % ContentList )
      allocate ( GroupName ( size ( GI % Stream % ContentList ) ) )
      GroupName = GI % Stream % ContentList
      do iG = 1, GI % nVariableGroups
        if ( len_trim ( GroupName ( iG ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( GroupName ( iG ) )
        call GI % Stream % ListContents &
               ( ContentTypeOption = 'UnstructuredGridVariable' )
        if ( allocated ( VariableName ) ) deallocate ( VariableName )
!        allocate ( VariableName, source = GI % Stream % ContentList )
        allocate ( VariableName ( size ( GI % Stream % ContentList ) ) )
        VariableName = GI % Stream % ContentList 
        nVariables = size ( GI % Stream % ContentList )
        if ( nVariables == 0 ) then
          GI % nVariableGroups = 0
        else
          call GI % VariableGroup ( iG ) % Initialize &
                 ( [ GI % nTotalCells, nVariables ], &
                   VariableOption = VariableName, & 
                   NameOption = GroupName ( iG ) )
        end if
        if ( len_trim ( GroupName ( iG ) ) > 0 ) &
          call GI % Stream % ChangeDirectory ( '..' )
      end do
    end if
      
    call GI % ReadVariableGroup ( )
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory )

  end subroutine Read 
  
  
  subroutine ClearGrid ( UGI )

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI

    if ( allocated ( UGI % NodeCoordinate_3 ) ) &
      deallocate ( UGI % NodeCoordinate_3 )
    if ( allocated ( UGI % NodeCoordinate_2 ) ) &
      deallocate ( UGI % NodeCoordinate_2 )
    if ( allocated ( UGI % NodeCoordinate_1 ) ) &
      deallocate ( UGI % NodeCoordinate_1 )

    if ( allocated ( UGI % NodeList ) ) deallocate ( UGI % NodeList )

  end subroutine ClearGrid


  impure elemental subroutine Finalize ( UGI )

    type ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI

    nullify ( UGI % Stream )

    if ( allocated ( UGI % VariableGroup ) ) deallocate ( UGI % VariableGroup )

    call UGI % ClearGrid ( )

  end subroutine Finalize
  
  
  subroutine WriteMesh ( UGI, TimeOption, CycleNumberOption )

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iD, &     !-- iDimension
      nTotalCells, &
      nSiloOptions, &
      SiloZoneType, &
      SiloOptionList, &
      Error
    
    select case ( UGI % nDimensions )
      case ( 1 )
        SiloZoneType = DB_F77NULL
      case ( 2 )
        SiloZoneType = DB_ZONETYPE_QUAD
      case ( 3 )
        SiloZoneType = DB_ZONETYPE_HEX
    end select
    
    if ( UGI % nTotalCells > 0 ) then
      Error = DBPUTZL2 (  &
                UGI % Stream % MeshBlockHandle, 'ZoneList', 8, &
                UGI % nTotalCells, UGI % nDimensions, UGI % NodeList, &
                UGI % NodeListSize, 1, 0, UGI % nGhostCells, &
                [ SiloZoneType ], [ 2 ** UGI % nDimensions ],&
                [ UGI % nTotalCells ], 1, DB_F77NULL, Error )
    else
      Error = DBPUTZL2 (  &
                UGI % Stream % MeshBlockHandle, 'ZoneList', 8, &
                1, UGI % nDimensions, UGI % NodeList, UGI % NodeListSize, &
                1, 0, UGI % nGhostCells, [ SiloZoneType ], &
                [ 2 ** UGI % nDimensions ], [ 1 ], 1, DB_F77NULL, Error )
    end if
    
    nSiloOptions = 3
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
    do iD = 1, 3
      if ( len_trim ( UGI % CoordinateUnit ( iD ) % Label ) > 0 ) &
        nSiloOptions = nSiloOptions + 1
    end do
    
    Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
    
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_XLABEL, &
                trim ( UGI % CoordinateLabel ( 1 ) ), &
                len_trim ( UGI % CoordinateLabel ( 1 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_YLABEL, &
                trim ( UGI % CoordinateLabel ( 2 ) ), &
                len_trim ( UGI % CoordinateLabel ( 2 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_ZLABEL, &
                trim ( UGI % CoordinateLabel ( 3 ) ), &
                len_trim ( UGI % CoordinateLabel ( 3 ) ) )
    
    if ( len_trim ( UGI % CoordinateUnit ( 1 ) % Label ) > 0 ) &
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_XUNITS, &
                  trim ( UGI % CoordinateUnit ( 1 ) % Label ), &
                  len_trim ( UGI % CoordinateUnit ( 1 ) % Label ) )
    if ( len_trim ( UGI % CoordinateUnit ( 2 ) % Label ) > 0 ) &
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_YUNITS, &
                  trim ( UGI % CoordinateUnit ( 2 ) % Label ), &
                  len_trim ( UGI % CoordinateUnit ( 2 ) % Label ) )
    if ( len_trim ( UGI % CoordinateUnit ( 3 ) % Label ) > 0 ) &
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_ZUNITS, &
                  trim ( UGI % CoordinateUnit ( 3 ) % Label ), &
                  len_trim ( UGI % CoordinateUnit ( 3 ) % Label ) )
    
    if ( present ( CycleNumberOption ) ) &
      Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
    if ( present ( TimeOption ) ) &
      Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
    
    if ( UGI % nTotalCells > 0 ) then
      nTotalCells = UGI % nTotalCells
    else 
      nTotalCells = 1
    end if
      
    Error &
      = DBPUTUM &
          ( UGI % Stream % MeshBlockHandle, 'Mesh', 4, UGI % nDimensions, &
            UGI % NodeCoordinate_1 / UGI % CoordinateUnit ( 1 ) % Number, &
            UGI % NodeCoordinate_2 / UGI % CoordinateUnit ( 2 ) % Number, &
            UGI % NodeCoordinate_3 / UGI % CoordinateUnit ( 3 ) % Number, &
            trim ( UGI % CoordinateLabel ( 1 ) ), &
            len_trim ( UGI % CoordinateLabel ( 1 ) ), &
            trim ( UGI % CoordinateLabel ( 2 ) ), &
            len_trim ( UGI % CoordinateLabel ( 2 ) ), &
            trim ( UGI % CoordinateLabel ( 3 ) ), &
            len_trim ( UGI % CoordinateLabel ( 3 ) ), &
            DB_DOUBLE, UGI % nNodes, nTotalCells, &
            'ZoneList', 8, DB_F77NULL, 0, SiloOptionList, Error )
      
    Error = DBFREEOPTLIST ( SiloOptionList ) 
    
    call UGI % WriteMultiMesh ( 'Mesh', TimeOption, CycleNumberOption )
    
  end subroutine WriteMesh
  
  
  subroutine WriteVariableGroup ( UGI, TimeOption, CycleNumberOption )

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
  
    integer ( KDI ) :: &
      iVrbl, &    !-- iVariable
      iS, &       !-- iSelected
      iG, &       !-- iGroup
      nSiloOptions, &
      Error, &
      SiloOptionList
    character ( LDF ) :: &
      MeshDirectory  
  
    MeshDirectory = UGI % Stream % CurrentDirectory 
      
    SiloOptionList = DB_F77NULL
  
    do iG = 1, UGI % nVariableGroups
  
      associate ( VG => UGI % VariableGroup ( iG ) )
    
      call Show ( 'Writing a VariableGroup (unstructured)', CONSOLE % INFO_5 )
      call Show ( iG, 'iGroup', CONSOLE % INFO_5 )
      call Show ( VG % Name, 'Name', CONSOLE % INFO_5 )

      call UGI % Stream % MakeDirectory ( VG % Name )
    
      do iS = 1, VG % nVariables
      
        iVrbl = VG % iaSelected ( iS )
      
        call Show ( 'Writing a Variable (unstructured)', CONSOLE % INFO_6 )
        call Show ( iS, 'iSelected', CONSOLE % INFO_6 )
        call Show ( VG % Variable ( iVrbl ), 'Name', CONSOLE % INFO_6 )

        if ( UGI % nTotalCells > 0 ) then
        
          nSiloOptions = 0
          if ( present ( TimeOption ) ) &
            nSiloOptions = nSiloOptions + 1
          if ( present ( CycleNumberOption ) ) &
            nSiloOptions = nSiloOptions + 1
          if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) &
            nSiloOptions = nSiloOptions + 1
          
          if ( nSiloOptions > 0 ) &
            Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )
          
          if ( present ( TimeOption ) ) &
            Error = DBADDDOPT &
                      ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
          if ( present ( CycleNumberOption ) ) &
            Error = DBADDIOPT &
                      ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
          if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) &
            Error = DBADDCOPT &
                      ( SiloOptionList, DBOPT_UNITS, &
                        trim ( VG % Unit ( iVrbl ) % Label ), &
                        len_trim ( VG % Unit ( iVrbl ) % Label ) )
                  
          call Show &
                 ( trim ( VG % Variable ( iVrbl ) ), 'Variable', &
                   CONSOLE % INFO_6 )
          call Show &
                 ( VG % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
          call Show &
                 ( trim ( MeshDirectory ) // 'Mesh', 'MeshDirectory', &
                   CONSOLE % INFO_6 )
          call Show &
                 ( len_trim ( MeshDirectory ) + 4, 'lDirectory', &
                   CONSOLE % INFO_6 )
          call Show ( nSiloOptions, 'nSiloOptions', CONSOLE % INFO_6 )
          if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) then
            call Show &
                   ( trim ( VG % Unit ( iVrbl ) % Label ), 'Unit', &
                     CONSOLE % INFO_6 )
            call Show &
                   ( len_trim ( VG % Unit ( iVrbl ) % Label ), 'lUnit', &
                     CONSOLE % INFO_6 )
          end if

          Error = DBPUTUV1 &
                    ( UGI % Stream % MeshBlockHandle, &
                      trim ( VG % Variable ( iVrbl ) ), &
                      VG % lVariable ( iVrbl ), &
                      trim ( MeshDirectory ) // 'Mesh', &
                      len_trim ( MeshDirectory ) + 4, &
                      VG % Value &
                        ( UGI % oValue + 1 : UGI % oValue + UGI % nTotalCells, &
                          iVrbl ) / VG % Unit ( iVrbl ) % Number, &
                      UGI % nTotalCells, DB_F77NULL, 0, DB_DOUBLE, &
                      DB_ZONECENT, SiloOptionList, Error )
        
        else !-- UGI % nTotalCells == 0
      
          call Show &
                 ( trim ( VG % Variable ( iVrbl ) ), 'Variable', &
                   CONSOLE % INFO_6 )
          call Show &
                 ( VG % lVariable ( iVrbl ), 'lVariable', CONSOLE % INFO_6 )
          call Show &
                 ( trim ( MeshDirectory ) // 'Mesh', 'MeshDirectory', &
                   CONSOLE % INFO_6 )
          call Show &
                 ( len_trim ( MeshDirectory ) + 4, 'lDirectory', &
                   CONSOLE % INFO_6 )
          call Show ( nSiloOptions, 'nSiloOptions', CONSOLE % INFO_6 )
          if ( len_trim ( VG % Unit ( iVrbl ) % Label ) > 0 ) then
            call Show &
                   ( trim ( VG % Unit ( iVrbl ) % Label ), 'Unit', &
                     CONSOLE % INFO_6 )
            call Show &
                   ( len_trim ( VG % Unit ( iVrbl ) % Label ), 'lUnit', &
                     CONSOLE % INFO_6 )
          end if

          Error = DBPUTUV1 &
                    ( UGI % Stream % MeshBlockHandle, &
                      trim ( VG % Variable ( iVrbl ) ), &
                  VG % lVariable ( iVrbl ), &
                  trim ( MeshDirectory ) // 'Mesh', &
                  len_trim ( MeshDirectory ) + 4, &
                  [ tiny ( 1.0_KDR ) ], 1, DB_F77NULL, 0, DB_DOUBLE, &
                  DB_ZONECENT, SiloOptionList, Error )
      
        end if
      
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
      
        call UGI % WriteMultiVariable  &
               ( VG % Variable ( iVrbl ), TimeOption, CycleNumberOption )
      
      end do
    
      call UGI % WriteVectorVariable ( VG )
    
      call UGI % Stream % ChangeDirectory ( '../' )
    
      end associate
  
    end do
   
  end subroutine WriteVariableGroup
  
  
  subroutine ReadMesh ( UGI )
    
    class ( UnstructuredGridImageForm ), intent ( inout ), target :: &
      UGI
    
    type ( c_ptr ) :: &
      DB_File, &
      DB_UM_Handle
    type ( DB_UnstructuredMeshType ), pointer :: &
      DB_UM
    type ( DB_ZoneListType ), pointer :: &
      DB_ZoneList
    real ( KDR ), dimension ( : ), pointer :: &
      NC_1, &
      NC_2, &
      NC_3
  
    DB_File &
      = UGI % Stream % AccessSiloPointer ( UGI % Stream % MeshBlockHandle )
    DB_UM_Handle &
      = DB_GetUnstructuredMesh ( DB_File, c_char_'Mesh' // c_null_char )
    call c_f_pointer ( DB_UM_Handle, DB_UM )
    
    UGI % nDimensions = DB_UM % nDimensions
    UGI % nNodes = DB_UM % nNodes
    call c_f_pointer ( DB_UM % ZoneList, DB_ZoneList )
    UGI % nTotalCells = DB_ZoneList % nzones
    UGI % nGhostCells &
      = UGI % nTotalCells &
          - ( DB_ZoneList % MaximumIndex - DB_ZoneList % MinimumIndex + 1 )
    call c_f_pointer ( DB_UM % NodeCoordinate ( 1 ), NC_1, [ UGI % nNodes ] )
    call c_f_pointer ( DB_UM % NodeCoordinate ( 2 ), NC_2, [ UGI % nNodes ] )
    if ( UGI % nDimensions == 3 ) &
      call c_f_pointer ( DB_UM % NodeCoordinate ( 3 ), NC_3, [ UGI % nNodes ] )
    
    if ( allocated ( UGI % NodeCoordinate_1 ) ) &
      deallocate ( UGI % NodeCoordinate_1 ) 
    if ( allocated ( UGI % NodeCoordinate_2 ) ) &
      deallocate ( UGI % NodeCoordinate_2 ) 
    if ( allocated ( UGI % NodeCoordinate_3 ) ) &
      deallocate ( UGI % NodeCoordinate_3 )

    allocate ( UGI % NodeCoordinate_1 ( size ( NC_1 ) ) )
    allocate ( UGI % NodeCoordinate_2 ( size ( NC_2 ) ) )
    UGI % NodeCoordinate_1 = NC_1
    UGI % NodeCoordinate_2 = NC_2
    if ( UGI % nDimensions == 3 ) then
      allocate ( UGI % NodeCoordinate_3 ( size ( NC_3 ) ) )
      UGI % NodeCoordinate_3 = NC_3
    end if
      
    call DB_FreeUnstructuredMesh ( DB_UM_Handle )
    
    DB_UM_Handle = c_null_ptr
    DB_File      = c_null_ptr
    
    nullify ( NC_3, NC_2, NC_1 )
    nullify ( DB_ZoneList )
    nullify ( DB_UM )
  
  end subroutine ReadMesh
  
  
  subroutine ReadVariableGroup ( UGI ) 

    class ( UnstructuredGridImageForm ), intent ( inout ) :: &
      UGI
    
    integer ( KDI ) :: &
      iV, &      !-- iValue
      iVrbl, &   !-- iVariable
      iG, &      !-- iSelected
      iS, &      !-- iGroup
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
      DB_UV_Handle
    type ( c_ptr ), dimension ( : ), pointer :: &
      ValueArrays
    type ( DB_UnstructuredVariableType ), pointer :: &
      DB_UV

    MeshDirectory = UGI % Stream % CurrentDirectory
  
    do iG = 1, UGI % nVariableGroups
    
      associate ( VG => UGI % VariableGroup ( iG ), &
                  nDims => UGI % nDimensions )
  
      call Show ( 'Reading a VariableGroup', CONSOLE % INFO_5 )
      call Show ( iG, 'iGroup', CONSOLE % INFO_5 )
      call Show ( VG % Name, 'Name', CONSOLE % INFO_5 )

      if ( len_trim ( VG % Name ) > 0 ) &
        call UGI % Stream % ChangeDirectory ( VG  % Name ) 
    
      DB_File &
        = UGI % Stream % AccessSiloPointer ( UGI % Stream % MeshBlockHandle )
    
      do iS = 1, VG % nVariables
      
        iVrbl = VG % iaSelected ( iS )
        VariableName = trim ( VG % Variable ( iVrbl ) ) // c_null_char

        DB_UV_Handle = DB_GetUnstructuredVariable ( DB_File, VariableName )
        call c_f_pointer ( DB_UV_Handle, DB_UV )
      
        call c_f_pointer ( DB_UV % Value, ValueArrays, [ DB_UV % nValues ] )
        
        oV = 0
        do iA = 1, DB_UV % nValues
          !-- FIXME: An assumption is made that the unit used to write
          !          and read are the same. A better way would be to read
          !          the unit directly from Silo file.
          call c_f_pointer &
                 ( ValueArrays ( iA ), VariableValue, [ DB_UV % nElements ] )
          
          !-- Assume ghost cells is always after proper cells (see call to
          !   DBPUTZL2 in WriteMesh where this assumption is also expressed 
          !   in the argument list)
          call Copy &
                 ( VariableValue ( oV + 1 : ) * VG % Unit ( iVrbl ) % Number, &
                   VG % Value &
                     ( UGI % oValue + oV + 1 &
                         : UGI % oValue + oV + DB_UV % nElements, iVrbl ) )
        end do
      
        call DB_FreeUnstructuredVariable ( DB_UV_Handle )
      
      end do
  
      if ( len_trim ( VG % Name ) > 0 ) &
        call UGI % Stream % ChangeDirectory ( '../' ) 
    
      end associate 
  
    end do

    nullify ( VariableValue )
    nullify ( ValueArrays )
    nullify ( DB_UV )
    

  end subroutine ReadVariableGroup


end module UnstructuredGridImage_Form
