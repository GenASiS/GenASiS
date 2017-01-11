!-- StructuredGridImageSiloForm is a class for handling structured grid 
!   in Silo format, suitable for visualization tool such as VisIt.

module StructuredGridImage_Form
  
  use iso_c_binding
  use VariableManagement
  use Display
  use GridImageBasics
  
  implicit none 
  private
  
  include 'silo_f9x.inc'
  
  type, public, extends ( GridImageSiloTemplate ) &
    :: StructuredGridImageForm 
      integer ( KDI ) :: &
        RECTILINEAR = DB_COLLINEAR, &
        CURVILINEAR = DB_NONCOLLINEAR
      integer ( KDI ) :: &
        MeshType
      integer ( KDI ), dimension ( 3 ) :: &
        oProperCellNodeLow  = 0, &
        oProperCellNodeHigh = 0, &
        nNodes              = 0, &
        nCells              = 0, &
        oValueInner      = 0, &
        oValueOuter      = 0
      integer ( KDI ), dimension ( : ), allocatable :: & 
        PillarHash
  contains
    procedure, public, pass :: &
      SetPillarHash
    procedure, public, pass :: &
      SetGridUnigrid
    procedure, public, pass :: &
      SetGridRefinable
    generic :: &
      SetGrid => SetGridUnigrid, SetGridRefinable
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
  end type StructuredGridImageForm
    
contains

  
  subroutine SetPillarHash ( SGI, BrickHash )
  
    class ( StructuredGridImageForm ), intent ( inout ), target :: &
      SGI 
    integer ( KDI ), dimension ( :, :, : ), intent ( in ) :: &
      BrickHash
    
    allocate ( SGI % PillarHash ( size ( BrickHash ) ) )
    SGI % PillarHash = reshape ( BrickHash, shape = [ size ( BrickHash ) ] )
  
  end subroutine SetPillarHash


  subroutine SetGridUnigrid &
               ( SGI, Directory, Edge, nCells, nGhostInner, &
                 nGhostOuter, oValueInner, oValueOuter, nDimensions, &
                 nProperCells, nGhostCells, CoordinateLabelOption, &
                 CoordinateUnitOption )

    class ( StructuredGridImageForm ), intent ( inout ), target :: &
      SGI 
    character ( * ), intent ( in ) :: &
      Directory    
    type ( Real_1D_Form ), dimension ( : ), intent ( in ) :: &
      Edge
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells, &
      nGhostInner, &
      nGhostOuter, &
      oValueInner, &
      oValueOuter
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nProperCells, &
      nGhostCells
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption

    integer ( KDI ) :: &
      iLB, &  !-- iLowerBound
      iUB     !-- iUpperBound
    
    SGI % oValue      = 0
    SGI % nDimensions = nDimensions
    SGI % nTotalCells = nProperCells + nGhostCells
    SGI % nGhostCells = nGhostCells
    SGI % lDirectory  = len_trim ( Directory )

    SGI % MeshType = SGI % RECTILINEAR
    SGI % MultiMeshType = DB_QUAD_RECT
    SGI % MultiVariableType = DB_QUADVAR

    iLB  =  oValueInner ( 1 )  +  1
    iUB  =  size ( Edge ( 1 ) % Value )  -  oValueOuter ( 1 )
    SGI % nNodes ( 1 ) = size ( Edge ( 1 ) % Value ( iLB : iUB ) )
    allocate ( SGI % NodeCoordinate_1 ( SGI % nNodes ( 1 ) ) )
    SGI % NodeCoordinate_1 = Edge ( 1 ) % Value ( iLB : iUB )

    iLB  =  oValueInner ( 2 )  +  1
    iUB  =  size ( Edge ( 2 ) % Value )  -  oValueOuter ( 2 )
    SGI % nNodes ( 2 ) = size ( Edge ( 2 ) % Value ( iLB : iUB ) )
    allocate ( SGI % NodeCoordinate_2 ( SGI % nNodes ( 2 ) ) )
    SGI % NodeCoordinate_2 = Edge ( 2 ) % Value ( iLB : iUB )
    
    if ( SGI % nDimensions > 2 ) then
      iLB  =  oValueInner ( 3 )  +  1
      iUB  =  size ( Edge ( 3 ) % Value )  -  oValueOuter ( 3 )
      SGI % nNodes ( 3 ) = size ( Edge ( 3 ) % Value ( iLB : iUB ) )
      allocate ( SGI % NodeCoordinate_3 ( SGI % nNodes ( 3 ) ) )
      SGI % NodeCoordinate_3 = Edge ( 3 ) % Value ( iLB : iUB )
    else
      SGI % nNodes ( 3 ) = 0
      allocate ( SGI % NodeCoordinate_3 ( 0 ) )
    end if

    SGI % oProperCellNodeLow  = 0
    SGI % oProperCellNodeHigh = 0
    SGI % oProperCellNodeLow  ( : SGI % nDimensions ) &
      = nGhostInner ( : SGI % nDimensions )
    SGI % oProperCellNodeHigh ( : SGI % nDimensions ) &
      = nGhostOuter ( : SGI % nDimensions )

    SGI % nCells = reshape ( nCells, shape = [ 3 ], pad = [ 1 ] )
    
    SGI % oValueInner = oValueInner
    SGI % oValueOuter = oValueOuter

    SGI % Directory = Directory
    if ( trim ( SGI % Directory ) == '/' ) SGI % Directory = ''
    
    if ( present ( CoordinateUnitOption ) ) &
      SGI % CoordinateUnit ( 1 : size ( CoordinateUnitOption ) ) &
        = CoordinateUnitOption
    
    if ( present ( CoordinateLabelOption ) ) &
      SGI % CoordinateLabel ( 1 : size ( CoordinateLabelOption ) ) &
        = CoordinateLabelOption
        
  end subroutine SetGridUnigrid


  subroutine SetGridRefinable &
               ( SGI, Directory, NodeCoordinate, nCells, nDimensions, &
                 nProperCells, nGhostCells, oValue, &
                 CoordinateUnitOption, CoordinateLabelOption, &
                 ProperCellCoordinateBoundaryOption, MeshTypeOption )
                 
    class ( StructuredGridImageForm ), intent ( inout ), target :: &
      SGI 
    character ( * ), intent ( in ) :: &
      Directory
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      NodeCoordinate
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      nCells
    integer ( KDI ), intent ( in ) :: &
      nDimensions, &
      nProperCells, &
      nGhostCells, &
      oValue
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    real ( KDR ), dimension ( :, : ), intent ( in ), optional :: &
      ProperCellCoordinateBoundaryOption
    integer ( KDI ), intent ( in ), optional :: &
      MeshTypeOption
      
    integer ( KDI ) :: &
      iD, &  !-- iDim
      iN     !-- iNode
    real ( KDR ) :: &
      MinProperCellCoordinate, &
      MaxProperCellCoordinate
    real ( KDR ), dimension ( size ( NodeCoordinate, dim = 1 ) ), target :: &
      ScratchNode_1, &
      ScratchNode_2, &
      ScratchNode_3
    real ( KDR ), dimension ( : ), pointer :: &
      NC, &
      SN
      
    SGI % oValue      = oValue
    SGI % nDimensions = nDimensions
    SGI % nTotalCells = nProperCells + nGhostCells
    SGI % nGhostCells = nGhostCells
    SGI % lDirectory  = len_trim ( Directory )
    
    SGI % MeshType = SGI % RECTILINEAR
    if ( present ( MeshTypeOption ) ) SGI % MeshType = MeshTypeOption
    if ( SGI % MeshType == SGI % RECTILINEAR ) then
      SGI % MultiMeshType = DB_QUAD_RECT
    else if ( SGI % MeshType == SGI % CURVILINEAR ) then
      SGI % MultiMeshType = DB_QUAD_CURV
    end if
    SGI % MultiVariableType = DB_QUADVAR
    
    SGI % nNodes = size ( NodeCoordinate, dim = 2 ) 
!-- FIXME: NAG 5.3.1 should support sourced allocation
!    allocate ( SGI % NodeCoordinate_1, source = NodeCoordinate ( :, 1 ) )
!    allocate ( SGI % NodeCoordinate_2, source = NodeCoordinate ( :, 2 ) )
!    allocate ( SGI % NodeCoordinate_3, source = NodeCoordinate ( :, 3 ) )
    allocate ( SGI % NodeCoordinate_1 ( size ( NodeCoordinate ( :, 1 ) ) ) )
    allocate ( SGI % NodeCoordinate_2 ( size ( NodeCoordinate ( :, 2 ) ) ) )
    allocate ( SGI % NodeCoordinate_3 ( size ( NodeCoordinate ( :, 3 ) ) ) )
    SGI % NodeCoordinate_1 = NodeCoordinate ( :, 1 )
    SGI % NodeCoordinate_2 = NodeCoordinate ( :, 2 )
    SGI % NodeCoordinate_3 = NodeCoordinate ( :, 3 )

    if ( SGI % MeshType == SGI % RECTILINEAR ) then
      
      !-- Remove redundant nodes for rectilinear quadmesh
    
      call Sort ( SGI % NodeCoordinate_1 ( : ) )
      call Sort ( SGI % NodeCoordinate_2 ( : ) )
      call Sort ( SGI % NodeCoordinate_3 ( : ) )
      
      ScratchNode_1 ( 1 ) = SGI % NodeCoordinate_1 ( 1 )
      ScratchNode_2 ( 1 ) = SGI % NodeCoordinate_2 ( 1 )
      ScratchNode_3 ( 1 ) = SGI % NodeCoordinate_3 ( 1 )

      SGI % nNodes = 1
      
      associate ( nN => SGI % nNodes ( : ) )
        
        do iD = 1, nDimensions
          
          select case ( iD )
            case ( 1 )
              NC => SGI % NodeCoordinate_1
              SN => ScratchNode_1
            case ( 2 )
              NC => SGI % NodeCoordinate_2
              SN => ScratchNode_2
            case ( 3 )
              NC => SGI % NodeCoordinate_3
              SN => ScratchNode_3
          end select 
          
          do iN = 2, size ( NodeCoordinate, dim = 1 )
            if ( NC ( iN ) > SN ( nN ( iD ) ) ) then
              nN ( iD ) = nN ( iD ) + 1
              SN ( nN ( iD ) ) = NC ( iN )
            end if
          end do
        
        end do
        
        deallocate ( SGI % NodeCoordinate_1 )
!-- FIXME: Intel 12.1.2 did not like this
!        allocate &
!          ( SGI % NodeCoordinate_1, source = ScratchNode_1 ( 1 : nN ( 1 ) ) )
        allocate ( SGI % NodeCoordinate_1 ( nN ( 1 ) ) )
        SGI % NodeCoordinate_1 = ScratchNode_1 ( 1 : nN ( 1 ) )

        deallocate ( SGI % NodeCoordinate_2 )
!-- FIXME: Intel 12.1.2 did not like this
!        allocate &
!          ( SGI % NodeCoordinate_2, source = ScratchNode_2 ( 1 : nN ( 2 ) ) )
        allocate ( SGI % NodeCoordinate_2 ( nN ( 2 ) ) )
        SGI % NodeCoordinate_2 = ScratchNode_2 ( 1 : nN ( 2 ) )
        
        deallocate ( SGI % NodeCoordinate_3 )
!-- FIXME: Intel 12.1.2 did not like this
!        allocate &
!          ( SGI % NodeCoordinate_3, source = ScratchNode_3 ( 1 : nN ( 3 ) ) )
        allocate ( SGI % NodeCoordinate_3 ( nN ( 3 ) ) )
        SGI % NodeCoordinate_3 = ScratchNode_3 ( 1 : nN ( 3 ) )
        
      end associate
      
      do iD = 1, nDimensions
        
        select case ( iD )
          case ( 1 )
            NC  => SGI % NodeCoordinate_1
          case ( 2 )
            NC  => SGI % NodeCoordinate_2
          case ( 3 )
            NC  => SGI % NodeCoordinate_3
        end select
        
        if ( present ( ProperCellCoordinateBoundaryOption ) ) then
          MinProperCellCoordinate &
            = minval ( ProperCellCoordinateBoundaryOption ( :, iD ) )
          MaxProperCellCoordinate &
            = maxval ( ProperCellCoordinateBoundaryOption ( :, iD ) )
        else 
          MinProperCellCoordinate = minval ( NC )
          MaxProperCellCoordinate = maxval ( NC ) 
        end if
        
        do iN = 1, SGI % nNodes ( iD )
          if ( NC ( iN ) == MinProperCellCoordinate ) &
          then
            SGI % oProperCellNodeLow ( iD ) = iN - 1
            exit
          end if
        end do
        
        do iN = SGI % nNodes ( iD ), 1, -1 
          if ( NC ( iN ) == MaxProperCellCoordinate ) &
          then
            SGI % oProperCellNodeHigh ( iD ) = SGI % nNodes ( iD ) - iN
            exit
          end if
        end do
        
      end do
      
    end if
    
    SGI % nCells = reshape ( nCells, shape = [ 3 ], pad = [ 1 ] )
    
    SGI % Directory = Directory
    if ( trim ( SGI % Directory ) == '/' ) SGI % Directory = ''
    
    if ( present ( CoordinateUnitOption ) ) &
      SGI % CoordinateUnit ( 1 : size ( CoordinateUnitOption ) ) &
        = CoordinateUnitOption
    
    if ( present ( CoordinateLabelOption ) ) &
      SGI % CoordinateLabel ( 1 : size ( CoordinateLabelOption ) ) &
        = CoordinateLabelOption
        
    nullify ( SN )
    nullify ( NC )

  end subroutine SetGridRefinable
  
  
  subroutine SetReadAttributes ( GI, Directory, oValue )
  
    class ( StructuredGridImageForm ), intent ( inout ) :: &
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
  
    class ( StructuredGridImageForm ), intent ( inout ) :: &
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
    
    call GI % WriteVariableGroup ( TimeOption, CycleNumberOption ) 
    
    call GI % Stream % ChangeDirectory ( WorkingDirectory ) 
  
  end subroutine Write
  
  
  subroutine Read ( GI, TimeOption, CycleNumberOption ) 
    
    class ( StructuredGridImageForm ), intent ( inout ) :: &
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
               ( ContentTypeOption = 'StructuredGridVariable' )
        if ( allocated ( VariableName ) ) deallocate ( VariableName )
!        allocate ( VariableName, source = GI % Stream % ContentList )
        allocate ( VariableName ( size ( GI % Stream % ContentList ) ) )
        VariableName = GI % Stream % ContentList 
        nVariables = size ( GI % Stream % ContentList )
        if ( nVariables == 0 ) then
          GI % nVariableGroups = 0
        else
          call GI % VariableGroup ( iG ) % Initialize &
                 ( [ product ( GI % nNodes ( 1 : GI % nDimensions ) - 1 ), &
                     nVariables ], &
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

  
  subroutine ClearGrid ( SGI )

    class ( StructuredGridImageForm ), intent ( inout ) :: & 
      SGI 
    
    if ( allocated ( SGI % NodeCoordinate_3 ) ) &
      deallocate ( SGI % NodeCoordinate_3 )
    if ( allocated ( SGI % NodeCoordinate_2 ) ) &
      deallocate ( SGI % NodeCoordinate_2 )
    if ( allocated ( SGI % NodeCoordinate_1 ) ) &
      deallocate ( SGI % NodeCoordinate_1 )
      
    if ( allocated ( SGI % PillarHash ) ) &
      deallocate ( SGI % PillarHash ) 
      
  end subroutine ClearGrid


  impure elemental subroutine Finalize ( SGI )
    
    type ( StructuredGridImageForm ), intent ( inout ) :: & 
      SGI 
    
    nullify ( SGI % Stream )

    if ( allocated ( SGI % VariableGroup ) ) &
      deallocate ( SGI % VariableGroup )

    call SGI % ClearGrid ( )

  end subroutine Finalize 
  
  
  subroutine WriteMesh ( SGI, TimeOption, CycleNumberOption )
    
    class ( StructuredGridImageForm ), intent ( inout ) :: &
      SGI
    type ( MeasuredValueForm ), intent ( in ) , optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
    
    integer ( KDI ) :: &
      iD, &     !-- iDimension
      nSiloOptions, &
      SiloOptionList, &
      CycleNumber, &
      Error
      
    nSiloOptions = 5
    if ( present ( TimeOption ) ) &
      nSiloOptions = nSiloOptions + 1
    if ( present ( CycleNumberOption ) ) &
      nSiloOptions = nSiloOptions + 1
    do iD = 1, 3
      if ( len_trim ( SGI % CoordinateUnit ( iD ) % Label ) > 0 ) &
        nSiloOptions = nSiloOptions + 1
    end do
      
    Error = DBMKOPTLIST ( nSiloOptions, SiloOptionList )

!-- FIXME: Tweak for NAG 5.3.1    
!    Error = DBADDIOPT &
!              ( SiloOptionList, DBOPT_LO_OFFSET, SGI % oProperCellNodeLow )
!    Error = DBADDIOPT &
!              ( SiloOptionList, DBOPT_HI_OFFSET, SGI % oProperCellNodeHigh )
    Error = DBADDIOPT &
              ( SiloOptionList, DBOPT_LO_OFFSET, &
                SGI % oProperCellNodeLow ( 1 ) )
    Error = DBADDIOPT &
              ( SiloOptionList, DBOPT_HI_OFFSET, &
                SGI % oProperCellNodeHigh ( 1 ) )
    
    Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_XLABEL, &
                  trim ( SGI % CoordinateLabel ( 1 ) ), &
                  len_trim ( SGI % CoordinateLabel ( 1 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_YLABEL, &
                trim ( SGI % CoordinateLabel ( 2 ) ), &
                len_trim ( SGI % CoordinateLabel ( 2 ) ) )
    Error = DBADDCOPT &
              ( SiloOptionList, DBOPT_ZLABEL, &
                trim ( SGI % CoordinateLabel ( 3 ) ), &
                len_trim ( SGI % CoordinateLabel ( 3 ) ) )
    
    if ( len_trim ( SGI % CoordinateUnit ( 1 ) % Label ) > 0 ) &
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_XUNITS, &
                  trim ( SGI % CoordinateUnit ( 1 ) % Label ), &
                  len_trim ( SGI % CoordinateUnit ( 1 ) % Label ) )
    if ( len_trim ( SGI % CoordinateUnit ( 2 ) % Label ) > 0 ) &   
      Error = DBADDCOPT & 
                ( SiloOptionList, DBOPT_YUNITS, &
                  trim ( SGI % CoordinateUnit ( 2 ) % Label ), &
                  len_trim ( SGI % CoordinateUnit ( 2 ) % Label ) )
    if ( len_trim ( SGI % CoordinateUnit ( 3 ) % Label ) > 0 ) &   
      Error = DBADDCOPT &
                ( SiloOptionList, DBOPT_ZUNITS, &
                  trim ( SGI % CoordinateUnit ( 3 ) % Label ), &
                  len_trim ( SGI % CoordinateUnit ( 3 ) % Label ) )

    if ( present ( CycleNumberOption ) ) &
      Error = DBADDIOPT ( SiloOptionList, DBOPT_CYCLE, CycleNumberOption )
    if ( present ( TimeOption ) ) &
      Error = DBADDDOPT ( SiloOptionList, DBOPT_DTIME, TimeOption % Number )
    
    Error = DBPUTQM &
              ( SGI % Stream % MeshBlockHandle, 'Mesh', 4, &
                SGI % CoordinateLabel ( 1 ), &
                len_trim ( SGI % CoordinateLabel ( 1 ) ), &
                SGI % CoordinateLabel ( 2 ), &
                len_trim ( SGI % CoordinateLabel ( 2 ) ), &
                SGI % CoordinateLabel ( 3 ), &
                len_trim ( SGI % CoordinateLabel ( 3 ) ), &
                SGI % NodeCoordinate_1 / SGI % CoordinateUnit ( 1 ) % Number, &
                SGI % NodeCoordinate_2 / SGI % CoordinateUnit ( 2 ) % Number, &
                SGI % NodeCoordinate_3 / SGI % CoordinateUnit ( 3 ) % Number, &
                SGI % nNodes, SGI % nDimensions, DB_DOUBLE, SGI % MeshType, &
                SiloOptionList, Error )
    
    Error = DBFREEOPTLIST ( SiloOptionList )
    
    call SGI % WriteMultiMesh ( 'Mesh', TimeOption, CycleNumberOption )
  
  end subroutine WriteMesh
  
  
  subroutine ReadMesh ( SGI )
    
    class ( StructuredGridImageForm ), intent ( inout ) :: &
      SGI
    
    real ( c_double ), dimension ( : ), pointer :: &
      NC_1, &
      NC_2, &
      NC_3
    type ( c_ptr ) :: &
      DB_File, &
      DB_QM_Handle
    type ( DB_QuadMeshType ), pointer :: &
      DB_QM
    
    DB_File &
      = SGI % Stream % AccessSiloPointer ( SGI % Stream % MeshBlockHandle )
    DB_QM_Handle = DB_GetQuadMesh ( DB_File, c_char_'Mesh' // c_null_char ) 
    call c_f_pointer ( DB_QM_Handle, DB_QM )
    
    SGI % nDimensions = DB_QM % nDimensions
    SGI % MeshType = DB_QM % CoordinateType
    if ( SGI % MeshType == SGI % RECTILINEAR ) then
      SGI % MultiMeshType = DB_QUAD_RECT
    else if ( SGI % MeshType == SGI % CURVILINEAR ) then
      SGI % MultiMeshType = DB_QUAD_CURV
    end if
    SGI % MultiVariableType = DB_QUADVAR
    
    associate ( nDims => SGI % nDimensions )
    
    SGI % nNodes = 1
    SGI % nNodes ( 1 : nDims ) = DB_QM % nNodesPerDimension ( 1 : nDims )
    
    SGI % oProperCellNodeLow ( 1 : nDims ) = DB_QM % MinimumIndex ( 1 : nDims ) 
    SGI % oProperCellNodeHigh ( 1 : nDims ) &
      = SGI % nNodes ( 1 : nDims ) - 1 - DB_QM % MaximumIndex ( 1 : nDims ) 
    
    !-- nCells, nTotalCells, and nGhostCells are set in ReadVariableGroup ( )
    
    end associate
    
    if ( DB_QM % nDimensions >= 1 ) then
      call c_f_pointer &
             ( DB_QM % NodeCoordinate ( 1 ), NC_1, &
               [ DB_QM % nNodesPerDimension ( 1 ) ] )
      if ( allocated ( SGI % NodeCoordinate_1 ) ) &
        deallocate ( SGI % NodeCoordinate_1 )
      allocate &
        ( SGI % NodeCoordinate_1 ( DB_QM % nNodesPerDimension ( 1 ) ), &
          Source = NC_1 ( 1 : DB_QM % nNodesPerDimension ( 1 ) ) )
    else 
      allocate ( SGI % NodeCoordinate_1 ( 0 ) )
    end if
    
    if ( DB_QM % nDimensions >= 2 ) then
      call c_f_pointer &
             ( DB_QM % NodeCoordinate ( 2 ), NC_2, &
               [ DB_QM % nNodesPerDimension ( 2 ) ] )
      if ( allocated ( SGI % NodeCoordinate_2 ) ) &
        deallocate ( SGI % NodeCoordinate_2 )
      allocate &
        ( SGI % NodeCoordinate_2 ( DB_QM % nNodesPerDimension ( 2 ) ), &
          Source = NC_2 ( 1 : DB_QM % nNodesPerDimension ( 2 ) ) )
    else 
      allocate ( SGI % NodeCoordinate_2 ( 0 ) )
    end if
    
    if ( DB_QM % nDimensions == 3 ) then
      call c_f_pointer &
             ( DB_QM % NodeCoordinate ( 3 ), NC_3, &
               [ DB_QM % nNodesPerDimension ( 3 ) ] )
      if ( allocated ( SGI % NodeCoordinate_3 ) ) &
        deallocate ( SGI % NodeCoordinate_3 )
      allocate &
        ( SGI % NodeCoordinate_3 ( DB_QM % nNodesPerDimension ( 3 ) ), &
          Source = NC_3 ( 1 : DB_QM % nNodesPerDimension ( 3 ) ) )
    else 
      allocate ( SGI % NodeCoordinate_3 ( 0 ) )
    end if
    
    call DB_FreeQuadMesh ( DB_QM_Handle )
    
    !-- FIXME: Here we make the assumption that the CoordinateUnit for reading 
    !          is the same as the ones used to write. A better way would be to
    !          read the unit directly from Silo attribute (but that's harder
    !          to do esp. with bugs on C compatibility for CCE compiler, hence
    !          more wrappers would had to be written)
    SGI % NodeCoordinate_1 &
      = SGI % NodeCoordinate_1 * SGI % CoordinateUnit ( 1 ) % Number
    SGI % NodeCoordinate_2 &
      = SGI % NodeCoordinate_2 * SGI % CoordinateUnit ( 2 ) % Number
    SGI % NodeCoordinate_3 &
      = SGI % NodeCoordinate_3 * SGI % CoordinateUnit ( 3 ) % Number
    
    DB_QM_Handle = c_null_ptr
    DB_File      = c_null_ptr
    
    nullify ( NC_3, NC_2, NC_1 )
    nullify ( DB_QM )

  end subroutine ReadMesh
  
  
  subroutine WriteVariableGroup ( SGI, TimeOption, CycleNumberOption ) 

    class ( StructuredGridImageForm ), intent ( inout ) :: &
      SGI
    type ( MeasuredValueForm ), intent ( in ) , optional :: &
      TimeOption
    integer ( KDI ), intent ( in ), optional :: &
      CycleNumberOption
  
    integer ( KDI ) :: &
      iV, &      !-- iValue
      iVrbl, &   !-- iVariable
      iS, &      !-- iSelected
      iG, &      !-- iGroup
      nSiloOptions, &
      Centering, &
      SiloOptionList, &
      Error
    real ( KDR ), dimension ( : ), allocatable, target :: &
      Value
    real ( KDR ), dimension ( :, :, : ), pointer :: &
      Value_PGF, &  !-- Value_ProperGhostFull
      Value_PG      !-- Value_ProperGhost
    character ( LDF ) :: &
      MeshDirectory
  
    associate &
      ( nC  => SGI % nCells, &
        oVI => SGI % oValueInner, &
        oVO => SGI % oValueOuter )

    if ( allocated ( SGI % PillarHash ) ) then
      allocate ( Value ( size ( SGI % PillarHash ) ) )
    else
      allocate ( Value ( product ( nC ) ) )
    end if

    if ( all ( SGI % nNodes == SGI % nCells ) ) then
      Centering = DB_NODECENT
    else
      Centering = DB_ZONECENT
    end if
  
    SiloOptionList = DB_F77NULL
  
    MeshDirectory = SGI % Stream % CurrentDirectory
  
    do iG = 1, SGI % nVariableGroups
    
      associate ( VG => SGI % VariableGroup ( iG ) )
    
      call Show ( 'Writing a VariableGroup (structured)', CONSOLE % INFO_5 )
      call Show ( iG, 'iGroup', CONSOLE % INFO_5 )
      call Show ( VG % Name, 'Name', CONSOLE % INFO_5 )

      call SGI % Stream % MakeDirectory ( VG % Name ) 
    
      do iS = 1, VG % nVariables
        
        iVrbl = VG % iaSelected ( iS )
        
        call Show ( 'Writing a Variable (structured)', CONSOLE % INFO_6 )
        call Show ( iS, 'iSelected', CONSOLE % INFO_6 )
        call Show ( VG % Variable ( iVrbl ), 'Name', CONSOLE % INFO_6 )

        if ( allocated ( SGI % PillarHash ) ) then
          do iV = 1, size ( SGI % PillarHash )
            Value ( iV ) &
              = VG % Value ( SGI % oValue + SGI % PillarHash ( iV ), iVrbl )
          end do
        else
          Value_PGF ( -oVI ( 1 ) + 1 : nC ( 1 ) + oVO ( 1 ), &
                      -oVI ( 2 ) + 1 : nC ( 2 ) + oVO ( 2 ), &
                      -oVI ( 3 ) + 1 : nC ( 3 ) + oVO ( 3 ) ) &
            => VG % Value ( SGI % oValue + 1 :, iVrbl )
          Value_PG ( 1 : nC ( 1 ), 1 : nC ( 2 ), 1 : nC ( 3 ) ) &
            => Value 
          call Copy ( Value_PGF ( 1 : nC ( 1 ), 1 : nC ( 2 ), 1 : nC ( 3 ) ), &
                      Value_PG ( 1 : nC ( 1 ), 1 : nC ( 2 ), 1 : nC ( 3 ) ) )
        end if
        
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

        Error = DBPUTQV1 &
                  ( SGI % Stream % MeshBlockHandle, &
                    trim ( VG % Variable ( iVrbl ) ), &
                    VG % lVariable ( iVrbl ), &
                    trim ( MeshDirectory ) // 'Mesh', &
                    len_trim ( MeshDirectory ) + 4, &
                    Value / VG % Unit ( iVrbl ) % Number, SGI % nCells, &
                    SGI % nDimensions, DB_F77NULL, 0, DB_DOUBLE, Centering, &
                    SiloOptionList, Error )

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
          
        call SGI % WriteMultiVariable &
               ( VG % Variable ( iVrbl ), TimeOption, CycleNumberOption )
        
      end do
      
      call SGI % WriteVectorVariable ( VG )
      
      call SGI % Stream % ChangeDirectory ( '../' )
      
      end associate 
  
    end do
  
    end associate !-- nC, etc.

    nullify ( Value_PG, Value_PGF )

  end subroutine WriteVariableGroup 
  
  
  subroutine ReadVariableGroup ( SGI ) 

    class ( StructuredGridImageForm ), intent ( inout ) :: &
      SGI
    
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
    logical ( KDL ) :: &
      SetSizes
    character ( LDF ) :: &
      MeshDirectory
    character ( kind = c_char, len = LDL ) :: &
      VariableName
    type ( c_ptr ) :: &
      DB_File, &
      DB_QV_Handle
    type ( c_ptr ), dimension ( : ), pointer :: &
      ValueArrays
    type ( DB_QuadVariableType ), pointer :: &
      DB_QV

    MeshDirectory = SGI % Stream % CurrentDirectory
  
    do iG = 1, SGI % nVariableGroups
    
      associate ( VG => SGI % VariableGroup ( iG ), &
                  nDims => SGI % nDimensions )
  
      call Show ( 'Reading a VariableGroup', CONSOLE % INFO_5 )
      call Show ( iG, 'iGroup', CONSOLE % INFO_5 )
      call Show ( VG % Name, 'Name', CONSOLE % INFO_5 )

      if ( len_trim ( VG % Name ) > 0 ) &
        call SGI % Stream % ChangeDirectory ( VG  % Name ) 
    
      DB_File &
        = SGI % Stream % AccessSiloPointer ( SGI % Stream % MeshBlockHandle )
      
      SetSizes = .true. 
      
      do iS = 1, VG % nVariables
      
        iVrbl = VG % iaSelected ( iS )
        VariableName = trim ( VG % Variable ( iVrbl ) ) // c_null_char

        DB_QV_Handle = DB_GetQuadVariable ( DB_File, VariableName )
        call c_f_pointer ( DB_QV_Handle, DB_QV )
      
        call c_f_pointer ( DB_QV % Value, ValueArrays, [ DB_QV % nValues ] )
        
        if ( SetSizes ) then
          SGI % nCells ( 1 : nDims ) &
            = DB_QV % nElementsPerDimension ( 1 : nDims ) 
          SGI % nTotalCells = product ( SGI % nCells )
          if ( DB_QV % Centering == DB_NODECENT ) then
            nProperCells &
              = product ( DB_QV % MaximumIndex ( 1 : nDims ) &
                          - DB_QV % MinimumIndex ( 1 : nDims ) + 1 )
          elseif ( DB_QV % Centering == DB_ZONECENT ) then
            nProperCells &
              = product ( DB_QV % MaximumIndex ( 1 : nDims ) &
                          - DB_QV % MinimumIndex ( 1 : nDims ) ) 
          end if
          SGI % nGhostCells = SGI % nTotalCells - nProperCells 
          SetSizes = .false.
        end if
      
        oV = 0
        do iA = 1, DB_QV % nValues
          !-- FIXME: An assumption is made that the unit used to write
          !          and read are the same. A better way would be to read
          !          the unit directly from Silo file.
          call c_f_pointer &
                 ( ValueArrays ( iA ), VariableValue, [ DB_QV % nElements ] )
          if ( allocated ( SGI % PillarHash ) ) then
            do iV = 1, DB_QV % nElements
              VG % Value &
                ( SGI % oValue + SGI % PillarHash ( oV + iV ), iVrbl ) &
                = VariableValue ( oV + iV ) * VG % Unit ( iVrbl ) % Number
            end do
          else
            call Copy &
                   ( VariableValue ( oV + 1 : ) &
                       * VG % Unit ( iVrbl ) % Number, &
                     VG % Value ( SGI % oValue + oV + 1 &
                                  : SGI % oValue + oV + DB_QV % nElements, &
                                  iVrbl ) )                                
          end if
        end do
      
        call DB_FreeQuadVariable ( DB_QV_Handle )
      
      end do
  
      if ( len_trim ( VG % Name ) > 0 ) &
        call SGI % Stream % ChangeDirectory ( '../' ) 
    
      end associate 
  
    end do
    
    DB_File      = c_null_ptr
    DB_QV_Handle = c_null_ptr
    
    nullify ( VariableValue )
    nullify ( ValueArrays )
    nullify ( DB_QV )

  end subroutine ReadVariableGroup


end module StructuredGridImage_Form
