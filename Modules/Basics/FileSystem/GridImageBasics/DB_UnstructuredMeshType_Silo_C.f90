!-- Bindings for the Silo C structure DBucdmesh, DBzonelist,  and 
!   C functions DBGetUcdmesh and DBFreeUcdmesh

module DB_UnstructuredMeshType_Silo_C

  use iso_c_binding
  
  implicit none 
  private
  
  type, public, bind ( c ) :: DB_UnstructuredMeshType
    integer ( c_int ) :: &
      Identifier, &           !-- Identifier for this object
      BlockNumber, &          !-- Block number for this mesh
      GroupNumber             !-- Block group number for this mesh
    type ( c_ptr ) :: &                     
      Name                    !-- Name associated with mesh
    integer ( c_int ) :: &                  
      CycleNumber, &          !-- Problem cycle number
      CoordinateSystem, &     !-- Coordinate system
      TopologicalDimension    !-- Topological dimension.
    type ( c_ptr ), dimension ( 3 ) :: &
      Unit, &                 !-- Units for variable, e.g, 'mm/ms'
      Label, &                !-- Label associated with each dimension
      NodeCoordinate          !-- Mesh node coordinates
    integer ( c_int ) :: &                  
      DataType                !-- Type of coordinate arrays (double,float)
    real ( c_float ) :: &
      Time                    !-- Problem time
    real ( c_double ) :: &
      TimeInDoublePrecision   !-- Problem time, double data type
    real ( c_float ), dimension ( 6 ) :: &
      MinimumExtent, &
      MaximumExtent
    integer ( c_int ) :: &
      nDimensions, &          !-- Number of computational dimensions
      nNodes, &               !-- Total number of nodes
      Origin                  !-- '0' or '1'
    type ( c_ptr ) :: &
      FaceList, &             !-- Data structure describing mesh faces
      ZoneList, &             !-- Data structure describing mesh zones
      EdgeList, &             !-- Data struct describing mesh edges
      GlobalNodeNumber, &     !-- [nnodes] global node number of each node
      NodeNumber, &           !-- [nnodes] node number of each node
      PolyhedralZoneList      !-- Data structure describing mesh zones
    integer ( c_int ) :: &
      GUI_HideFlag            !--  Flag to hide from post-processor's GUI
    type ( c_ptr ) :: &
      MRG_TreeName            !--  optional name of assoc. mrgtree object
    integer ( c_int ) :: &
      TimeVaryingConnectivityFlag, &
      DisjointMode, &
      GlobalNode_ID_DataType  !--  datatype for global node/zone ids
  end type DB_UnstructuredMeshType
  
  type, public, bind ( c ) :: DB_ZoneListType
    integer ( c_int ) :: &
      nDimensions, &          !-- Number of dimensions (2,3)
      nZones, &               !-- Number of zones in list
      nShapes                 !-- Number of zone shapes
    type ( c_ptr ) :: &
      ShapeCount, &           !-- [nshapes] occurences of each shape
      ShapeSize, &            !-- [nshapes] Number of nodes per shape
      ShapeType, &            !-- [nshapes] Type of shape
      NodeList                !-- Sequent lst of nodes which comprise zones
    integer ( c_int ) :: &
      nNodeLists, &           !-- Number of nodes in nodelist
      Origin, &               !-- '0' or '1'
      MinimumIndex, &         !-- Index of first real zone
      MaximumIndex            !-- Index of last real zone
  end type DB_ZoneListType
  
  interface 
    
    type ( c_ptr ) function DB_GetUnstructuredMesh ( DB_File, MeshName ) &
                              bind ( c, name = 'DBGetUcdmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        MeshName ( * )
    end function DB_GetUnstructuredMesh
    
    subroutine DB_FreeUnstructuredMesh ( DB_UnstructuredMesh ) &
                 bind ( c, name = 'DBFreeUcdmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_UnstructuredMesh
    end subroutine DB_FreeUnstructuredMesh
  
  end interface
  
  public :: &
    DB_GetUnstructuredMesh, &
    DB_FreeUnstructuredMesh
  
end module DB_UnstructuredMeshType_Silo_C
