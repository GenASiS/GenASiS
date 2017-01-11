!-- Bindings for the Silo C structure DBquadmesh and C functions 
!   DBGetQuadmesh and DBFreeQuadmesh

module DB_QuadMeshType_Silo_C

  use iso_c_binding
  
  implicit none
  private
  
  type, public, bind ( c ) :: DB_QuadMeshType
    integer ( c_int ) :: &
      Identifier, &
      BlockNumber, &
      GroupNumber
    type ( c_ptr ) :: &
      Name
    integer ( c_int ) :: &
      CycleNumber, &
      CoordynateSystem, &
      MajorOrder
    integer ( c_int ), dimension ( 3 ) :: &
      Stride
    integer ( c_int ) :: &
      CoordinateType, &
      FaceType, &
      Planar
    type ( c_ptr ), dimension ( 3 ) :: &
      NodeCoordinate
    integer ( c_int ) :: &
      DataType
    real ( c_float ) :: &
      Time
    real ( c_double ) :: &
      TimeInDoublePrecision
    real ( c_float ), dimension ( 6 ) :: &
      MinimumExtent, &
      MaximumExtent
    type ( c_ptr ), dimension ( 3 ) :: &
      Label, &
      Unit
    integer ( c_int ) :: &
      nDimensions, &
      nSpaces, &
      nNodes
    integer ( c_int ), dimension ( 3 ) :: &
      nNodesPerDimension
    integer ( c_int ) :: &
      Origin
    integer ( c_int ), dimension ( 3 ) :: &
      MinimumIndex, &
      MaximumIndex, &
      BaseIndex, &
      StartIndex, &
      SizeIndex
    integer ( c_int ) :: &
      GUI_HideFlag
    type ( c_ptr ) :: &
      MRG_TreeName
  end type DB_QuadMeshType

  public :: &
    DB_GetQuadMesh, &
    DB_FreeQuadMesh

  interface 
  
    type ( c_ptr ) function DB_GetQuadMesh ( DB_File, MeshName ) &
                              bind ( c, name = 'DBGetQuadmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        MeshName ( * )
    end function DB_GetQuadMesh
  
    subroutine DB_FreeQuadMesh ( DB_QuadMesh ) &
                 bind ( c, name = 'DBFreeQuadmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_QuadMesh      
    end subroutine DB_FreeQuadMesh

  end interface

end module DB_QuadMeshType_Silo_C
