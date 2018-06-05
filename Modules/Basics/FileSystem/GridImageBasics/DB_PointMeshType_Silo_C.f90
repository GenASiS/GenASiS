!-- Bindings for the Silo C structure DBpointmesh and C functions 
!   DBGetPointmesh and DBFreePointmesh

module DB_PointMeshType_Silo_C

  use iso_c_binding
  
  implicit none
  private
  
  type, public, bind ( c ) :: DB_PointMeshType
    integer ( c_int ) :: &
      Identifier, &
      BlockNumber, &
      StorageNumber
    type ( c_ptr ) :: &
      Name
    integer ( c_int ) :: &
      CycleNumber
    type ( c_ptr ), dimension ( 3 ) :: &
      Unit, &
      Label
    type ( c_ptr ) :: &
      Title
    type ( c_ptr ), dimension ( 3 ) :: &
      NodeCoordinate
    real ( c_float ) :: &
      Time
    real ( c_double ) :: &
      TimeInDoublePrecision
    real ( c_float ), dimension ( 6 ) :: &
      MinimumExtent, &
      MaximumExtent
    integer ( c_int ) :: &
      DataType, &
      nDimensions, &
      nElements, &
      Origin, &
      GUI_HideFlag
    type ( c_ptr ) :: &
      GlobalNodeNumber, &
      MRG_TreeName
    integer ( c_int ) :: &
      GlobalNodeDatatype
  end type DB_PointMeshType
  
  public :: &
    DB_GetPointMesh, &
    DB_FreePointMesh

  interface 
  
    type ( c_ptr ) function DB_GetPointMesh ( DB_File, MeshName ) &
                              bind ( c, name = 'DBGetPointmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        MeshName ( * )
    end function DB_GetPointMesh
  
    subroutine DB_FreePointMesh ( DB_PointMesh ) &
                 bind ( c, name = 'DBFreePointmesh' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_PointMesh      
    end subroutine DB_FreePointMesh

  end interface

end module DB_PointMeshType_Silo_C
