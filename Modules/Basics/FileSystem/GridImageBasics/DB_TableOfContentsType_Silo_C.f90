!-- Bindings for the Silo C structure DBtoc and C function DBGetToc

module DB_TableOfContentsType_Silo_C

  use iso_c_binding

  implicit none
  private

  type, public, bind ( c ) :: DB_TableOfContentsType
    type ( c_ptr ) :: &
      Curve
    integer ( c_int ) :: &
      nCurves
    type ( c_ptr ) :: &
      MultiMesh
    integer ( c_int ) :: &
      nMultiMeshes
    type ( c_ptr ) :: &
      MultiMeshAdjacency
    integer ( c_int ) :: &
      nMultiMeshAdjacency 
    type ( c_ptr ) :: &
      MultiVariable 
    integer ( c_int ) :: &
      nMultiVariables
    type ( c_ptr ) :: &
      MultiMaterial
    integer ( c_int ) :: &
      nMultiMaterials
    type ( c_ptr ) :: & 
      MultiMaterialSpecies;
    integer ( c_int ) :: & 
      nMultiMaterialSpecies;
    type ( c_ptr ) :: & 
      CSG_Mesh;
    integer ( c_int ) :: & 
      n_CSG_Meshes;
    type ( c_ptr ) :: & 
      CSG_Variable;
    integer ( c_int ) :: & 
      n_CSG_Variables
    type ( c_ptr ) :: & 
      DefinedVariable;
    integer ( c_int ) :: & 
      nDefinedVariables
    type ( c_ptr ) :: &
      QuadMesh
    integer ( c_int ) :: &
      nQuadMeshes
    type ( c_ptr ) :: &
      QuadVariable
    integer ( c_int ) :: &
      nQuadVariables
    type ( c_ptr ) :: &
      UnstructuredMesh
    integer ( c_int ) :: &
      nUnstructuredMeshes
    type ( c_ptr ) :: &
      UnstructuredVariable
    integer ( c_int ) :: &
      nUnstructuredVariables
    type ( c_ptr ) :: & 
      PointMesh
    integer ( c_int ) :: & 
      nPointMeshes;
    type ( c_ptr ) :: & 
      PointVariable
    integer ( c_int ) :: & 
      nPointVariables;
    type ( c_ptr ) :: & 
      Material;
    integer ( c_int ) :: & 
      nMaterials;
    type ( c_ptr ) :: & 
      MaterialSpecies
    integer ( c_int ) :: & 
      nMaterialSpecies;
    type ( c_ptr ) :: & 
      Variable
    integer ( c_int ) :: & 
      nVariables
    type ( c_ptr ) :: & 
      Object
    integer ( c_int ) :: & 
      nObjects
    type ( c_ptr ) :: & 
      Directory
    integer ( c_int ) :: & 
      nDirectories;
    type ( c_ptr ) :: & 
      Array
    integer ( c_int ) :: & 
      nArrays
    type ( c_ptr ) :: & 
      MRG_Tree
    integer ( c_int ) :: & 
      n_MRG_Trees
    type ( c_ptr ) :: & 
      StorageelMap;
    integer ( c_int ) :: & 
      nStorageelMaps;
    type ( c_ptr ) :: & 
      MRG_Variable;
    integer ( c_int ) :: & 
      n_MRG_Variables;
  end type DB_TableOfContentsType

  public :: &
    DB_Get_TOC

  interface

    type ( c_ptr ) function DB_Get_TOC ( DB_File ) &
                              bind ( c, name = 'DBGetToc' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
    end function DB_Get_TOC
    
  end interface

end module DB_TableOfContentsType_Silo_C
