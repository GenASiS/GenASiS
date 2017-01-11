!-- Bindings for the Silo C structure DBmeshvar and C functions 
!   DBGetPointvar and DBFreeMeshvar

module DB_MeshVariableType_Silo_C

  use iso_c_binding
  
  implicit none
  private

  type, public, bind ( c ) :: DB_MeshVariableType
    integer ( c_int ) :: &
      Identifier
    type ( c_ptr ) :: &
      Name, &
      Unit, &
      Label
    integer ( c_int ) :: &
      CycleNumber, &
      MeshIdentifier
    type ( c_ptr ) :: &
      Value
    integer ( c_int ) :: &
      Datatype, &
      nElements, &
      nValues, &
      nSpaces, &
      nDimensions, &
      Origin, &
      Centering
    real ( c_float ) :: &
      Time
    real ( c_double ) :: &
      TimeInDoublePrecision
    real ( c_float ), dimension ( 6 ) :: &
      Alignment
    integer ( c_int ), dimension ( 3 ) :: &
      nElementsPerDimension
    integer ( c_int ) :: &
      MajorOrder
    integer ( c_int ), dimension ( 3 ) :: &
      Stride
    integer ( c_int ), dimension ( 6 ) :: &  
      MinimumIndex, &
      MaximumIndex
    integer ( c_int ) :: &
      ASCII_LableFlag
    type ( c_ptr ) :: &
      MeshName
    integer ( c_int ) :: &
      GUI_HideFlag
    type ( c_ptr ) :: &
      RegionPathName
    integer ( c_int ) :: &
      ConservedFlag, &
      ExtensiveFlag
  end type DB_MeshVariableType

  public :: &
    DB_GetPointVariable, &
    DB_FreeMeshVariable

  interface 
  
    type ( c_ptr ) function DB_GetPointVariable ( DB_File, VariableName ) &
                              bind ( c, name = 'DBGetPointvar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        VariableName ( * )
    end function DB_GetPointVariable
  
    subroutine DB_FreeMeshVariable ( DB_MV ) &
                 bind ( c, name = 'DBFreeMeshvar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_MV      
    end subroutine DB_FreeMeshVariable

  end interface

end module DB_MeshVariableType_Silo_C
