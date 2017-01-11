!-- Bindings for the Silo C structure DBquadvar and C functions 
!   DBGetQuadvar and DBFreeQuadvar

module DB_QuadVariableType_Silo_C

  use iso_c_binding
  
  implicit none
  private

  type, public, bind ( c ) :: DB_QuadVariableType
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
      nDimensions
    integer ( c_int ), dimension ( 3 ) :: &
      nElementsPerDimension
    integer ( c_int ) :: &
      MajorOrder
    integer ( c_int ), dimension ( 3 ) :: &
      Stride, &
      MinimumIndex, &
      MaximumIndex
    integer ( c_int ) :: &
      Origin
    real ( c_float ) :: &
      Time
    real ( c_double ) :: &
      TimeInDoublePrecision
    real ( c_float ), dimension ( 6 ) :: &
      Alignment
    type ( c_ptr ) :: &
      MixValue
    integer ( c_int ) :: &
      nMixElements, &
      UseSpeciesMassFractionFlag, &
      ASCII_LableFlag
    type ( c_ptr ) :: &
      MeshName
    integer ( c_int ) :: &
      GUI_HideFlag
    type ( c_ptr ) :: &
      RegionPathName
    integer ( c_int ) :: &
      ConservedFlag, &
      ExtensiveFlag, &
      Centering
  end type DB_QuadVariableType

  public :: &
    DB_GetQuadVariable, &
    DB_FreeQuadVariable

  interface 
  
    type ( c_ptr ) function DB_GetQuadVariable ( DB_File, VariableName ) &
                              bind ( c, name = 'DBGetQuadvar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        VariableName ( * )
    end function DB_GetQuadVariable
  
    subroutine DB_FreeQuadVariable ( DB_QV ) &
                 bind ( c, name = 'DBFreeQuadvar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_QV      
    end subroutine DB_FreeQuadVariable

  end interface

end module DB_QuadVariableType_Silo_C
