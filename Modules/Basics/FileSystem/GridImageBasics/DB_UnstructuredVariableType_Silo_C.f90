!-- Bindings for the Silo C structure DBucdvar and C functions 
!   DBGetUcdvar and DBFreeUcdvar

module DB_UnstructuredVariableType_Silo_C

  use iso_c_binding
  
  implicit none
  private
  
  type, public, bind ( c ) :: DB_UnstructuredVariableType
    integer ( c_int ) :: &
      Identifier
    type ( c_ptr ) :: &
      Name
    integer ( c_int ) :: &
      CycleNumber
    type ( c_ptr ) :: &
      Unit, &
      Label
    real ( c_float ) :: &
      Time
    real ( c_double ) :: &
      TimeInDoublePrecision         !-- Problem time, double data type
    integer ( c_int ) :: &
      MeshIdentifier
    type ( c_ptr ) :: &
      Value
    integer ( c_int ) :: &
      DataType, &
      nElements, &
      nValues, &
      nDimensions, &
      Origin, &
      Centering
    type ( c_ptr ) :: &
      MixValue
    integer ( c_int ) :: &
      nMixElements, &
      UseSpeciesMassFractionFlag, &
      ASCII_LabelFlag               !-- Treat variable values as ASCII values
    type ( c_ptr ) :: &
      MeshName
    integer ( c_int ) :: &
      GUI_HideFlag
    type ( c_ptr ) :: &
      RegionPathName
    integer ( c_int ) :: &
      ConservedFlag, &
      ExtensiveFlag
  end type DB_UnstructuredVariableType

  public :: &
    DB_GetUnstructuredVariable, &
    DB_FreeUnstructuredVariable
  
  interface
  
    type ( c_ptr ) function DB_GetUnstructuredVariable &
                              ( DB_File, VariableName ) &
                              bind ( c, name = 'DBGetUcdvar' )
      use iso_c_binding
      implicit none
      type ( c_ptr ), value :: &
        DB_File
      character ( kind = c_char ) :: &
        VariableName ( * )
    end function DB_GetUnstructuredVariable
    
    subroutine DB_FreeUnstructuredVariable ( DB_UnstructuredVariable ) &
                   bind ( c, name = 'DBFreeUcdvar' )
        use iso_c_binding
        implicit none
        type ( c_ptr ), value :: &
          DB_UnstructuredVariable
    end subroutine DB_FreeUnstructuredVariable

  end interface

end module DB_UnstructuredVariableType_Silo_C
