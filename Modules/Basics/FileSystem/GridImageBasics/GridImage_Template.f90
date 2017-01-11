!-- GridImageTemplate is an abstract class whose members will be needed and 
!   inherited by any implementation of GridImage in the form of its extension.

module GridImage_Template

  use VariableManagement
  use Display
  use GridImageStream_Template
        
  implicit none
  private
  
  integer, private, parameter :: &
    MAX_VARIABLE_GROUPS = 16 

  type, public, abstract :: GridImageTemplate
    integer ( KDI ) :: &
      MAX_VARIABLE_GROUPS = MAX_VARIABLE_GROUPS, &
      oValue              = 0, &
      nDimensions         = 0, &
      nTotalCells         = 0, &
      nGhostCells         = 0, &
      nVariableGroups     = 0
    real ( KDR ), dimension ( : ), allocatable :: &
      NodeCoordinate_1, &
      NodeCoordinate_2, &
      NodeCoordinate_3
    character ( LDL ), dimension ( 3 ) :: &
      CoordinateLabel = [ 'X                             ', &
                          'Y                             ', &
                          'Z                             ' ]
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    type ( VariableGroupForm ), dimension ( : ), allocatable :: &
      VariableGroup
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      AddVariableGroup
    procedure, public, pass :: &
      ClearVariableGroups
    !-- FIXME: Removed the following deferred procedures (and their 
    !          interfaces) to avoid multiple definition error during linking
    !          with Cray compiler (CCE). These deferredn procesdures will only
    !          be needed when there is another backend for GridImage
    !          (beside Silo)
    !procedure ( WriteMultiMeshInterface ), public, pass, deferred :: &
    !  WriteMultiMesh
    !procedure ( WriteMultiVariableInterface ), public, pass, deferred :: &
    !  WriteMultiVariable
  end type GridImageTemplate
  
  !-- See FIXME above
  !abstract interface 
  !  
  !  subroutine WriteMultiMeshInterface &
  !               ( GI, Name, TimeOption, CycleNumberOption, HideOption )
  !    use VariableManagement
  !    import GridImageTemplate
  !    class ( GridImageTemplate ), intent ( inout ) :: &
  !      GI
  !    character ( * ), intent ( in ) :: &
  !      Name
  !    type ( MeasuredValueForm ), intent ( in ), optional :: &
  !      TimeOption
  !    integer ( KDI ), intent ( in ), optional :: &
  !      CycleNumberOption
  !    logical ( KDL ), intent ( in ), optional :: &
  !      HideOption
  !  end subroutine WriteMultiMeshInterface
  !  
  !  
  !  subroutine WriteMultiVariableInterface &
  !               ( GI, Name, TimeOption, CycleNumberOption )
  !    use VariableManagement
  !    import GridImageTemplate
  !    class ( GridImageTemplate ), intent ( inout ) :: &
  !      GI
  !    character ( * ), intent ( in ) :: &
  !      Name
  !    type ( MeasuredValueForm ), intent ( in ), optional :: &
  !      TimeOption
  !    integer ( KDI ), intent ( in ), optional :: &
  !      CycleNumberOption
  !  end subroutine WriteMultiVariableInterface
  !
  !end interface
  
  
contains

  
  subroutine InitializeTemplate ( GI )
  
    class ( GridImageTemplate ), intent ( inout ) :: &
      GI

    allocate ( GI % VariableGroup ( GI % MAX_VARIABLE_GROUPS ) )

  end subroutine InitializeTemplate


  subroutine AddVariableGroup ( GI, VG )

    class ( GridImageTemplate ), intent ( inout ), target :: &
      GI
    class ( VariableGroupForm ), intent ( in ) :: &
      VG

    GI % nVariableGroups = GI % nVariableGroups + 1

    if ( GI % nVariableGroups > size ( GI % VariableGroup ) ) then
      call Show &
             ( 'Too many VariableGroups in GridImage', CONSOLE % ERROR )
      return
    end if

    call GI % VariableGroup ( GI % nVariableGroups ) % Initialize ( VG )

  end subroutine AddVariableGroup


  subroutine ClearVariableGroups ( GI )

    class ( GridImageTemplate ), intent ( inout ), target :: &
      GI
    
    if ( allocated ( GI % VariableGroup ) ) deallocate ( GI % VariableGroup )

    allocate ( GI % VariableGroup ( GI % MAX_VARIABLE_GROUPS ) )

    GI % nVariableGroups = 0

  end subroutine ClearVariableGroups


end module GridImage_Template
