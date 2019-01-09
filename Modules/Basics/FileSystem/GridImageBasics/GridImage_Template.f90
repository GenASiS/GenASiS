!-- GridImageTemplate is an abstract class whose members will be needed and 
!   inherited by any implementation of GridImage in the form of its extension.

module GridImage_Template

  use Specifiers
  use DataManagement
  use Display
  use GridImageStream_Template
        
  implicit none
  private
  
  integer, private, parameter :: &
    MAX_VARIABLE_GROUPS = 96

  type, public, abstract :: GridImageTemplate
    integer ( KDI ) :: &
      MAX_VARIABLE_GROUPS = MAX_VARIABLE_GROUPS, &
      oValue              = 0, &
      nDimensions         = 0, &
      nTotalCells         = 0, &
      nGhostCells         = 0, &
      nStorages     = 0
    real ( KDR ), dimension ( : ), allocatable :: &
      NodeCoordinate_1, &
      NodeCoordinate_2, &
      NodeCoordinate_3
    character ( LDL ), dimension ( 3 ) :: &
      CoordinateLabel = [ 'X', 'Y', 'Z' ]
    type ( MeasuredValueForm ), dimension ( 3 ) :: &
      CoordinateUnit
    type ( StorageForm ), dimension ( : ), allocatable :: &
      Storage
  contains
    procedure, public, pass :: &
      InitializeTemplate
    procedure, public, pass :: &
      AddStorage
    procedure, public, pass :: &
      ClearStorages
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
  !    use Specifiers
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
  !    use Specifiers
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

    allocate ( GI % Storage ( GI % MAX_VARIABLE_GROUPS ) )

  end subroutine InitializeTemplate


  subroutine AddStorage ( GI, S )

    class ( GridImageTemplate ), intent ( inout ), target :: &
      GI
    class ( StorageForm ), intent ( in ) :: &
      S

    GI % nStorages = GI % nStorages + 1

    if ( GI % nStorages > size ( GI % Storage ) ) then
      call Show &
             ( 'Too many Storages in GridImage', CONSOLE % ERROR )
      return
    end if

    call GI % Storage ( GI % nStorages ) % Initialize ( S )

  end subroutine AddStorage


  subroutine ClearStorages ( GI )

    class ( GridImageTemplate ), intent ( inout ), target :: &
      GI
    
    if ( allocated ( GI % Storage ) ) deallocate ( GI % Storage )

    allocate ( GI % Storage ( GI % MAX_VARIABLE_GROUPS ) )

    GI % nStorages = 0

  end subroutine ClearStorages


end module GridImage_Template
