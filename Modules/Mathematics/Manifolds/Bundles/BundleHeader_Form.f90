!-- BundleHeader handles metadata of a Bundle.

module BundleHeader_Form

  use Basics
  use Atlases

  implicit none
  private

  type, public :: BundleHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFields = 0
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    type ( GridImageStreamForm ), dimension ( : ), allocatable :: &
      GridImageStream
    type ( Atlas_SC_Form ), allocatable :: &
      FiberMaster
    class ( AtlasHeaderForm ), pointer :: &
      Base => null ( )
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      SetBoundaryConditionsFace
    procedure, private, pass :: &
      ShowHeader
    generic, public :: &
      Show => ShowHeader
    final :: &
      Finalize
  end type BundleHeaderForm

contains


  subroutine Initialize ( B, Base, Name )

    class ( BundleHeaderForm ), intent ( inout ) :: &
      B
    class ( AtlasHeaderForm ), intent ( inout ), target :: &
      Base
    character ( * ), intent ( in )  :: &
      Name

    B % IGNORABILITY = CONSOLE % INFO_1

    if ( B % Type == '' ) &
      B % Type = 'a Bundle' 

    B % Name = Name

    call Show ( 'Initializing ' // trim ( B % Type ), B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

    B % Base => Base

    allocate ( B % FiberMaster )
    call B % FiberMaster % Initialize ( B % Name, iDimensionalityOption = 2 )

    allocate ( B % GridImageStream ( ATLAS % MAX_STREAMS ) )

  end subroutine Initialize


  subroutine SetBoundaryConditionsFace ( B, BoundaryCondition, iDimension )

    class ( BundleHeaderForm ), intent ( inout ) :: &
      B
    character ( * ), dimension ( 2 ), intent ( in ) :: &
      BoundaryCondition  !-- [ Inner, Outer ]
    integer ( KDI ), intent ( in ) :: &
      iDimension

    call B % FiberMaster % SetBoundaryConditionsFace &
           ( BoundaryCondition, iDimension )

  end subroutine SetBoundaryConditionsFace


  subroutine ShowHeader ( B )

    class ( BundleHeaderForm ), intent ( inout ) :: &
      B

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( B % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )
    call Show ( B % Base % Name, 'Base', B % IGNORABILITY )

  end subroutine ShowHeader


  impure elemental subroutine Finalize ( B )

    type ( BundleHeaderForm ), intent ( inout ) :: &
      B

    nullify ( B % Base )

    if ( allocated ( B % FiberMaster ) ) &
      deallocate ( B % FiberMaster )
    if ( allocated ( B % GridImageStream ) ) &
      deallocate ( B % GridImageStream )

    if ( B % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( B % Type ), B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

  end subroutine Finalize

  
end module BundleHeader_Form
