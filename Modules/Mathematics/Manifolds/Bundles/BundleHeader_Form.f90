!-- BundleHeader handles metadata of a Bundle.

module BundleHeader_Form

  use Basics
  use Atlases

  implicit none
  private

  type, public :: BundleHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nFields = 0, &
      nFibersWrite
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
    final :: &
      Finalize
  end type BundleHeaderForm

contains


  subroutine Initialize ( B, Base, NameBase )

    class ( BundleHeaderForm ), intent ( inout ) :: &
      B
    class ( AtlasHeaderForm ), intent ( inout ), target :: &
      Base
    character ( * ), intent ( in )  :: &
      NameBase

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    B % IGNORABILITY = CONSOLE % INFO_1

    if ( B % Type == '' ) &
      B % Type = 'a Bundle' 

    call Split ( B % Type, ' ', TypeWord )
    B % Name = trim ( NameBase ) // '_' // trim ( TypeWord ( 2 ) )

    call Show ( 'Initializing ' // trim ( B % Type ), B % IGNORABILITY )
    call Show ( B % Name, 'Name', B % IGNORABILITY )

    B % Base => Base

    allocate ( B % FiberMaster )
    call B % FiberMaster % Initialize ( B % Name, iDimensionalityOption = 2 )

    B % nFibersWrite = 5
    call PROGRAM_HEADER % GetParameter ( B % nFibersWrite, 'nFibersWrite' )

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
