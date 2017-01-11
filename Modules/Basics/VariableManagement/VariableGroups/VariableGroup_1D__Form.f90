!-- VariableGroup_1D_Form allows bundling of an array of 
!   VariableGroup and stores information about each element of the array.

module VariableGroup_1D__Form

  use Specifiers
  use VariableGroup_Form
  
  implicit none
  private

  type, public :: VariableGroup_1D_Form
    integer ( KDI ) :: &
      nGroups         = 0, &
      nVariablesTotal = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      nVariables
    type ( VariableGroupForm ), dimension ( : ), allocatable :: &
      VariableGroup
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type VariableGroup_1D_Form

contains

  
  subroutine Initialize ( VG_1D, VG )

    class ( VariableGroup_1D_Form ), intent ( inout ) :: &
      VG_1D
    class ( VariableGroupForm ), dimension ( : ), intent ( in ), target :: &
      VG

    integer ( KDI ) :: &
      iG  !-- iGroup

    VG_1D % nGroups = size ( VG )

    allocate ( VG_1D % nVariables ( VG_1D % nGroups ) )
    VG_1D % nVariables &
      = [ ( VG ( iG ) % nVariables, iG = 1, VG_1D % nGroups ) ]           

    VG_1D % nVariablesTotal = sum ( VG_1D % nVariables )
    
    allocate ( VG_1D % VariableGroup ( VG_1D % nGroups ) )
    do iG = 1, VG_1D % nGroups
      call VG_1D % VariableGroup ( iG ) % Initialize ( VG ( iG ) )
    end do

  end subroutine Initialize


  impure elemental subroutine Finalize ( VG_1D )

    type ( VariableGroup_1D_Form ), intent ( inout ) :: &
      VG_1D
    
    if ( allocated ( VG_1D % VariableGroup ) ) &
      deallocate ( VG_1D % VariableGroup )
    
    if ( allocated ( VG_1D % nVariables ) ) deallocate ( VG_1D % nVariables )

  end subroutine Finalize


end module VariableGroup_1D__Form
