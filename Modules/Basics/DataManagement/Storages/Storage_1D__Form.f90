!-- Storage_1D_Form allows bundling of an array of 
!   Storage and stores information about each element of the array.

module Storage_1D__Form

  use Specifiers
  use Storage_Form
  
  implicit none
  private

  type, public :: Storage_1D_Form
    integer ( KDI ) :: &
      nStorages         = 0, &
      nVariablesTotal = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      nVariables
    type ( StorageForm ), dimension ( : ), allocatable :: &
      Storage
  contains
    procedure, public, pass :: &
      Initialize
    final :: &
      Finalize
  end type Storage_1D_Form

contains

  
  subroutine Initialize ( S_1D, S )

    class ( Storage_1D_Form ), intent ( inout ) :: &
      S_1D
    class ( StorageForm ), dimension ( : ), intent ( in ), target :: &
      S

    integer ( KDI ) :: &
      iStrg  !-- iStorage

    S_1D % nStorages = size ( S )

    allocate ( S_1D % nVariables ( S_1D % nStorages ) )
    S_1D % nVariables &
      = [ ( S ( iStrg ) % nVariables, iStrg = 1, S_1D % nStorages ) ]           

    S_1D % nVariablesTotal = sum ( S_1D % nVariables )
    
    allocate ( S_1D % Storage ( S_1D % nStorages ) )
    do iStrg = 1, S_1D % nStorages
      call S_1D % Storage ( iStrg ) % Initialize ( S ( iStrg ) )
    end do

  end subroutine Initialize


  impure elemental subroutine Finalize ( S_1D )

    type ( Storage_1D_Form ), intent ( inout ) :: &
      S_1D
    
    if ( allocated ( S_1D % Storage ) ) &
      deallocate ( S_1D % Storage )
    
    if ( allocated ( S_1D % nVariables ) ) deallocate ( S_1D % nVariables )

  end subroutine Finalize


end module Storage_1D__Form
