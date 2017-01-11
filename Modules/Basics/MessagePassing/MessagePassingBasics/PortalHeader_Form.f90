!-- PortalHeaderForm handles basic information about an infrastructure for
!   message passing

module PortalHeader_Form

  use VariableManagement
  use Display

  implicit none
  private

  type, public :: PortalHeaderForm
    integer ( KDI ) :: &
      nSources = 0, &
      nTargets = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      nChunksFrom, &
      nChunksTo, &
      Source, &
      Target
  contains
    procedure, public, pass :: &
      Initialize
    procedure, public, pass :: &
      Show => Show_PH
    final :: &
      Finalize
  end type PortalHeaderForm
  
contains

 
  subroutine Initialize ( PH, Source, Target, nChunksFrom, nChunksTo )

    class ( PortalHeaderForm ), intent ( inout ) :: &
      PH
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      Source, &
      Target, &
      nChunksFrom, &
      nChunksTo
    
    PH % nSources = size ( Source )
    PH % nTargets = size ( Target )

    allocate ( PH % nChunksFrom ( PH % nSources ) )
    PH % nChunksFrom = nChunksFrom
    
    allocate ( PH % nChunksTo ( PH % nTargets ) )
    PH % nChunksTo = nChunksTo
    
    allocate ( PH % Source ( PH % nSources ) )
    PH % Source = Source
      
    allocate ( PH % Target ( PH % nTargets ) )
    PH % Target = Target

  end subroutine Initialize


  subroutine Show_PH ( PH, Label, IgnorabilityOption )

    class ( PortalHeaderForm ), intent ( in ) :: &
      PH
    character ( * ), intent ( in ) :: &
      Label
    integer ( KDI ), intent ( in ), optional :: &
      IgnorabilityOption
  
    call Show ( Label, IgnorabilityOption )
    call Show ( PH % Source, 'Source', IgnorabilityOption )
    call Show ( PH % nChunksFrom, 'nChunksFrom', IgnorabilityOption )
    call Show ( PH % Target, 'Target', IgnorabilityOption )
    call Show ( PH % nChunksTo, 'nChunksTo', IgnorabilityOption )

  end subroutine Show_PH


  elemental subroutine Finalize ( PH )

    type ( PortalHeaderForm ), intent ( inout ) :: &
      PH
    
    if ( allocated ( PH % Target ) ) deallocate ( PH % Target )
    if ( allocated ( PH % Source ) ) deallocate ( PH % Source )
    if ( allocated ( PH % nChunksTo ) ) deallocate ( PH % nChunksTo )
    if ( allocated ( PH % nChunksFrom ) ) deallocate ( PH % nChunksFrom )
    
  end subroutine Finalize


end module PortalHeader_Form
