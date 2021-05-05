module FluxSet_C__Form

  !-- FluxSet_Chart_Form

  use Basics
  use FieldSets
  use CurrentSet_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_C_Form ) :: FluxSet_C_Form
    integer ( KDI ) :: &
      iTimer = 0
    class ( CurrentSet_C_Form ), pointer :: &
      CurrentSet_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type FluxSet_C_Form


contains


  subroutine InitializeAllocate_F ( FSC, CSC, NameOption )

    class ( FluxSet_C_Form ), intent ( inout ) :: &
      FSC
    class ( CurrentSet_C_Form ), intent ( in ), target :: &
      CSC
    character ( * ), intent ( in ), optional :: &
      NameOption
 
    integer ( KDI ) :: &
      iB, &  !-- iBalanced
      iF     !-- iField
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( FSC % Type  ==  '' ) &
      FSC % Type  =  'a FluxSet_C' 
    
    Name  =  'FS_' // trim ( CSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    FSC % CurrentSet_C  =>  CSC

    associate &
      ( nB  =>  CSC % nBalanced, &
        DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate ) 

    allocate ( Field ( nB ) )
    do iB  =  1,  nB
      iF  =  CSC % iaBalanced ( iB )
      Field ( iB )  =  CSC % Field ( iF )
    end do !-- iS

    call FSC % FieldSet_C_Form % Initialize &
           ( CSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nB, &
             IgnorabilityOption = CSC % IGNORABILITY )

    end associate !-- nB, etc.

  end subroutine InitializeAllocate_F


  subroutine Compute ( FSC, iD, TimerLevelOption )

    class ( FluxSet_C_Form ), intent ( inout ) :: &
      FSC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  FSC % iTimer )
    if ( iT == 0 ) then
      TimerName  =  FSC % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( FSC % iTimer )
    call T % Start ( )

    call Show ( 'Computing ' // trim ( FSC % Type ), FSC % IGNORABILITY + 4 )
    call Show ( FSC % Name, 'Name', FSC % IGNORABILITY + 4 )

    associate ( CSC  =>  FSC % CurrentSet_C )
    call CSC % ComputeFluxes ( FSC, iD )
    end associate !-- CSC

    call T % Stop ( )

  end subroutine Compute


  impure elemental subroutine Finalize ( FSC )

    type ( FluxSet_C_Form ), intent ( inout ) :: &
      FSC

    nullify ( FSC % CurrentSet_C )

  end subroutine Finalize


end module FluxSet_C__Form
