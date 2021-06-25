module FluxSet_Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: FluxSetForm
    integer ( KDI ) :: &
      iTimer = 0
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FxS
    generic, public :: &
      Initialize => InitializeAllocate_FxS
    procedure, public, pass :: &
      TimerFlux
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type FluxSetForm


contains


  subroutine InitializeAllocate_FxS ( FS, CS, NameOption )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
 
    character ( LDL ) :: &
      Name

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a FluxSet' 
    
    Name  =  'FS_' // trim ( CS % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    FS % CurrentSet  =>  CS

    call FS % FieldSetForm % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_FxS


  function TimerFlux ( FS, TimerLevelOption ) result ( TF )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption
    type ( TimerForm ), pointer :: &
      TF

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  FS % iTimer )

    if ( iT == 0 ) then
      TimerName  =  'Flux_' // trim ( FS % CurrentSet % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    TF  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function TimerFlux


  subroutine Compute ( FS, iD, TimerLevelOption )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    call Show ( 'Computing ' // trim ( FS % Type ), FS % IGNORABILITY + 3 )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY + 3 )

    associate ( CS  =>  FS % CurrentSet )
    call CS % ComputeFluxes ( FS, iD )
    end associate !-- CS

  end subroutine Compute


  impure elemental subroutine Finalize ( FS )

    type ( FluxSetForm ), intent ( inout ) :: &
      FS

    nullify ( FS % CurrentSet )

  end subroutine Finalize


end module FluxSet_Form
