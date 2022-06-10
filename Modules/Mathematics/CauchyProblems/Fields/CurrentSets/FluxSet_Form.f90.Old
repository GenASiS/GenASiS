module FluxSet_Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

  type, public, extends ( FieldSetForm ) :: FluxSetForm
    integer ( KDI ) :: &
      iTimer = 0
    class ( FieldSetForm ), pointer :: &
      FieldSet_CS => null ( )
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_FxS
    generic, public :: &
      Initialize => InitializeAllocate_FxS
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type FluxSetForm


contains


  subroutine InitializeAllocate_FxS ( FS, CS, FS_CS, PrefixOption )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    class ( FieldSetForm ), intent ( in ), target :: &
      FS_CS
    character ( * ), intent ( in ), optional :: &
      PrefixOption
 
    character ( LDL ) :: &
      Name

    if ( FS % Type  ==  '' ) &
      FS % Type  =  'a FluxSet' 
    
    Name  =  'F_' // trim ( FS_CS % Name )
    if ( present ( PrefixOption ) ) &
      Name  =  trim ( PrefixOption ) // '_' // trim ( FS_CS % Name )

    FS % FieldSet_CS  =>  FS_CS
    FS % CurrentSet   =>  CS

    call FS % FieldSetForm % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_FxS


  function Timer ( FS, LevelOption ) result ( T )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  FS % iTimer )

    if ( iT == 0 ) then
      TimerName  =  FS % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  subroutine Compute ( FS, iC, iD )

    class ( FluxSetForm ), intent ( inout ) :: &
      FS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions

    call Show ( 'Computing ' // trim ( FS % Type ), FS % IGNORABILITY + 3 )
    call Show ( FS % Name, 'Name', FS % IGNORABILITY + 3 )

    associate &
      ( CS     =>  FS % CurrentSet, &
        FS_CS  =>  FS % FieldSet_CS )
    call CS % ComputeFluxes ( FS, FS_CS, iC, iD )
    end associate !-- CS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( FS )

    type ( FluxSetForm ), intent ( inout ) :: &
      FS

    nullify ( FS % CurrentSet )
    nullify ( FS % FieldSet_CS )

  end subroutine Finalize


end module FluxSet_Form
