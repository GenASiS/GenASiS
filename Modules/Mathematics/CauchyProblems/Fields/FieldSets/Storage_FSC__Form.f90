module Storage_FSC__Form

  !-- Storage_FieldSetChart_Form

  use Basics
  use Manifolds

  implicit none
  private

  type, public :: Storage_FSC_Form
    integer ( KDI ) :: &
      iTimerUpdateDevice = 0, &
      iTimerUpdateHost   = 0
    logical ( KDL ) :: &
      DeviceMemory, &
      PinnedMemory, &
      DevicesCommunicate
    class ( StorageForm ), allocatable :: &
      Storage
  contains
    procedure, private, pass :: &
      InitializeAllocate
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeAllocate, InitializeClone
    procedure, public, pass :: &
      UpdateDevice => UpdateDevice_SFS
    procedure, public, pass :: &
      UpdateHost => UpdateHost_SFS
    final :: &
      Finalize
  end type Storage_FSC_Form

contains


  subroutine InitializeAllocate &
               ( SFSC, C, Field, Vector, Name, Unit, VectorIndices, nFields, &
                 DeviceMemoryOption, PinnedMemoryOption )

    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC
    class ( Chart_H_Form ), intent ( in ) :: &
      C
    character ( * ), dimension ( : ), intent ( in ) :: &
      Field, &
      Vector
    character ( * ), intent ( in ) :: &
      Name
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ) :: &
      Unit
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ) :: &
      VectorIndices
    integer ( KDI ), intent ( in ) :: &
      nFields
    logical ( KDL ), intent ( in ), optional :: &
      DeviceMemoryOption, &
      PinnedMemoryOption

    SFSC % DeviceMemory  =  .false.
    if ( present ( DeviceMemoryOption ) ) &
      SFSC % DeviceMemory  =  DeviceMemoryOption
    
    SFSC % PinnedMemory  =  .false.
    if ( present ( PinnedMemoryOption ) ) &
      SFSC % PinnedMemory  =  PinnedMemoryOption

    allocate ( SFSC % Storage )
    associate ( S  =>  SFSC % Storage )

    select type ( C )
    class is ( Chart_GS_Form )
      call S % Initialize &
             ( [ C % nCellsLocal, nFields ], &
               VariableOption = Field, &
               VectorOption = Vector, &
               NameOption = Name, &
               ClearOption = .true., &
               PinnedOption = SFSC % PinnedMemory, &
               UnitOption = Unit, &
               VectorIndicesOption = VectorIndices )
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Storage_FSC_Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

      if ( SFSC % DeviceMemory ) &
        call S % AllocateDevice ( )
    end associate !-- S

  end subroutine InitializeAllocate


  subroutine InitializeClone &
               ( SFSC_T, SFSC_S, Vector, Name, VectorIndices, iaSelected )

    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC_T
    class ( Storage_FSC_Form ), intent ( in ) :: &
      SFSC_S
    character ( * ), dimension ( : ), intent ( in ) :: &
      Vector
    character ( * ), intent ( in ) :: &
      Name
    type ( Integer_1D_Form ), dimension ( : ), intent ( in ) :: &
      VectorIndices
    integer ( KDI ), dimension ( : ), intent ( in ) :: &
      iaSelected

    SFSC_T % DeviceMemory        =  SFSC_S % DeviceMemory    
    SFSC_T % PinnedMemory        =  SFSC_S % PinnedMemory
    SFSC_T % DevicesCommunicate  =  SFSC_S % DevicesCommunicate

    allocate ( SFSC_T % Storage )
    associate ( S  =>  SFSC_T % Storage )
    call S % Initialize &
           ( SFSC_S % Storage, &
             VectorOption = Vector, &
             NameOption = Name, &
             VectorIndicesOption = VectorIndices, &
             iaSelectedOption = iaSelected )
    end associate !-- S

  end subroutine InitializeClone


  subroutine UpdateDevice_SFS ( SFSC, TimerLevelOption )

    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  SFSC % iTimerUpdateDevice )
    if ( iT == 0 ) then
      TimerName  =  'UpdateDevice ' // trim ( SFSC % Storage % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( SFSC % iTimerUpdateDevice )

    call T % Start ( )
    call SFSC % Storage % UpdateDevice ( )
    call T % Stop ( )

  end subroutine UpdateDevice_SFS
  

  subroutine UpdateHost_SFS ( SFSC, TimerLevelOption )

    class ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  SFSC % iTimerUpdateHost )
    if ( iT == 0 ) then
      TimerName  =  'UpdateHost ' // trim ( SFSC % Storage % Name )
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( SFSC % iTimerUpdateHost )

    call T % Start ( )
    call SFSC % Storage % UpdateHost ( )
    call T % Stop ( )

  end subroutine UpdateHost_SFS


  subroutine Finalize ( SFSC )

    type ( Storage_FSC_Form ), intent ( inout ) :: &
      SFSC

    if ( allocated ( SFSC % Storage ) ) &
      deallocate ( SFSC % Storage )

  end subroutine Finalize


end module Storage_FSC__Form
