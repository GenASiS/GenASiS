module Eigenspeeds_F_C__Form

  !-- Eigenspeeds_Fast_Chart_Form

  use Basics
  use FieldSets
  use CurrentSet_C__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 2, &
      N_VECTORS_F = 0

  type, public, extends ( FieldSet_C_Form ) :: Eigenspeeds_F_C_Form
    integer ( KDI ) :: &
      iTimer = 0
    integer ( KDI ) :: &
      N_FIELDS_F  = N_FIELDS_F, &
      N_VECTORS_F = N_VECTORS_F
    integer ( KDI ) :: &
      EIGENSPEED_FAST_PLUS_U  = 0, &
      EIGENSPEED_FAST_MINUS_U = 0
    class ( CurrentSet_C_Form ), pointer :: &
      CurrentSet_C => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_E
    generic, public :: &
      Initialize => InitializeAllocate_E
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Eigenspeeds_F_C_Form


contains


  subroutine InitializeAllocate_E &
               ( EC, CSC, FieldOption, NameOption, nFieldsOption )

    class ( Eigenspeeds_F_C_Form ), intent ( inout ) :: &
      EC
    class ( CurrentSet_C_Form ), intent ( in ), target :: &
      CSC
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
 
    integer ( KDI ) :: &
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( EC % Type  ==  '' ) &
      EC % Type  =  'an Eigenspeeds_F_C' 
    
    Name  =  'E_' // trim ( CSC % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    EC % CurrentSet_C  =>  CSC

    associate &
      ( DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
        PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate ) 

    !-- Field indices

    EC % EIGENSPEED_FAST_PLUS_U   =  1
    EC % EIGENSPEED_FAST_MINUS_U  =  2

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      nFields  =  EC % N_FIELDS_F
    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : EC % N_FIELDS_F ) &
      =  [ 'Eigenspeed_F_Plus_U ', &
           'Eigenspeed_F_Minus_U' ]
          
    !-- FieldSet

    call EC % FieldSet_C_Form % Initialize &
           ( CSC % Chart, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CSC % IGNORABILITY )

    end associate !-- DeviceMemory, etc.

  end subroutine InitializeAllocate_E


  subroutine Compute ( EC, iD, TimerLevelOption )

    class ( Eigenspeeds_F_C_Form ), intent ( inout ) :: &
      EC
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    character ( LDL ) :: &
      TimerName
    type ( TimerForm ), pointer :: &
      T 

    associate ( iT  =>  EC % iTimer )
    if ( iT == 0 ) then
      TimerName  =  EC % Name
      if ( present ( TimerLevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, TimerLevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if
    end associate !-- iT

    T  =>  PROGRAM_HEADER % TimerPointer ( EC % iTimer )
    call T % Start ( )

    call Show ( 'Computing ' // trim ( EC % Type ), EC % IGNORABILITY + 4 )
    call Show ( EC % Name, 'Name', EC % IGNORABILITY + 4 )

    associate ( CSC  =>  EC % CurrentSet_C )
    call CSC % ComputeEigenspeeds &
           ( EC, &
             [ EC % EIGENSPEED_FAST_PLUS_U, EC % EIGENSPEED_FAST_MINUS_U ], &
             iD )
    end associate !-- CSC

    call T % Stop ( )

  end subroutine Compute


  impure elemental subroutine Finalize ( EC )

    type ( Eigenspeeds_F_C_Form ), intent ( inout ) :: &
      EC

    nullify ( EC % CurrentSet_C )

  end subroutine Finalize


end module Eigenspeeds_F_C__Form
