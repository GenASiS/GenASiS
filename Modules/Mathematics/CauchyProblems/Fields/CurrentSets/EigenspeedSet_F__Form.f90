module EigenspeedSet_F__Form

  !-- EigenspeedSet_Fast__Form

  use Basics
  use FieldSets
  use CurrentSet_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_FIELDS_F  = 2, &
      N_VECTORS_F = 0

  type, public, extends ( FieldSetForm ) :: EigenspeedSet_F_Form
    integer ( KDI ) :: &
      iTimer = 0
    integer ( KDI ) :: &
      N_FIELDS_F  = N_FIELDS_F, &
      N_VECTORS_F = N_VECTORS_F
    integer ( KDI ) :: &
      EIGENSPEED_FAST_PLUS_U  = 0, &
      EIGENSPEED_FAST_MINUS_U = 0
    class ( FieldSetForm ), pointer :: &
      FieldSet_CS => null ( )
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_ES
    generic, public :: &
      Initialize => InitializeAllocate_ES
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type EigenspeedSet_F_Form


contains


  subroutine InitializeAllocate_ES &
               ( ES, CS, FS_CS, FieldOption, SuffixOption, nFieldsOption )

    class ( EigenspeedSet_F_Form ), intent ( inout ) :: &
      ES
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    class ( FieldSetForm ), intent ( in ), target :: &
      FS_CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      SuffixOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption
 
    integer ( KDI ) :: &
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( ES % Type  ==  '' ) &
      ES % Type  =  'an EigenspeedSet_F' 
    
    Name  =  trim ( FS_CS % Name ) // '_Egnspd'
    if ( present ( SuffixOption ) ) &
      Name  =  trim ( FS_CS % Name ) // '_Egnspd_' // trim ( SuffixOption )

    ES % FieldSet_CS  =>  FS_CS
    ES % CurrentSet   =>  CS

    !-- Field indices

    ES % EIGENSPEED_FAST_PLUS_U   =  1
    ES % EIGENSPEED_FAST_MINUS_U  =  2

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      nFields  =  ES % N_FIELDS_F
    end if

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( 1 : ES % N_FIELDS_F ) &
      =  [ 'Eigenspeed_F_Plus_U ', &
           'Eigenspeed_F_Minus_U' ]
          
    !-- FieldSet

    call ES % FieldSetForm % Initialize &
           ( CS % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

  end subroutine InitializeAllocate_ES


  function Timer ( ES, Level ) result ( T )

    class ( EigenspeedSet_F_Form ), intent ( inout ) :: &
      ES
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = ES % iTimer, &
               Name = ES % Name, &
               Level = Level )

  end function Timer


  subroutine Compute ( ES, iC, iD )

    class ( EigenspeedSet_F_Form ), intent ( inout ) :: &
      ES
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions

    call Show ( 'Computing ' // trim ( ES % Type ), ES % IGNORABILITY + 3 )
    call Show ( ES % Name, 'Name', ES % IGNORABILITY + 3 )

    associate &
      ( CS     =>  ES % CurrentSet, &
        FS_CS  =>  ES % FieldSet_CS )
    call CS % ComputeEigenspeeds &
           ( ES, FS_CS, &
             [ ES % EIGENSPEED_FAST_PLUS_U, ES % EIGENSPEED_FAST_MINUS_U ], &
             iC, iD )
    end associate !-- CS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( ES )

    type ( EigenspeedSet_F_Form ), intent ( inout ) :: &
      ES

    nullify ( ES % CurrentSet )
    nullify ( ES % FieldSet_CS )

  end subroutine Finalize


end module EigenspeedSet_F__Form
