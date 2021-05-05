module Eigenspeeds_F_A__Form

  !-- Eigenspeeds_Fast_Atlas_Form

  use Basics
  use Manifolds
  use FieldSets
  use CurrentSet_C__Form
  use CurrentSet_A__Form
  use Eigenspeeds_F_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: Eigenspeeds_F_A_Form
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_F
    generic, public :: &
      Initialize => InitializeAllocate_F
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Eigenspeeds_F_A_Form


contains


  subroutine InitializeAllocate_F &
               ( EA, CSA, FieldOption, NameOption, nFieldsOption )

    class ( Eigenspeeds_F_A_Form ), intent ( inout ) :: &
      EA
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( EA % Type  ==  '' ) &
      EA % Type  =  'an Eigenspeeds_F_A'

    Name  =  'E_' // trim ( CSA % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    EA % CurrentSet_A  =>  CSA

    associate ( nC  =>  CSA % Atlas % nCharts )

    if ( allocated ( EA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( EA % FieldSet_C ( nC ) )
    end if

    call EA % FieldSet_A_Form % Initialize &
           ( CSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = CSA % IGNORABILITY )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Eigenspeeds_F_C_Form :: EA % FieldSet_C ( iC ) % Element ) 
        select type ( EC  =>  EA % FieldSet_C ( iC ) % Element )
        class is ( Eigenspeeds_F_C_Form )

        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )
        call EC % Initialize ( CSC, FieldOption, NameOption, nFieldsOption )
        end select !-- CSC

        end select !-- EC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_F


  subroutine Compute ( EA, iD, TimerLevelOption )

    class ( Eigenspeeds_F_A_Form ), intent ( inout ) :: &
      EA
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( EA % FieldSet_C )
      select type ( EC  =>  EA % FieldSet_C ( iC ) % Element )
      class is ( Eigenspeeds_F_C_Form )
      call EC % Compute ( iD, TimerLevelOption )
      end select !-- EC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( EA )

    type ( Eigenspeeds_F_A_Form ), intent ( inout ) :: &
      EA

    nullify ( EA % CurrentSet_A )

  end subroutine Finalize


end module Eigenspeeds_F_A__Form
