module FluxSet_A__Form

  !-- FluxSet_Atlas_Form

  use Basics
  use FieldSets
  use CurrentSet_C__Form
  use CurrentSet_A__Form
  use FluxSet_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: FluxSet_A_Form
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
  end type FluxSet_A_Form


contains


  subroutine InitializeAllocate_F ( FSA, CSA, NameOption )

    class ( FluxSet_A_Form ), intent ( inout ) :: &
      FSA
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( FSA % Type  ==  '' ) &
      FSA % Type  =  'a FluxSet_A'

    Name  =  'FS_' // trim ( CSA % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    FSA % CurrentSet_A  =>  CSA

    associate ( nC  =>  CSA % Atlas % nCharts )

    if ( allocated ( FSA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( FSA % FieldSet_C ( nC ) )
    end if

    call FSA % FieldSet_A_Form % Initialize &
           ( CSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = CSA % IGNORABILITY )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( FluxSet_C_Form :: FSA % FieldSet_C ( iC ) % Element ) 
        select type ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
          class is ( FluxSet_C_Form )
        select type ( CSC  =>  CSA % FieldSet_C ( iC ) % Element )
          class is ( CurrentSet_C_Form )
  
        call FSC % Initialize ( CSC, NameOption )
  
        end select !-- CSC
        end select !-- FSC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_F


  subroutine Compute ( FSA, iD, TimerLevelOption )

    class ( FluxSet_A_Form ), intent ( inout ) :: &
      FSA
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( FSA % FieldSet_C )
      select type ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
        class is ( FluxSet_C_Form )
      call FSC % Compute ( iD, TimerLevelOption )
      end select !-- FSC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( FSA )

    type ( FluxSet_A_Form ), intent ( inout ) :: &
      FSA

    nullify ( FSA % CurrentSet_A )

  end subroutine Finalize


end module FluxSet_A__Form
