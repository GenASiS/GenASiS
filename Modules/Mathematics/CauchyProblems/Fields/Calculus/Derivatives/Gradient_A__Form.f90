module Gradient_A__Form

  !-- Gradient_Atlas_Form

  use Basics
  use FieldSets
  use Streams
  use Geometries
  use Gradient_C__Form

  implicit none
  private

  type, public, extends ( FieldSet_A_Form ) :: Gradient_A_Form
    class ( FieldSet_A_Form ), pointer :: &
      FieldSet_A  => null ( )
    class ( Geometry_F_A_Form ), pointer :: &
      Geometry_A => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_G
    generic, public :: &
      Initialize => InitializeAllocate_G
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type Gradient_A_Form


contains


  subroutine InitializeAllocate_G ( GA, GyA, FSA, NameOption )

    class ( Gradient_A_Form ), intent ( inout ) :: &
      GA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GyA
    class ( FieldSet_A_Form ), intent ( in ), target :: &
      FSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iC !-- iChart
    logical :: &
      PreviouslyAllocated
    character ( LDL ) :: &
      Name

    if ( GA % Type  ==  '' ) &
      GA % Type  =  'a Gradient_A' 
    
    Name  =  'Grad_' // trim ( FSA % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    GA % FieldSet_A  =>  FSA
    GA % Geometry_A  =>  GyA

    associate ( nC  =>  FSA % Atlas % nCharts )

    if ( allocated ( GA % FieldSet_C ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( GA % FieldSet_C ( nC ) )
    end if

    call GA % FieldSet_A_Form % Initialize &
           ( FSA % Atlas, &
             NameOption = Name, &
             IgnorabilityOption = FSA % IGNORABILITY )

    if ( .not. PreviouslyAllocated ) then
      do iC  =  1,  nC

        allocate ( Gradient_C_Form :: GA % FieldSet_C ( iC ) % Element ) 
        select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
          class is ( Gradient_C_Form )
        select type ( GyC  =>  GyA % FieldSet_C ( iC ) % Element )
          class is ( Geometry_F_C_Form )
        select type ( FSC  =>  FSA % FieldSet_C ( iC ) % Element )
          class is ( FieldSet_C_Form )
  
        call GC % Initialize ( GyC, FSC, NameOption )

        end select !-- FSC
        end select !-- GyC
        end select !-- GC

      end do !-- iC
    end if !-- PreviouslyAllocated

    end associate !-- nC

  end subroutine InitializeAllocate_G


  subroutine Compute ( GA, iD, TimerLevelOption )

    class ( Gradient_A_Form ), intent ( inout ) :: &
      GA
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    integer ( KDI ), intent ( in ), optional :: &
      TimerLevelOption

    integer ( KDI ) :: &
      iC  !-- iChart

    do iC  =  1, size ( GA % FieldSet_C )
      select type ( GC  =>  GA % FieldSet_C ( iC ) % Element )
        class is ( Gradient_C_Form )
      call GC % Compute ( iD, TimerLevelOption )
      end select !-- GC
    end do !-- iC

  end subroutine Compute


  impure elemental subroutine Finalize ( GA )

    type ( Gradient_A_Form ), intent ( inout ) :: &
      GA

    nullify ( GA % Geometry_A )
    nullify ( GA % FieldSet_A )

  end subroutine Finalize


end module Gradient_A__Form
