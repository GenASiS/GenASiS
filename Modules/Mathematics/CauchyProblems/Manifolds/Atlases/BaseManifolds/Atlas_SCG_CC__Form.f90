module Atlas_SCG_CC__Form

  !-- Atlas_SingleChartGrid_CentralCore_Form

  use Basics
  use Charts
  use Atlas_SCG__Form

  implicit none
  private

  type, public, extends ( Atlas_SCG_Form ) :: Atlas_SCG_CC_Form
    class ( Chart_GS_CC_Form ), pointer :: &
      Chart_GS_CC => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG_CC
    generic, public :: &
      Initialize => Initialize_SCG_CC
    final :: &
      Finalize
  end type Atlas_SCG_CC_Form


contains


  subroutine Initialize_SCG_CC &
               ( A, RadiusMax, RadiusCore, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption )

    class ( Atlas_SCG_CC_Form ), intent ( inout ), target :: &
      A
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusCore
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption

    logical :: &
      PreviouslyAllocated

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SCG_CC'

    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( 1 ) )
    
    if ( allocated ( A % Chart ( 1 ) % Element ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( Chart_GS_CC_Form :: A % Chart ( 1 ) % Element )
    end if

    call A % Atlas_SCG_Form % Initialize ( )

    if ( .not. PreviouslyAllocated ) then

      select type ( C  =>  A % Chart ( 1 ) % Element )
      class is ( Chart_GS_CC_Form )

      call C % Initialize &
             ( RadiusMax, RadiusCore, &
               CommunicatorOption = CommunicatorOption, &
               NameOption = NameOption, &
               CoordinateUnitOption = CoordinateUnitOption, &
               RadialRatioOption = RadialRatioOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsPolarOption = nCellsPolarOption, &
               nEqualOption = nEqualOption )

      end select !--  C

    end if !-- PreviouslyAllocated  

    select type ( C  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_CC_Form )
      A % Chart_GS_CC  =>  C
    end select !-- C
      
  end subroutine Initialize_SCG_CC


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_CC_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS_CC )

  end subroutine Finalize


end module Atlas_SCG_CC__Form
