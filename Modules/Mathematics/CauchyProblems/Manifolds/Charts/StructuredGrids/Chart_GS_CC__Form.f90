module Chart_GS_CC__Form
  
  !-- Chart_GridStructured_CentralCore_Form

  use Basics
  use Chart_GS_C__Form

  implicit none
  private

  type, public, extends ( Chart_GS_C_Form ) :: Chart_GS_CC_Form
    integer ( KDI ) :: &
      nCellsCore
    real ( KDR ) :: &
      RadiusCore
  contains
    procedure, private, pass :: &
      Initialize_GS_CC
    generic, public :: &
      Initialize => Initialize_GS_CC
    procedure, private, pass :: &
      Show_C
    final :: &
      Finalize
    procedure, private, pass :: &
      SetCore
  end type Chart_GS_CC_Form


contains


  subroutine Initialize_GS_CC &
               ( C, RadiusMax, RadiusCore, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption, nDimensionsOption )

    class ( Chart_GS_CC_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusCore
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption, &
      nDimensionsOption

    real ( KDR ) :: &
      RadialRatio

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS_CC'

    call C % SetCore ( RadiusCore, nCellsPolarOption )

    RadialRatio  =  2.45_KDR
    if ( present ( RadialRatioOption ) )  &
      RadialRatio  =  RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    call C % Chart_GS_C_Form % Initialize &
           ( RadiusMin = 0.0_KDR, &
             RadiusMax = RadiusMax, &
             RadiusScale = RadiusCore, &
             RadialRatio = RadialRatio, &
             CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             CoordinateUnitOption = CoordinateUnitOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nEqualOption = C % nCellsCore, &
             nDimensionsOption = nDimensionsOption )

  end subroutine Initialize_GS_CC


  subroutine Show_C ( C )

    class ( Chart_GS_CC_Form ), intent ( in ) :: &
      C

    call C % Chart_GS_C_Form % Show ( )

    call Show ( 'Chart_GS_CC Proper Parameters' )
    call Show ( C % RadiusCore, C % CoordinateUnit ( 1 ), 'RadiusCore' )
    call Show ( C % nCellsCore, 'nCellsCore' )

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_CC_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize


  subroutine SetCore ( C, RadiusCore, nCellsPolarOption )

    class ( Chart_GS_CC_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusCore
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    call C % SetPolar ( nCellsPolarOption )

    ! if ( .not. any ( C % nCellsPolar &
    !                    == [ 32, 64, 128, 256, 512, 1024, 2048, 4096 ] ) ) &
    ! then 
    !   call Show ( 'nCellsPolar must be a power of 2 between 32 and 4096', &
    !               CONSOLE % ERROR )
    !   call Show ( 'SetCore', 'subroutine', CONSOLE % ERROR )
    !   call Show ( 'Chart_GS_CC__Form', 'module', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end if

    if ( .not. ( mod ( C % nCellsPolar, 32 )  ==  0 ) ) then 
      call Show ( '32 must evenly divide nCellsPolar', &
                  CONSOLE % ERROR )
      call Show ( C % nCellsPolar, 'nCellsPolar', CONSOLE % ERROR )
      call Show ( 'SetCore', 'subroutine', CONSOLE % ERROR )
      call Show ( 'Chart_GS_CC__Form', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

!-- 40 nCellsCore for 128 nCellsPolar, aspect ratio close to 1
!    C % nCellsCore  =  10 * ( C % nCellsPolar / 32 )

!-- 100 nCellsCore for 128 nCellsPolar, aspect ratio close to 2.5  
    C % nCellsCore  =  25 * ( C % nCellsPolar / 32 )  

    C % RadiusCore  =  RadiusCore

  end subroutine SetCore


end module Chart_GS_CC__Form
