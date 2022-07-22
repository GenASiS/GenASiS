module Chart_GS_CE__Form
  
  !-- Chart_GridStructured_CentralExcision_Form

  use Basics
  use Chart_GS_C__Form

  implicit none
  private

  type, public, extends ( Chart_GS_C_Form ) :: Chart_GS_CE_Form
    real ( KDR ) :: &
      RadiusExcision
  contains
    procedure, private, pass :: &
      Initialize_GS_CE
    generic, public :: &
      Initialize => Initialize_GS_CE
    procedure, private, pass :: &
      Show_C
    final :: &
      Finalize
  end type Chart_GS_CE_Form


contains


  subroutine Initialize_GS_CE &
               ( C, RadiusMax, RadiusExcision, CommunicatorOption, &
                 NameOption, CoordinateUnitOption, RadialRatioOption, &
                 nGhostLayersOption, nCellsPolarOption, nDimensionsOption )

    class ( Chart_GS_CE_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusExcision
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
      nDimensionsOption

    real ( KDR ) :: &
      RadialRatio

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS_CE'

    C % RadiusExcision  =  RadiusExcision

    RadialRatio  =  1.0_KDR
    if ( present ( RadialRatioOption ) )  &
      RadialRatio  =  RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    call C % Chart_GS_C_Form % Initialize &
           ( RadiusMin = RadiusExcision, &
             RadiusMax = RadiusMax, &
             RadiusScale = RadiusExcision, &
             RadialRatio = RadialRatio, &
             CommunicatorOption = CommunicatorOption, &
             NameOption = NameOption, &
             CoordinateUnitOption = CoordinateUnitOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nCellsPolarOption = nCellsPolarOption, &
             nDimensionsOption = nDimensionsOption )

  end subroutine Initialize_GS_CE


  subroutine Show_C ( C )

    class ( Chart_GS_CE_Form ), intent ( in ) :: &
      C

    call C % Chart_GS_C_Form % Show ( )

    call Show ( 'Chart_GS_CE Proper Parameters' )
    call Show ( C % RadiusExcision, C % CoordinateUnit ( 1 ), 'RadiusExcision' )

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_CE_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize


end module Chart_GS_CE__Form
