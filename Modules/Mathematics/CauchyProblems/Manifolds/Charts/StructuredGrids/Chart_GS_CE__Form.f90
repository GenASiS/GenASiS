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
                 CoordinateUnitOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nEqualOption  )

    class ( Chart_GS_CE_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusExcision
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS_CE'

    C % RadiusExcision  =  RadiusExcision

    call C % Chart_GS_C_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             RadiusMin = RadiusExcision, &
             RadiusMax = RadiusMax, &
             RadiusScale = RadiusExcision, &
             CoordinateUnitOption = CoordinateUnitOption, &
             RadialRatioOption = RadialRatioOption, &
             nGhostLayersOption = nGhostLayersOption, &
             nCellsPolarOption = nCellsPolarOption )

  end subroutine Initialize_GS_CE


  subroutine Show_C ( C )

    class ( Chart_GS_CE_Form ), intent ( in ) :: &
      C

    call C % Chart_GS_C_Form % Show ( )

    call Show ( 'Chart_GS_CE parameters' )
    call Show ( C % RadiusExcision, C % CoordinateUnit ( 1 ), 'RadiusExcision' )

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_CE_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize


end module Chart_GS_CE__Form
