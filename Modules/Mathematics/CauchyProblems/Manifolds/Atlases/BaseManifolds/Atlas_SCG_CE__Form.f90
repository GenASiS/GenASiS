module Atlas_SCG_CE__Form

  !-- Atlas_SingleChartGrid_CentralExcision_Form

  use Basics
  use Charts
  use Atlas_SCG__Form

  implicit none
  private

  type, public, extends ( Atlas_SCG_Form ) :: Atlas_SCG_CE_Form
    class ( Chart_GS_CE_Form ), pointer :: &
      Chart_GS_CE => null ( )
  contains
    procedure, private, pass :: &
      Initialize_SCG_CE
    procedure, private, pass :: &
      Initialize_SCG_CE_SA  !-- SphericalAverage
    generic, public :: &
      Initialize => Initialize_SCG_CE, Initialize_SCG_CE_SA
    final :: &
      Finalize
  end type Atlas_SCG_CE_Form


contains


  subroutine Initialize_SCG_CE &
               ( A, RadiusMax, RadiusExcision, CommunicatorOption, NameOption, &
                 CoordinateUnitOption, RadialRatioOption, nGhostLayersOption, &
                 nCellsPolarOption, nDimensionsOption )

    class ( Atlas_SCG_CE_Form ), intent ( inout ), target :: &
      A
    real ( KDR ), intent ( in ) :: &
      RadiusMax, &
      RadiusExcision
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
      nDimensionsOption

    logical :: &
      PreviouslyAllocated

    if ( A % Type  ==  '' ) &
      A % Type  =  'an Atlas_SCG_CE'

    if ( .not. allocated ( A % Chart ) ) &
      allocate ( A % Chart ( 1 ) )
    
    if ( allocated ( A % Chart ( 1 ) % Element ) ) then
      PreviouslyAllocated  =  .true.
    else
      PreviouslyAllocated  =  .false.
      allocate ( Chart_GS_CE_Form :: A % Chart ( 1 ) % Element )
    end if

    call A % Atlas_SCG_Form % Initialize ( NameOption = NameOption )

    if ( .not. PreviouslyAllocated ) then

      select type ( C  =>  A % Chart ( 1 ) % Element )
      class is ( Chart_GS_CE_Form )

      call C % Initialize &
             ( RadiusMax, RadiusExcision, &
               CommunicatorOption = CommunicatorOption, &
               NameOption = NameOption, &
               CoordinateUnitOption = CoordinateUnitOption, &
               RadialRatioOption = RadialRatioOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nCellsPolarOption = nCellsPolarOption, &
               nDimensionsOption = nDimensionsOption )

      end select !--  C

    end if !-- PreviouslyAllocated  

    select type ( C  =>  A % Chart ( 1 ) % Element )
    class is ( Chart_GS_CE_Form )
      A % Chart_GS_CE  =>  C
    end select !-- C
      
  end subroutine Initialize_SCG_CE


  subroutine Initialize_SCG_CE_SA ( A, A_S )

    class ( Atlas_SCG_CE_Form ), intent ( inout ) :: &
      A
    class ( Atlas_SCG_CE_Form ), intent ( in ) :: &
      A_S  !-- Source

    associate ( C_S  =>  A_S % Chart_GS_CE )

    call A % Initialize &
           ( RadiusMax = C_S % RadiusMax, &
             RadiusExcision = C_S % RadiusScale, &
             CommunicatorOption = C_S % Communicator, &
             NameOption = trim ( A_S % Name ) // '_SA', &
             CoordinateUnitOption = C_S % CoordinateUnit, &
             RadialRatioOption = C_S % RadialRatio, &
             nGhostLayersOption = C_S % nGhostLayers, &
             nCellsPolarOption = C_S % nCellsPolar, &
             nDimensionsOption = 1 )

    end associate !-- C_S

  end subroutine Initialize_SCG_CE_SA


  impure elemental subroutine Finalize ( A )

    type ( Atlas_SCG_CE_Form ), intent ( inout ) :: &
      A

    nullify ( A % Chart_GS_CE )

  end subroutine Finalize


end module Atlas_SCG_CE__Form
