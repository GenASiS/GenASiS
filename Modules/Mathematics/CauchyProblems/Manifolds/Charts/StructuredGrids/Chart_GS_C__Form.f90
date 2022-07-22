module Chart_GS_C__Form
  
  !-- Chart_GridStructured_Central_Form

  use Basics
  use Chart_GS__Form

  implicit none
  private

  type, public, extends ( Chart_GS_Form ) :: Chart_GS_C_Form
    integer ( KDI ) :: &
      nCellsPolar = 0
    real ( KDR ) :: &
      RadiusMax, &    !-- should be set by a descendant
      RadiusScale, &  !-- should be set by a descendant
      RadialRatio, &  !-- nCellsRadial / nCellsPolar
      MinWidth
  contains
    procedure, private, pass :: &
      Initialize_GS_C
    generic, public :: &
      Initialize => Initialize_GS_C
    procedure, private, pass :: &
      Show_C
    final :: &
      Finalize
    procedure, public, pass :: &
      SetPolar
  end type Chart_GS_C_Form


contains


  subroutine Initialize_GS_C &
               ( C, RadiusMin, RadiusMax, RadiusScale, RadialRatio, &
                 CommunicatorOption, NameOption, CoordinateUnitOption, &
                 nGhostLayersOption, nCellsPolarOption, nEqualOption, &
                 nDimensionsOption )

    class ( Chart_GS_C_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMin, &
      RadiusMax, &
      RadiusScale, &
      RadialRatio
    type ( CommunicatorForm ), intent ( in ), optional :: &
      CommunicatorOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    type ( QuantityForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption, &
      nDimensionsOption

    integer ( KDI ) :: &
      nCellsRadial, &
      nCellsPolar, &
      nCellsAzimuthal
    integer ( KDI ), dimension ( 3 ) :: &
      nCells, &
      nBricks
    real ( KDR ) :: &
      Pi
    real ( KDR ), dimension ( 3 ) :: &
      MinCoordinate, &
      MaxCoordinate, &
      Ratio, &
      Scale
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    if ( C % Type  ==  '' ) &
      C % Type  =  'a Chart_GS_C'

    Pi  =  CONSTANT % PI

    CoordinateSystem  =  'SPHERICAL'

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    C % RadiusMax    =  RadiusMax
    C % RadiusScale  =  RadiusScale

    MinCoordinate  =  [ RadiusMin,     0.0_KDR,      0.0_KDR ]
    MaxCoordinate  =  [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]

    if ( C % nCellsPolar  ==  0 ) &
      call C % SetPolar ( nCellsPolarOption )

    C % RadialRatio  =  RadialRatio

    nCellsRadial     =  C % RadialRatio * C % nCellsPolar !-- Aim for RadiusMax
    nCellsPolar      =  C % nCellsPolar
    nCellsAzimuthal  =  2 * nCellsPolar
 
    nCells  =  [ nCellsRadial, nCellsPolar, nCellsAzimuthal ]

    if ( present ( nEqualOption ) ) then
      C % MinWidth  =  C % RadiusScale  /  nEqualOption
    else
      C % MinWidth  =  C % RadiusScale  *  Pi / nCellsPolar
    end if

    Ratio  =  0.0_KDR
    if ( present ( nEqualOption ) ) then
      Ratio ( 1 )  =  1.0_KDR / nEqualOption
    else
      Ratio ( 1 )  =  Pi / nCellsPolar  !-- dTheta
    end if

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusScale

    nBricks  =  [ 1, 1, 1 ]
    if ( present ( CommunicatorOption ) ) &
      nBricks ( 1 )  =  CommunicatorOption % Size  !-- spherical shells

    call C % Chart_GS_Form % Initialize &
           ( CommunicatorOption = CommunicatorOption, &
             SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             NameOption = NameOption, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption, &
             nBricksOption = nBricks, &
             nEqualOption = nEqualOption, &
             nDimensionsOption = nDimensionsOption )

    if ( C % nBricks ( 2 )  /=  1  .or.  C % nBricks ( 3 ) /= 1 ) then
      call Show ( 'Decomposition in angle not allowed', CONSOLE % ERROR )
      call Show ( 'Do not use nBricks command line option', CONSOLE % ERROR )
      call Show ( 'Chart_GS_C__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Initialize_GS_C', 'module', CONSOLE % ERROR )
      call PROGRAM_HEADER % Communicator % Synchronize ( )
      call PROGRAM_HEADER % Abort ( )
    end if

    ! if ( C % nDimensions > 1 .and. C % nCells ( 2 ) /= C % nCellsPolar ) then
    !   call Show ( 'Choose nBricks such that nCells ( 2 ) need not be', &
    !               CONSOLE % ERROR )
    !   call Show ( 'changed from requested nCellsPolar', &
    !               CONSOLE % ERROR )
    !   call Show ( C % nCellsPolar, 'nCellsPolar', CONSOLE % ERROR )
    !   call Show ( C % nBricks ( 2 ), 'nBricks ( 2 )', CONSOLE % ERROR )
    !   call Show ( mod ( C % nCellsPolar, C % nBricks ( 2 ) ), &
    !               'mod ( nCellsPolar, nBricks ( 2 ) )', CONSOLE % ERROR )
    !   call Show ( C % nCells ( 2 ), 'nCells ( 2 )', CONSOLE % ERROR )
    !   call Show ( 'InitializeTemplate_C', 'subroutine', CONSOLE % ERROR )
    !   call Show ( 'Chart_SLD_C__Template', 'module', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Communicator % Synchronize ( )
    !   call PROGRAM_HEADER % Abort ( )
    ! end if
      
  end subroutine Initialize_GS_C


  subroutine Show_C ( C )

    class ( Chart_GS_C_Form ), intent ( in ) :: &
      C

    call C % Chart_GS_Form % Show ( )

    call Show ( 'Chart_GS_C Proper Parameters' )
    call Show ( C % nCellsPolar, 'nCellsPolar', C % IGNORABILITY )
    call Show ( C % RadiusScale, C % CoordinateUnit ( 1 ), 'RadiusScale', &
                C % IGNORABILITY )
    call Show ( C % MinWidth, C % CoordinateUnit ( 1 ), 'MinWidth', &
                C % IGNORABILITY )
    call Show ( C % RadialRatio, 'RadialRatio' )
    call Show ( C % RadiusMax, C % CoordinateUnit ( 1 ), &
                'RadiusMax requested', C % IGNORABILITY )
    call Show ( C % MaxCoordinate ( 1 ), C % CoordinateUnit ( 1 ), &
                'RadiusMax actual', C % IGNORABILITY )

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_GS_C_Form ), intent ( inout ) :: &
      C

  end subroutine Finalize


  subroutine SetPolar ( C, nCellsPolarOption )

    class ( Chart_GS_C_Form ), intent ( inout ) :: &
      C
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption

    C % nCellsPolar  =  128
    if ( present ( nCellsPolarOption ) ) &
      C % nCellsPolar  =  nCellsPolarOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsPolar, 'nCellsPolar' )

  end subroutine SetPolar


end module Chart_GS_C__Form
