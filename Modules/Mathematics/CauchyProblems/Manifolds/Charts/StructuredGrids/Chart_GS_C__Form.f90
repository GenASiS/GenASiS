module Chart_GS_C__Form
  
  !-- Chart_GridStructured_Central_Form

  use Basics
  use Chart_GS__Form

  implicit none
  private

  type, public, extends ( Chart_GS_Form ) :: Chart_GS_C_Form
    integer ( KDI ) :: &
      nCellsPolar
    real ( KDR ) :: &
      RadiusMax = 0.0_KDR, &    !-- should be set by a descendant
      RadiusScale = 0.0_KDR, &  !-- should be set by a descendant
      RadialRatio, &            !-- nCellsRadial / nCellsPolar
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
    procedure, private, pass :: &
      SetCore
  end type Chart_GS_C_Form


contains


  subroutine Initialize_GS_C &
               ( C, RadiusMin, CoordinateUnitOption, RadiusMaxOption, &
                 RadialRatioOption, nGhostLayersOption, nCellsPolarOption, &
                 nEqualOption  )

    class ( Chart_GS_C_Form ), intent ( inout ) :: &
      C
    real ( KDR ), intent ( in ) :: &
      RadiusMin
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), intent ( in ), optional :: &
      RadiusMaxOption, &
      RadialRatioOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nCellsPolarOption, &
      nEqualOption

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
    logical ( KDL ), dimension ( 3 ) :: &
      Periodic
    character ( LDL ) :: &
      CoordinateSystem
    character ( LDL ), dimension ( 3 ) :: &
      Spacing

    Pi  =  CONSTANT % PI

    CoordinateSystem  =  'SPHERICAL'

    Spacing        =  'EQUAL'
    Spacing ( 1 )  =  'PROPORTIONAL'
    
    Periodic  =  .false.
    Periodic ( 3 )  =  .true.

    if ( C % RadiusMax  ==  0.0_KDR ) &
      C % RadiusMax  =  10.0_KDR
    if ( C % RadiusScale  ==  0.0_KDR ) &
      C % RadiusScale  =  C % RadiusMax  /  8.0_KDR

    MinCoordinate  =  [ RadiusMin,     0.0_KDR,      0.0_KDR ]
    MaxCoordinate  =  [ C % RadiusMax,      Pi, 2.0_KDR * Pi ]

    C % nCellsPolar  =  128
    if ( present ( nCellsPolarOption ) ) &
      C % nCellsPolar = nCellsPolarOption
    call PROGRAM_HEADER % GetParameter ( C % nCellsPolar, 'nCellsPolar' )

    call C % SetCore ( )

    C % RadialRatio  =  1
    if ( present ( RadialRatioOption ) ) &
      C % RadialRatio = RadialRatioOption
    call PROGRAM_HEADER % GetParameter ( C % RadialRatio, 'RadialRatio' )

    nCellsRadial     =  C % RadialRatio * C % nCellsPolar !-- Aim for RadiusMax
    nCellsPolar      =  C % nCellsPolar
    nCellsAzimuthal  =  2 * nCellsPolar
 
    C % MinWidth  =  C % RadiusScale  *  Pi / nCellsPolar

    nCells  =  [ nCellsRadial, 1, 1 ]
    if ( C % nDimensions  >  1 ) &
      nCells ( 2 )  =  nCellsPolar
    if ( C % nDimensions > 2 ) &
      nCells ( 3 )  =  nCellsAzimuthal

    Ratio        =  0.0_KDR
    Ratio ( 1 )  =  Pi / nCellsPolar  !-- dTheta

    Scale        =  0.0_KDR
    Scale ( 1 )  =  C % RadiusScale

    nBricks  =  [ 1, 1, 1 ]
    if ( C % Distributed ) &
      nBricks ( 1 )  =  C % Communicator % Size  !-- spherical shells

    call C % Chart_GS_Form % Initialize &
           ( SpacingOption = Spacing, &
             CoordinateSystemOption = CoordinateSystem, &
             PeriodicOption = Periodic, &
             CoordinateUnitOption = CoordinateUnitOption, &
             MinCoordinateOption = MinCoordinate, &
             MaxCoordinateOption = MaxCoordinate, &
             RatioOption = Ratio, &
             ScaleOption = Scale, &
             nCellsOption = nCells, &
             nGhostLayersOption = nGhostLayersOption, &
             nBricksOption = nBricks, &
             nEqualOption = nEqualOption )

    if ( C % nBricks ( 2 )  /=  1  .or.  C % nBricks ( 3 ) /= 1 ) then
      call Show ( 'Decomposition in angle not allowed', CONSOLE % ERROR )
      call Show ( 'Do not use nBricks command line option', CONSOLE % ERROR )
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

    call Show ( 'Chart_GS_C parameters' )
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


  subroutine SetCore ( C )

    class ( Chart_GS_C_Form ), intent ( inout ) :: &
      C

  end subroutine SetCore


end module Chart_GS_C__Form
