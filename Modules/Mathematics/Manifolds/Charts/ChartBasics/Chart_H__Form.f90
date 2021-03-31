module Chart_H__Form

  !-- Chart_Header_Form

  use Basics
  use ManifoldBasics

  implicit none
  private

  type, public :: Chart_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iChart       = 0, &
      nDimensions  = 0
    type ( MeasuredValueForm ), dimension ( : ), pointer :: &
      CoordinateUnit => null ( )
    logical ( KDL ) :: &
      Distributed = .false., &
      AllocatedValues = .false.
    logical ( KDL ), dimension ( : ), pointer :: &
      Periodic => null ( )
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( ), &
      CoordinateSystem => null ( )
    character ( LDL ), dimension ( : ), pointer :: &
      CoordinateLabel => null ( )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    class ( Manifold_H_Form ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic_H
    generic, public :: &
      Initialize_H => InitializeBasic_H
    procedure, private, pass :: &
      Show_C
    generic, public :: &
      Show => Show_C
    final :: &
      Finalize
  end type Chart_H_Form

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = MANIFOLD % MAX_DIMENSIONS

    private :: &
      SetCoordinateSystem


contains


  subroutine InitializeBasic_H &
               ( C, M, Name, Periodic, CommunicatorOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption, nDimensionsOption )

    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    class ( Manifold_H_Form ), intent ( inout ), target :: &
      M
    character ( * ), intent ( in )  :: &
      Name    
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Periodic
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    C % IGNORABILITY     =   M % IGNORABILITY
    C % AllocatedValues  =  .true.

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart'
    end if

    allocate ( C % Name )
    C % Name  =  Name

    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    M % nCharts  =  M % nCharts  +  1
    C % iChart   =  M % nCharts

    if ( present ( CommunicatorOption ) ) then
      C % Distributed  =   .true.
      C % Communicator   =>  CommunicatorOption
    else
      C % Distributed  =   M % Distributed
      C % Communicator   =>  M % Communicator
    end if !-- present Communicator 

    if ( present ( nDimensionsOption ) ) then
      C % nDimensions  =  nDimensionsOption
    else
      C % nDimensions  =  M % nDimensions
    end if

    C % Manifold  =>  M

    call SetCoordinateSystem &
           ( C, Periodic, CoordinateLabelOption, CoordinateSystemOption, &
             CoordinateUnitOption )

  end subroutine InitializeBasic_H


  subroutine Show_C ( C )

    class ( Chart_H_Form ), intent ( in ) :: &
      C

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    associate &
      (  M  =>  C % Manifold, & 
        nD  =>  C % nDimensions )

    call Show ( M % Name,        'Manifold',    C % IGNORABILITY )
    call Show ( C % iChart,      'iChart',      C % IGNORABILITY )
    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )

    call Show ( C % Distributed,       'Distributed', C % IGNORABILITY )
    call Show ( C % Periodic ( : nD ), 'Periodic',    C % IGNORABILITY )

    call Show ( C % CoordinateSystem, 'CoordinateSystem', C % IGNORABILITY )
    call Show ( C % CoordinateLabel ( : nD ), 'CoordinateLabel', &
                C % IGNORABILITY )

    end associate !-- M, etc.

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( Chart_H_Form ), intent ( inout ) :: &
      C

    nullify ( C % Manifold )
    nullify ( C % Communicator )

    if ( .not. associated ( C % Name ) ) &
      return
    if ( C % Name == '' ) &
      return

    if ( C % AllocatedValues ) then
      deallocate ( C % CoordinateLabel )
      deallocate ( C % CoordinateSystem )
      deallocate ( C % Periodic )
      deallocate ( C % CoordinateUnit )
    else
      nullify ( C % CoordinateLabel )
      nullify ( C % CoordinateSystem )
      nullify ( C % Periodic )
      nullify ( C % CoordinateUnit )
    end if !-- AllocatedValues

    call Show ( 'Finalizing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    if ( C % AllocatedValues ) then
      deallocate ( C % Name )
      deallocate ( C % Type )
    else
      nullify ( C % Name )
      nullify ( C % Type )
    end if !-- AllocatedValues

  end subroutine Finalize


  subroutine SetCoordinateSystem &
               ( C, Periodic, CoordinateLabelOption, CoordinateSystemOption, &
                 CoordinateUnitOption )

    class ( Chart_H_Form ), intent ( inout ) :: &
      C
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      Periodic
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption

    associate ( nD => C % nDimensions )

    allocate ( C % Periodic ( MAX_DIMENSIONS ) )
    C % Periodic = .false.
    C % Periodic ( : nD ) = Periodic ( : nD )

    allocate ( C % CoordinateSystem )
    C % CoordinateSystem = 'RECTANGULAR'
    if ( present ( CoordinateSystemOption ) ) &
      C % CoordinateSystem = CoordinateSystemOption
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateSystem, 'CoordinateSystem' )

    allocate ( C % CoordinateLabel ( MAX_DIMENSIONS ) )
    C % CoordinateLabel  =  [ 'X', 'Y', 'Z' ]
    select case ( trim ( C % CoordinateSystem ) )
    case ( 'CYLINDRICAL' )
      C % CoordinateLabel  =  [ 'R_Perp', 'Z     ', 'Phi   ' ]
    case ( 'SPHERICAL' )
      C % CoordinateLabel  =  [ 'R    ', 'Theta', 'Phi  ' ]
    end select !-- CoordinateSystem
    if ( present ( CoordinateLabelOption ) ) &
      C % CoordinateLabel ( : nD ) = CoordinateLabelOption ( : nD )
    call PROGRAM_HEADER % GetParameter &
           ( C % CoordinateLabel ( : nD ), 'CoordinateLabel' )

    allocate ( C % CoordinateUnit ( MAX_DIMENSIONS ) )
    C % CoordinateUnit = [ UNIT % IDENTITY, UNIT % IDENTITY, UNIT % IDENTITY ]
    if ( present ( CoordinateUnitOption ) ) &
      C % CoordinateUnit ( : nD ) = CoordinateUnitOption ( : nD )

    end associate !-- nD

  end subroutine SetCoordinateSystem


end module Chart_H__Form
