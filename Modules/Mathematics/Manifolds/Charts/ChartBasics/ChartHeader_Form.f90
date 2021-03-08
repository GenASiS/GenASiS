module ChartHeader_Form

  use Basics
  use ManifoldBasics

  implicit none
  private

  type, public :: ChartHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iChart       = 0, &
      nDimensions  = 0
    logical ( KDL ) :: &
      IsDistributed = .false., &
      AllocatedValues = .false.
    logical ( KDL ), dimension ( : ), pointer :: &
      IsPeriodic => null ( )
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( ), &
      CoordinateSystem => null ( )
    character ( LDL ), dimension ( : ), pointer :: &
      CoordinateLabel => null ( )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    class ( ManifoldHeaderForm ), pointer :: &
      Manifold => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    generic, public :: &
      Initialize => InitializeBasic
    procedure, private, pass :: &
      Show_C
    generic, public :: &
      Show => Show_C
    final :: &
      Finalize
  end type ChartHeaderForm

    integer ( KDI ), private, parameter :: &
      MAX_DIMENSIONS = MANIFOLD % MAX_DIMENSIONS

    private :: &
      SetCoordinateSystem


contains


  subroutine InitializeBasic &
               ( C, M, IsPeriodic, iChart, CommunicatorOption, &
                 CoordinateLabelOption, CoordinateSystemOption, &
                 nDimensionsOption )

    class ( ChartHeaderForm ), intent ( inout ) :: &
      C
    class ( ManifoldHeaderForm ), intent ( in ), target :: &
      M
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    character ( 2 ) :: &
      ChartNumber

    C % IGNORABILITY  =   M % IGNORABILITY
        C % Manifold  =>  M
          C % iChart  =   iChart

    C % AllocatedValues = .true.

    if ( .not. associated ( C % Type ) ) then
      allocate ( C % Type )
      C % Type = 'a Chart'
    end if

    allocate ( C % Name )
    write ( ChartNumber, fmt = '(i2.2)' ) iChart
    C % Name = 'Chart_' // ChartNumber // '_' // trim ( M % Name ) 

    call Show ( 'Initializing ' // trim ( C % Type ), C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    if ( present ( CommunicatorOption ) ) then
      C % IsDistributed  =   .true.
      C % Communicator   =>  CommunicatorOption
    else
      C % IsDistributed  =   M % IsDistributed
      C % Communicator   =>  M % Communicator
    end if !-- present Communicator 

    if ( present ( nDimensionsOption ) ) then
      C % nDimensions  =  nDimensionsOption
    else
      C % nDimensions  =  M % nDimensions
    end if

    call SetCoordinateSystem &
           ( C, IsPeriodic, CoordinateLabelOption, CoordinateSystemOption )

  end subroutine InitializeBasic


  subroutine Show_C ( C )

    class ( ChartHeaderForm ), intent ( in ) :: &
      C

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    associate ( nD => C % nDimensions )

    call Show ( C % IsDistributed, 'IsDistributed', C % IGNORABILITY )
    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )

    call Show ( C % IsPeriodic ( : nD ), 'IsPeriodic', C % IGNORABILITY )

    call Show ( C % CoordinateSystem, 'CoordinateSystem', C % IGNORABILITY )
    call Show ( C % CoordinateLabel ( : nD ), 'CoordinateLabel', &
                C % IGNORABILITY )

    end associate !-- nD

  end subroutine Show_C


  impure elemental subroutine Finalize ( C )

    type ( ChartHeaderForm ), intent ( inout ) :: &
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
      deallocate ( C % IsPeriodic )
    else
      nullify ( C % CoordinateLabel )
      nullify ( C % CoordinateSystem )
      nullify ( C % IsPeriodic )
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
               ( C, IsPeriodic, CoordinateLabelOption, CoordinateSystemOption )

    class ( ChartHeaderForm ), intent ( inout ) :: &
      C
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsPeriodic
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption

    associate ( nD => C % nDimensions )

    allocate ( C % IsPeriodic ( MAX_DIMENSIONS ) )
    C % IsPeriodic = .false.
    C % IsPeriodic ( : nD ) = IsPeriodic ( : nD )

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

    end associate !-- nD

  end subroutine SetCoordinateSystem


end module ChartHeader_Form
