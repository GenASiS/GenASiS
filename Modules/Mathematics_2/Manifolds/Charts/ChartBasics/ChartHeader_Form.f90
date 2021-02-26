!-- ChartHeader_Form contains the most basic functionality characterizing a
!   coordinate chart.

module ChartHeader_Form

  use Basics
  use ManifoldBasics

  implicit none
  private

  type, public :: ChartHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      iChart       = 0, &
      nDimensions  = 0, &
      nFields      = 0
    logical ( KDL ) :: &
      IsDistributed = .false., &
      AllocatedValues = .false.
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( )
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
      Show_CH
    generic, public :: &
      Show => Show_CH
    final :: &
      Finalize
  end type ChartHeaderForm


contains


  subroutine InitializeBasic &
               ( C, M, iChart, CommunicatorOption, nDimensionsOption )

    class ( ChartHeaderForm ), intent ( inout ) :: &
      C
    class ( ManifoldHeaderForm ), intent ( in ), target :: &
      M
    integer ( KDI ), intent ( in ) :: &
      iChart
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption

    character ( 2 ) :: &
      ChartNumber

    C % IGNORABILITY  =  M % IGNORABILITY
    C % Manifold  =>  M
    C % iChart  =  iChart

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

  end subroutine InitializeBasic


  subroutine Show_CH ( C )

    class ( ChartHeaderForm ), intent ( in ) :: &
      C

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( C % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', C % IGNORABILITY )
    call Show ( C % Name, 'Name', C % IGNORABILITY )

    call Show ( C % IsDistributed, 'IsDistributed', C % IGNORABILITY )
    call Show ( C % nDimensions, 'nDimensions', C % IGNORABILITY )
    call Show ( C % nFields, 'nFields', C % IGNORABILITY )

  end subroutine Show_CH


  impure elemental subroutine Finalize ( C )

    type ( ChartHeaderForm ), intent ( inout ) :: &
      C

    nullify ( C % Manifold )
    nullify ( C % Communicator )

    if ( .not. associated ( C % Name ) ) return

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


end module ChartHeader_Form
