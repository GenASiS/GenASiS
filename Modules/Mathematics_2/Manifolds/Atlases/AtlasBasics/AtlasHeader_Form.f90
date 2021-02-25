!-- AtlasHeader handles metadata of an Atlas.

module AtlasHeader_Form

  use Basics

  implicit none
  private

  type, public :: AtlasHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nDimensions = 0
    logical ( KDL ) :: &
      IsDistributed = .false., &
      AllocatedValues = .false.
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    generic, public :: &
      Initialize => InitializeBasic
    procedure, private, pass :: &
      Show_AH
    generic, public :: &
      Show => Show_AH
    final :: &
      Finalize
  end type AtlasHeaderForm

    private :: &
      SetDimensionality


contains


  subroutine InitializeBasic &
               ( A, Name, CommunicatorOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    character ( * ), intent ( in )  :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    A % IGNORABILITY  =  CONSOLE % INFO_1

    A % AllocatedValues  =  .true.

    if ( .not. associated ( A % Type ) ) then
      allocate ( A % Type )
      A % Type  =  'an Atlas' 
    end if

    allocate ( A % Name )
    A % Name  =  Name

    call Show ( 'Initializing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    if ( present ( CommunicatorOption ) ) then
      A % IsDistributed  =   .true.
      A % Communicator   =>  CommunicatorOption
    end if !-- present Communicator 

    call SetDimensionality ( A, nDimensionsOption, iDimensionalityOption )

  end subroutine InitializeBasic


  subroutine Show_AH ( A )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( A % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    call Show ( A % nDimensions, 'nDimensions', A % IGNORABILITY )
    call Show ( A % IsDistributed, 'IsDistributed', A % IGNORABILITY )

  end subroutine Show_AH


  impure elemental subroutine Finalize ( A )

    type ( AtlasHeaderForm ), intent ( inout ) :: &
      A

    nullify ( A % Communicator )

    if ( A % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    if ( A % AllocatedValues ) then
      deallocate ( A % Name )
      deallocate ( A % Type )
    else
      nullify ( A % Name )
      nullify ( A % Type )
    end if !-- AllocatedValues

  end subroutine Finalize


  subroutine SetDimensionality ( A, nDimensionsOption, iDimensionalityOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    integer ( KDI ) :: &
      iDimensionality
    character ( LDL ), dimension ( : ), allocatable :: &
      Dimensionality

    if ( present ( nDimensionsOption ) ) then
      A % nDimensions = nDimensionsOption
      call Show ( A % nDimensions, 'nDimensions', A % IGNORABILITY )
      return
    end if

    !-- Allow for specification of base manifold and bundle 
    !   dimensionalities; take the first element here, the dimensionality
    !   of the base manifold
    call Split ( PROGRAM_HEADER % Dimensionality, '_', Dimensionality )

    iDimensionality = 1
    if ( present ( iDimensionalityOption ) ) &
      iDimensionality = iDimensionalityOption

    if ( iDimensionality > size ( Dimensionality ) ) then
      call Show ( 'Too few dimensionalities specified', CONSOLE % ERROR )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    select case ( trim ( Dimensionality ( iDimensionality ) ) )
    case ( '1D' )
      A % nDimensions = 1
    case ( '2D' )
      A % nDimensions = 2
    case ( '3D' )
      A % nDimensions = 3
    case default
      call Show ( 'PROGRAM_HEADER % Dimensionality not recognized', &
                  CONSOLE % WARNING )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % WARNING )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % WARNING )
      call Show ( 'Defaulting to 3D', CONSOLE % WARNING )
      A % nDimensions = 3
    end select !-- Dimensionality

    call Show ( A % nDimensions, 'nDimensions', A % IGNORABILITY )

  end subroutine SetDimensionality


end module AtlasHeader_Form
