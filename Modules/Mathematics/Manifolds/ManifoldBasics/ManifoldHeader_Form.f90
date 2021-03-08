module ManifoldHeader_Form

  use Basics

  implicit none
  private

  type, public :: ManifoldHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nDimensions  = 0
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
      Show_M
    generic, public :: &
      Show => Show_M
    final :: &
      Finalize
  end type ManifoldHeaderForm

    private :: &
      SetDimensionality


contains


  subroutine InitializeBasic &
               ( M, Name, CommunicatorOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( ManifoldHeaderForm ), intent ( inout ) :: &
      M
    character ( * ), intent ( in )  :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    M % IGNORABILITY  =  CONSOLE % INFO_1

    M % AllocatedValues  =  .true.

    if ( .not. associated ( M % Type ) ) then
      allocate ( M % Type )
      M % Type  =  'a Manifold' 
    end if

    allocate ( M % Name )
    M % Name  =  Name

    call Show ( 'Initializing ' // trim ( M % Type ), M % IGNORABILITY )
    call Show ( M % Name, 'Name', M % IGNORABILITY )

    if ( present ( CommunicatorOption ) ) then
      M % IsDistributed  =   .true.
      M % Communicator   =>  CommunicatorOption
    end if !-- present Communicator 

    call SetDimensionality ( M, nDimensionsOption, iDimensionalityOption )

  end subroutine InitializeBasic


  subroutine Show_M ( M )

    class ( ManifoldHeaderForm ), intent ( in ) :: &
      M

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( M % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', M % IGNORABILITY )
    call Show ( M % Name, 'Name', M % IGNORABILITY )

    call Show ( M % IsDistributed, 'IsDistributed', M % IGNORABILITY )
    call Show ( M % nDimensions, 'nDimensions', M % IGNORABILITY )

  end subroutine Show_M


  impure elemental subroutine Finalize ( M )

    type ( ManifoldHeaderForm ), intent ( inout ) :: &
      M

    nullify ( M % Communicator )

    if ( .not. associated ( M % Name ) ) &
      return
    if ( M % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( M % Type ), M % IGNORABILITY )
    call Show ( M % Name, 'Name', M % IGNORABILITY )

    if ( M % AllocatedValues ) then
      deallocate ( M % Name )
      deallocate ( M % Type )
    else
      nullify ( M % Name )
      nullify ( M % Type )
    end if !-- AllocatedValues

  end subroutine Finalize


  subroutine SetDimensionality ( M, nDimensionsOption, iDimensionalityOption )

    class ( ManifoldHeaderForm ), intent ( inout ) :: &
      M
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      iDimensionalityOption

    integer ( KDI ) :: &
      iDimensionality
    character ( LDL ), dimension ( : ), allocatable :: &
      Dimensionality

    if ( present ( nDimensionsOption ) ) then
      M % nDimensions = nDimensionsOption
      call Show ( M % nDimensions, 'nDimensions', M % IGNORABILITY )
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
      call Show ( 'ManifoldHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    select case ( trim ( Dimensionality ( iDimensionality ) ) )
    case ( '1D' )
      M % nDimensions = 1
    case ( '2D' )
      M % nDimensions = 2
    case ( '3D' )
      M % nDimensions = 3
    case default
      call Show ( 'PROGRAM_HEADER % Dimensionality not recognized', &
                  CONSOLE % WARNING )
      call Show ( 'ManifoldHeader_Form', 'module', CONSOLE % WARNING )
      call Show ( 'SetDimensionality', 'subroutine', CONSOLE % WARNING )
      call Show ( 'Defaulting to 3D', CONSOLE % WARNING )
      M % nDimensions = 3
    end select !-- Dimensionality

  end subroutine SetDimensionality


end module ManifoldHeader_Form
