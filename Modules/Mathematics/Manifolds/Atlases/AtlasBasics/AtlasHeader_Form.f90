!-- AtlasHeader handles metadata of an Atlas.

module AtlasHeader_Form

  use Basics
  use ATLAS_Singleton
  use Connectivity_Form

  implicit none
  private

  type, public :: AtlasHeaderForm
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nDimensions = 0, &
      nBoundaries = 0, &
      nFields = 0
    logical ( KDL ) :: &
      IsDistributed = .false., &
      AllocatedValues = .false.
    character ( LDF ), pointer :: &
      Type => null ( ), &
      Name => null ( )
    character ( LDL ), dimension ( : ), pointer :: &
      BoundaryName => null ( )
    character ( LDL ), dimension ( :, : ), pointer :: &
      BoundaryCondition => null ( )
    type ( CommunicatorForm ), pointer :: &
      Communicator => null ( )
    type ( ConnectivityForm ), pointer :: &
      Connectivity => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, private, pass :: &
      InitializeClone
    generic, public :: &
      Initialize => InitializeBasic, InitializeClone
    procedure, public, pass :: &
      SetBoundaryConditionsFace
    procedure, public, pass :: &
      SetBoundaryConditionsEdge
    procedure, private, pass :: &
      ShowHeader
    generic, public :: &
      Show => ShowHeader
    final :: &
      Finalize
  end type AtlasHeaderForm

    private :: &
      SetDimensionality, &
      SetDefaultBoundaries, &
      ShowBoundaryConditions

contains


  subroutine InitializeBasic &
               ( A, Name, CommunicatorOption, IncludeFacesOption, &
                 IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    character ( * ), intent ( in )  :: &
      Name
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    logical ( KDL ), intent ( in ), optional :: &
      IncludeFacesOption, &
      IncludeEdgesOption
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption, &
      nDimensionsOption, &
      iDimensionalityOption

    A % IGNORABILITY = CONSOLE % INFO_1

    if ( .not. associated ( A % Type ) ) then
      allocate ( A % Type )
      A % Type = 'an Atlas' 
    end if

    allocate ( A % Name )
    A % Name = Name

    call Show ( 'Initializing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    if ( present ( CommunicatorOption ) ) then
      A % IsDistributed = .true.
      A % Communicator => CommunicatorOption
    end if !-- present Communicator 

    call SetDimensionality ( A, nDimensionsOption, iDimensionalityOption )

    A % AllocatedValues = .true.

    allocate ( A % Connectivity )
    call A % Connectivity % Initialize &
           ( A % nDimensions, IncludeFacesOption, IncludeEdgesOption )

    call SetDefaultBoundaries ( A, nExcisionsOption )

  end subroutine InitializeBasic


  subroutine InitializeClone ( A, A_Source )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A_Source

    A % IGNORABILITY         =  CONSOLE % INFO_5  !-- NOT COPIED!
    A % nDimensions          =  A_Source % nDimensions
    A % nBoundaries          =  A_Source % nBoundaries
    A % IsDistributed        =  A_Source % IsDistributed
    A % AllocatedValues      =  .false.  !-- NOT COPIED!
    A % Type                 => A_Source % Type
    A % Name                 => A_Source % Name
    A % BoundaryName         => A_Source % BoundaryName
    A % BoundaryCondition    => A_Source % BoundaryCondition
    A % Communicator         => A_Source % Communicator
    A % Connectivity         => A_Source % Connectivity

  end subroutine InitializeClone


  subroutine SetBoundaryConditionsFace &
               ( A, BoundaryCondition, iDimension, BoundaryNameOption, &
                 iBoundaryOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    character ( * ), dimension ( 2 ), intent ( in ) :: &
      BoundaryCondition  !-- [ Inner, Outer ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryNameOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    integer ( KDI ) :: &
      iB  !-- iBoundary

    if ( A % Connectivity % nFaces == 0 ) then
      call Show ( 'Faces not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( iDimension > A % nDimensions ) then
      call Show ( 'Selected iDimension > nDimensions', CONSOLE % ERROR )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsFace', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > A % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryNameOption ) ) then
          call Show ( 'BoundaryName not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryNameOption ) ) then
          call Show ( 'BoundaryName required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryNameOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when BoundaryName present', &
                    CONSOLE % ERROR )
        call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        A % BoundaryName ( iBoundaryOption ) = BoundaryNameOption
      end if
    end if

    associate &
      ( C  => A % Connectivity, &
        iD => iDimension )
    A % BoundaryCondition ( C % iaInner ( iD ), iB ) = BoundaryCondition ( 1 )
    A % BoundaryCondition ( C % iaOuter ( iD ), iB ) = BoundaryCondition ( 2 )
    end associate !-- C, etc.

  end subroutine SetBoundaryConditionsFace


  subroutine SetBoundaryConditionsEdge &
               ( A, BoundaryCondition, iDimension, BoundaryNameOption, &
                 iBoundaryOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    character ( * ), dimension ( 4 ), intent ( in ) :: &
      BoundaryCondition  !-- [ InnerInner, OuterInner, InnerOuter, OuterOuter ]
    integer ( KDI ), intent ( in ) :: &
      iDimension
    character ( * ), intent ( in ), optional :: &
      BoundaryNameOption
    integer ( KDI ), intent ( in ), optional :: &
      iBoundaryOption

    integer ( KDI ) :: &
      iD, jD, kD, &  !-- jDimension, etc.
      iB  !-- iBoundary

    iD = iDimension
    jD = mod ( iD, 3 ) + 1
    kD = mod ( jD, 3 ) + 1
 
    if ( A % Connectivity % nEdges == 0 ) then
      call Show ( 'Edges not included in Connectivity', CONSOLE % ERROR )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    if ( jD > A % nDimensions .or. kD > A % nDimensions ) then
      call Show ( 'Selected jDimension or kDimension > nDimensions', &
                  CONSOLE % ERROR )
      call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetBoundaryConditionsEdge', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end if

    iB = 1
    if ( present ( iBoundaryOption ) ) then
      if ( iBoundaryOption > A % nBoundaries ) then
        call Show ( 'Selected iBoundary > nBoundaries', CONSOLE % ERROR )
        call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end if
      if ( iBoundaryOption == 1 ) then
        if ( present ( BoundaryNameOption ) ) then
          call Show ( 'BoundaryName not allowed for iBoundary == 1', &
                      CONSOLE % ERROR )
          call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      else
        if ( .not. present ( BoundaryNameOption ) ) then
          call Show ( 'BoundaryName required for iBoundary > 1', &
                      CONSOLE % ERROR )
          call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
          call Show ( 'SetBoundaryConditionsEdge', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )           
        end if
      end if
      iB = iBoundaryOption
    end if

    if ( present ( BoundaryNameOption ) ) then
      if ( .not.present ( iBoundaryOption ) ) then
        call Show ( 'Argument iBoundary required when BoundaryName present', &
                    CONSOLE % ERROR )
        call Show ( 'AtlasHeader_Form', 'module', CONSOLE % ERROR )
        call Show ( 'SetBoundaryConditionsFace', 'subroutine', &
                    CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )         
      else
        A % BoundaryName ( iBoundaryOption ) = BoundaryNameOption
      end if
    end if

    associate ( C  => A % Connectivity )
    A % BoundaryCondition ( C % iaInnerInner ( iD ), iB ) &
      = BoundaryCondition ( 1 )
    A % BoundaryCondition ( C % iaOuterInner ( iD ), iB ) &
      = BoundaryCondition ( 2 )
    A % BoundaryCondition ( C % iaInnerOuter ( iD ), iB ) &
      = BoundaryCondition ( 3 )
    A % BoundaryCondition ( C % iaOuterOuter ( iD ), iB ) &
      = BoundaryCondition ( 4 )
    end associate !-- C, etc.

  end subroutine SetBoundaryConditionsEdge


  subroutine ShowHeader ( A )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A

    character ( LDL ), dimension ( : ), allocatable :: &
      TypeWord

    call Split ( A % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

    call Show ( A % nDimensions, 'nDimensions', A % IGNORABILITY )
    call Show ( A % IsDistributed, 'IsDistributed', A % IGNORABILITY )
    call ShowBoundaryConditions ( A )

  end subroutine ShowHeader


  impure elemental subroutine Finalize ( A )

    type ( AtlasHeaderForm ), intent ( inout ) :: &
      A

    nullify ( A % Communicator )

    if ( A % AllocatedValues ) then
      deallocate ( A % Connectivity )
      deallocate ( A % BoundaryCondition )
      deallocate ( A % BoundaryName )
    else
      nullify ( A % Connectivity )
      nullify ( A % BoundaryCondition )
      nullify ( A % BoundaryName )
    end if !-- AllocatedValues

    if ( A % Name == '' ) return

    call Show ( 'Finalizing ' // trim ( A % Type ), A % IGNORABILITY )
    call Show ( A % Name, 'Name', A % IGNORABILITY )

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
      call PROGRAM_HEADER % Communicator % Synchronize ( )
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


  subroutine SetDefaultBoundaries ( A, nExcisionsOption )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption

    A % nBoundaries = 1
    if ( present ( nExcisionsOption ) ) &
      A % nBoundaries = 1 + nExcisionsOption 

    associate ( C => A % Connectivity )
    allocate ( A % BoundaryCondition ( C % nConnections, A % nBoundaries ) )
    allocate ( A % BoundaryName ( A % nBoundaries ) )

    A % BoundaryName = ''
    A % BoundaryName ( 1 ) = 'Extent' 

    A % BoundaryCondition = ''
    A % BoundaryCondition ( :, 1 ) = 'PERIODIC'

    end associate !-- C

  end subroutine SetDefaultBoundaries


  subroutine ShowBoundaryConditions ( A )

    class ( AtlasHeaderForm ), intent ( inout ) :: &
      A

    integer ( KDI ) :: &
      iB, &  !-- iBoundary
      iD, jD, kD  !-- iDimension, etc.

    associate &
      ( C  => A % Connectivity, &
        BC => A % BoundaryCondition ( :, : ), &
        BN => A % BoundaryName )

    call Show ( 'Boundary conditions', A % IGNORABILITY )
    call Show ( A % nBoundaries, 'nBoundaries', A % IGNORABILITY )

    do iB = 1, A % nBoundaries
      call Show ( A % BoundaryName ( iB ), 'BoundaryName', A % IGNORABILITY )
      call Show ( iB, 'iBoundary', A % IGNORABILITY )
  
      if ( C % nFaces > 0 ) then
          do iD = 1, A % nDimensions
            call Show ( iD, 'Faces, iDimension', &
                        A % IGNORABILITY )
            associate &
              ( iaI => C % iaInner ( iD ), &
                iaO => C % iaOuter ( iD ) )
            call Show ( [ BC ( iaI, iB ), BC ( iaO, iB ) ], &
                        '[ Inner, Outer ]', A % IGNORABILITY )
            end associate !-- iaI, etc.
          end do !-- iD
      end if

      if ( C % nEdges > 0 ) then
          do iD = 1, A % nDimensions
            jD = mod ( iD, 3 ) + 1
            kD = mod ( jD, 3 ) + 1
            if ( jD > A % nDimensions .or. kD > A % nDimensions ) &
              cycle
            call Show ( iD, 'Edges parallel to iDimension', &
                        A % IGNORABILITY )
            associate &
              ( iaII => C % iaInnerInner ( iD ), &
                iaOI => C % iaOuterInner ( iD ), &
                iaIO => C % iaInnerOuter ( iD ), &
                iaOO => C % iaOuterOuter ( iD ) )
            call Show ( [ BC ( iaII, iB ), BC ( iaOI, iB ), &
                          BC ( iaIO, iB ), BC ( iaOO, iB ) ], &
                        '[ InnerInner, OuterInner, InnerOuter, OuterOuter ]', &
                        A % IGNORABILITY )
            end associate !-- iaII, etc.
          end do !-- iD
      end if

    end do !-- iB

    end associate !-- C, etc.

  end subroutine ShowBoundaryConditions


end module AtlasHeader_Form
