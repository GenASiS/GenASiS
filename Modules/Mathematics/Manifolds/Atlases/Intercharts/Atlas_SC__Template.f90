!-- Atlas_SC represents an Atlas with a single chart.

module Atlas_SC__Template

  !-- Atlas_SingleChart_Template

  use Basics
  use AtlasBasics
  use Charts
  use Field_ASC__Template

  implicit none
  private

  type, public, extends ( AtlasHeaderForm ), abstract :: Atlas_SC_Template
    class ( ChartTemplate ), allocatable :: &
      Chart
    type ( FieldAtlasPointer ), dimension ( : ), allocatable :: &
      Field
    procedure ( ABCCSLC ), pointer, nopass :: &
      Apply_BC_CSL_Custom => null ( )
  contains
    procedure, private, pass :: &
      InitializeBasic
    procedure, private, pass :: &
      InitializeClone
    procedure, public, pass :: &
      CreateChart
    procedure, public, pass :: &
      AddField
    procedure, public, pass :: &
      ApplyBoundaryConditions
    procedure, public, pass :: &
      ApplyBoundaryConditionsFaces
    procedure, public, pass :: &
      FinalizeTemplate
  end type Atlas_SC_Template

  abstract interface

    subroutine ABCCSLC ( CSL, F, iDimension, iConnection )
      use Basics
      use Charts
      class ( Chart_SL_Template ), intent ( inout ) :: &
        CSL
      class ( StorageForm ), intent ( inout ) :: &
        F  !-- Field
      integer ( KDI ), intent ( in ) :: &
        iDimension, &
        iConnection
    end subroutine ABCCSLC

  end interface

    private :: &
      Apply_BC_CSL_Outflow, &
      Apply_BC_CSL_Reflecting

contains


  subroutine InitializeBasic &
               ( A, Name, CommunicatorOption, IncludeFacesOption, &
                 IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
                 iDimensionalityOption )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
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

    if ( .not. associated ( A % Type ) ) then
      allocate ( A % Type )
      A % Type = 'an Atlas_SC' 
    end if

    call A % AtlasHeaderForm % Initialize &
           ( Name, CommunicatorOption, IncludeFacesOption, &
             IncludeEdgesOption, nExcisionsOption, nDimensionsOption, &
             iDimensionalityOption )

    allocate ( A % Field ( ATLAS % MAX_FIELDS ) )

  end subroutine InitializeBasic


  subroutine InitializeClone ( A, A_Source )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    class ( AtlasHeaderForm ), intent ( in ), target :: &
      A_Source

    call A % AtlasHeaderForm % Initialize ( A_Source )

    select type ( A_Source )
    class is ( Atlas_SC_Template )

    select type ( C_Source => A_Source % Chart )
    type is ( Chart_SLL_Form )
      allocate ( Chart_SLL_Form :: A % Chart )
      select type ( C => A % Chart )
      type is ( Chart_SLL_Form )
        call C % Initialize ( C_Source )
      end select !-- C
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeClone', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C_Source

    class default
      call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeClone', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- A_Source

    allocate ( A % Field ( ATLAS % MAX_FIELDS ) )

  end subroutine InitializeClone


  subroutine CreateChart &
               ( A, SpacingOption, CoordinateLabelOption, &
                 CoordinateSystemOption, CoordinateUnitOption, &
                 MinCoordinateOption, MaxCoordinateOption, RatioOption, &
                 ScaleOption, nCellsOption, nGhostLayersOption, &
                 nDimensionsOption, nEqualOption )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption, &
      CoordinateLabelOption
    character ( * ), intent ( in ), optional :: &
      CoordinateSystemOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      CoordinateUnitOption
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      MinCoordinateOption, &
      MaxCoordinateOption, &
      RatioOption, &
      ScaleOption
    integer ( KDI ), dimension ( : ), intent ( in ), optional :: &
      nCellsOption, &
      nGhostLayersOption
    integer ( KDI ), intent ( in ), optional :: &
      nDimensionsOption, &
      nEqualOption

    integer ( KDI ) :: &
      iD  !-- iDimension
    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic
    character ( LDL ) :: &
      ChartType

    ChartType = 'SINGLE_LEVEL'
    call PROGRAM_HEADER % GetParameter ( ChartType, 'ChartType' )

    select case ( trim ( ChartType ) )
    case ( 'SINGLE_LEVEL' )
      if ( A % IsDistributed ) then
        allocate ( Chart_SLD_Form :: A % Chart )
      else !-- .not. IsDistributed
        allocate ( Chart_SLL_Form :: A % Chart )
      end if !-- IsDistributed
    case default
      call Show ( 'ChartType not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'CreateChart', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- ChartType

    IsPeriodic = .false.
    do iD = 1, A % nDimensions
      associate &
        ( iaI => A % Connectivity % iaInner ( iD ), &
          iaO => A % Connectivity % iaOuter ( iD ) )
      IsPeriodic ( iD ) &
        = ( A % BoundaryCondition ( iaI, 1 ) == 'PERIODIC' ) &
          .and. ( A % BoundaryCondition ( iaO, 1 ) == 'PERIODIC' )
      end associate !-- iaI, iaO
    end do !-- iD

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      call C % Initialize &
             ( A, IsPeriodic, iChart = 1, SpacingOption = SpacingOption, & 
               CoordinateLabelOption = CoordinateLabelOption, &
               CoordinateSystemOption = CoordinateSystemOption, &
               CoordinateUnitOption = CoordinateUnitOption, &
               MinCoordinateOption = MinCoordinateOption, &
               MaxCoordinateOption = MaxCoordinateOption, &
               RatioOption = RatioOption, ScaleOption = ScaleOption, &
               nCellsOption = nCellsOption, &
               nGhostLayersOption = nGhostLayersOption, &
               nDimensionsOption = nDimensionsOption, &
               nEqualOption = nEqualOption )
    end select !-- C

  end subroutine CreateChart


  subroutine AddField ( A, FA )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    class ( Field_ASC_Template ), intent ( in ), target :: &
      FA

    A % nFields = A % nFields + 1
    A % Field ( A % nFields ) % Pointer => FA

    select type ( C => A % Chart )
    class is ( Chart_SL_Template )
      select type ( FC => FA % Chart )
      class is ( Field_CSL_Template )
        call C % AddField ( FC )
      end select !-- FC
    class default
      call Show ( 'Chart type not recognized', CONSOLE % ERROR )
      call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
      call Show ( 'AddField', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine AddField


  subroutine ApplyBoundaryConditions ( A, F, iDimension, iConnection )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    class ( StorageForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      iB  !-- iBoundary

    do iB = 1, A % nBoundaries
      !-- FIXME: GCC has trouble with associate to string here
      !associate ( BC => A % BoundaryCondition ( iConnection, iB ) )
      !select case ( trim ( BC ) )
      select case ( trim ( A % BoundaryCondition ( iConnection, iB ) ) )
      case ( 'PERIODIC', 'INFLOW' )

        cycle

      case ( 'OUTFLOW' )

        select type ( C => A % Chart )
        class is ( Chart_SL_Template )
          call Apply_BC_CSL_Outflow ( C, F, iDimension, iConnection )
        class default
          call Show ( 'Chart type not recognized', CONSOLE % ERROR )
          call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
          call Show ( 'ApplyBoundaryConditions', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- C

      case ( 'REFLECTING' )

        select type ( C => A % Chart )
        class is ( Chart_SL_Template )
          call Apply_BC_CSL_Reflecting ( C, F, iDimension, iConnection )
        class default
          call Show ( 'Chart type not recognized', CONSOLE % ERROR )
          call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
          call Show ( 'ApplyBoundaryConditions', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- C

      case ( 'CUSTOM' )

        select type ( C => A % Chart )
        class is ( Chart_SL_Template )
          if ( associated ( A % Apply_BC_CSL_Custom ) ) then
            call A % Apply_BC_CSL_Custom ( C, F, iDimension, iConnection )
          else
            call Show ( 'Custom boundary conditions not set', CONSOLE % ERROR )
            call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
            call Show ( 'ApplyBoundaryConditions', 'subroutine', &
                        CONSOLE % ERROR )
            call PROGRAM_HEADER % Abort ( )
          end if
        class default
          call Show ( 'Chart type not recognized', CONSOLE % ERROR )
          call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
          call Show ( 'ApplyBoundaryConditions', 'subroutine', &
                      CONSOLE % ERROR )
          call PROGRAM_HEADER % Abort ( )
        end select !-- C

      case default
        call Show ( 'BoundaryCondition not recognized', CONSOLE % ERROR )
        call Show ( A % BoundaryCondition ( iConnection, iB ), &
                    'BoundaryCondition', CONSOLE % ERROR )
        call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
        call Show ( 'ApplyBoundaryConditions', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- BC
      ! end associate !-- BC
    end do !-- iB

  end subroutine ApplyBoundaryConditions


  subroutine ApplyBoundaryConditionsFaces ( A, F )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    class ( StorageForm ), intent ( inout ) :: &
      F  !-- Field

    integer ( KDI ) :: &
      iD  !-- iDimension

    do iD = 1, A % nDimensions
      call A % ApplyBoundaryConditions &
             ( F, iD, A % Connectivity % iaInner ( iD ) )
      call A % ApplyBoundaryConditions &
             ( F, iD, A % Connectivity % iaOuter ( iD ) )
    end do

  end subroutine ApplyBoundaryConditionsFaces


  impure elemental subroutine FinalizeTemplate ( A )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A

    if ( allocated ( A % Field ) ) &
      deallocate ( A % Field )
    if ( allocated ( A % Chart ) ) &
      deallocate ( A % Chart )

  end subroutine FinalizeTemplate


  subroutine Apply_BC_CSL_Outflow ( CSL, F, iDimension, iConnection )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      CSL
    class ( StorageForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField

    do iS = 1, F % nVariables
      iF = F % iaSelected ( iS )
      call CSL % CopyBoundary ( F, iF, iDimension, iConnection )
    end do !-- iS

  end subroutine Apply_BC_CSL_Outflow


  subroutine Apply_BC_CSL_Reflecting ( CSL, F, iDimension, iConnection )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      CSL
    class ( StorageForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF, &  !-- iField
      iV     !-- iVector

    do iS = 1, F % nVariables
      iF = F % iaSelected ( iS )
      call CSL % CopyBoundary ( F, iF, iDimension, iConnection )
      do iV = 1, F % nVectors
        if ( iF == F % VectorIndices ( iV ) % Value ( iDimension ) ) &
          call CSL % ReverseBoundary ( F, iF, iDimension, iConnection )
      end do !-- iV
    end do !-- iS

  end subroutine Apply_BC_CSL_Reflecting


end module Atlas_SC__Template
