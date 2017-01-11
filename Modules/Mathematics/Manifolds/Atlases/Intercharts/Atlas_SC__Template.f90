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
    type ( Field_ASC_Pointer ), dimension ( : ), allocatable :: &
      Field
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
      FinalizeTemplate
  end type Atlas_SC_Template

    private :: &
      Apply_BC_CSL_Outflow, &
      Apply_BC_CSL_Reflecting

contains


  subroutine InitializeBasic &
               ( A, NameSuffix, CommunicatorOption, IncludeFacesOption, &
                 IncludeEdgesOption, nExcisionsOption, iDimensionalityOption )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    character ( * ), intent ( in )  :: &
      NameSuffix
    type ( CommunicatorForm ), intent ( in ), target, optional :: &
      CommunicatorOption
    logical ( KDL ), intent ( in ), optional :: &
      IncludeFacesOption, &
      IncludeEdgesOption
    integer ( KDI ), intent ( in ), optional :: &
      nExcisionsOption, &
      iDimensionalityOption

    logical ( KDL ), dimension ( 3 ) :: &
      IsPeriodic

    if ( .not. associated ( A % Type ) ) then
      allocate ( A % Type )
      A % Type = 'an Atlas_SC' 
    end if

    call A % AtlasHeaderForm % Initialize &
           ( NameSuffix, CommunicatorOption, IncludeFacesOption, &
             IncludeEdgesOption, nExcisionsOption, iDimensionalityOption )

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
               ( A, SpacingOption, CoordinateSystemOption, &
                 CoordinateUnitOption, MinCoordinateOption, &
                 MaxCoordinateOption, RatioOption, ScaleOption, &
                 nCellsOption, nGhostLayersOption, nDimensionsOption, &
                 nEqualOption )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      SpacingOption
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
      call Show ( 'AddField_CSL', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- C

  end subroutine AddField


  subroutine ApplyBoundaryConditions ( A, F, iDimension, iConnection )

    class ( Atlas_SC_Template ), intent ( inout ) :: &
      A
    class ( VariableGroupForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      iB  !-- iBoundary
    integer ( KDI ), dimension ( 2 ) :: &
      iaC

    do iB = 1, A % nBoundaries
      associate ( BC => A % BoundaryCondition ( iConnection, iB ) )
      select case ( trim ( BC ) )
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

      case default
        call Show ( 'BoundaryCondition not recognized', CONSOLE % ERROR )
        call Show ( BC, 'BoundaryCondition', CONSOLE % ERROR )
        call Show ( 'Atlas_SC__Template', 'module', CONSOLE % ERROR )
        call Show ( 'ApplyBoundaryConditions', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- BC
      end associate !-- BC
    end do !-- iB

  end subroutine ApplyBoundaryConditions


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
    class ( VariableGroupForm ), intent ( inout ) :: &
      F  !-- Field
    integer ( KDI ), intent ( in ) :: &
      iDimension, &
      iConnection

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iF     !-- iField

    do iS = 1, F % nVariables
      iF = F % iaSelected ( iS )
      call CSL % CopyBoundary ( F % Value ( :, iF ), iDimension, iConnection )
    end do !-- iS

  end subroutine Apply_BC_CSL_Outflow


  subroutine Apply_BC_CSL_Reflecting ( CSL, F, iDimension, iConnection )

    class ( Chart_SL_Template ), intent ( inout ) :: &
      CSL
    class ( VariableGroupForm ), intent ( inout ) :: &
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
      call CSL % CopyBoundary ( F % Value ( :, iF ), iDimension, iConnection )
      do iV = 1, F % nVectors
        if ( iF == F % VectorIndices ( iV ) % Value ( iDimension ) ) &
          call CSL % ReverseBoundary &
                 ( F % Value ( :, iF ), iDimension, iConnection )
      end do !-- iV
    end do !-- iS

  end subroutine Apply_BC_CSL_Reflecting


end module Atlas_SC__Template
