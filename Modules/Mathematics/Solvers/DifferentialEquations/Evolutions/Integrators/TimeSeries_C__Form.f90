!-- TimeSeries_C adds to a time series integrals of interest in the 
!   evolution of a conserved current.

module TimeSeries_C__Form

  !-- TimeSeries_Current_Form

  use Basics
  use Fields
  use Integrator_Template
  use TimeSeries_Form

  implicit none
  private

  type, public, extends ( TimeSeriesForm ) :: TimeSeries_C_Form
    type ( StorageForm ), allocatable :: &
      SeriesInterior, &
      SeriesBoundary, &  
      SeriesTotal, &
      SeriesChange
    class ( Tally_C_Form ), pointer :: &
      TallyInterior => null ( ), &
      TallyBoundary => null ( ), &
      TallyTotal => null ( ), &
      TallyChange => null ( )
  contains
    procedure, private, pass :: &
      Initialize_C
    generic, public :: &
      Initialize => Initialize_C
    procedure, public, pass :: &
      Record
    final :: &
      Finalize
  end type TimeSeries_C_Form

contains


  subroutine Initialize_C ( TS, I, CA )

    class ( TimeSeries_C_Form ), intent ( inout ) :: &
      TS
    class ( IntegratorTemplate ), intent ( in ) :: &
      I
    class ( Current_ASC_Template ), intent ( in ), target :: &
      CA

    integer ( KDI ) :: &
      iS  !-- iSelected
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( TS % Type == '' ) &
      TS % Type = 'a TimeSeries_C' 

    call TS % Initialize ( I )

    TS % TallyInterior  =>  CA % TallyInterior
    TS % TallyBoundary  =>  CA % TallyBoundaryGlobal ( 1 ) % Element
    TS % TallyTotal     =>  CA % TallyTotal
    TS % TallyChange    =>  CA % TallyChange

    associate &
      ( TT  => TS % TallyTotal, &
        iaS => TS % TallyTotal % iaSelected )
    allocate ( SeriesName ( TT % nSelected ) )
    allocate ( SeriesUnit ( TT % nSelected ) )
    do iS = 1, TT % nSelected
      SeriesName ( iS ) = TT % Variable ( iaS ( iS ) )
      SeriesUnit ( iS ) = TT % Unit ( iaS ( iS ) )
    end do !-- iS

    allocate ( TS % SeriesInterior )
    allocate ( TS % SeriesBoundary )
    allocate ( TS % SeriesTotal )
    allocate ( TS % SeriesChange )
    associate &
      ( SI => TS % SeriesInterior, &
        SB => TS % SeriesBoundary, &
        ST => TS % SeriesTotal, &
        SC => TS % SeriesChange, &
        SBsc => TS % SeriesBasic )
    call SI % Initialize &
           ( [ SBsc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = 'Interior_' // CA % Name, ClearOption = .true. )
    call SB % Initialize &
           ( [ SBsc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = 'Boundary_' // CA % Name, ClearOption = .true. )
    call ST % Initialize &
           ( [ SBsc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = 'Total_' // CA % Name, ClearOption = .true. )
    call SC % Initialize &
           ( [ SBsc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = 'Change_' // CA % Name, ClearOption = .true. )
    if ( allocated ( TS % CurveImage ) ) then
      associate ( CI => TS % CurveImage )
      call CI % AddStorage ( SI )
      call CI % AddStorage ( SB )
      call CI % AddStorage ( ST )
      call CI % AddStorage ( SC )
      end associate !-- CI
    end if
    end associate !-- SI, etc.

    end associate !-- TT, etc.

  end subroutine Initialize_C


  subroutine Record ( TS, MaxTime, MinTime, MeanTime )

    class ( TimeSeries_C_Form ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iS  !-- iSelected

    call TS % TimeSeriesForm % Record ( MaxTime, MinTime, MeanTime )

    associate &
      ( SIV => TS % SeriesInterior % Value, &
        SBV => TS % SeriesBoundary % Value, &
        STV => TS % SeriesTotal % Value, &
        SCV => TS % SeriesChange % Value, &
        iV  => TS % iTime, &
        TIV => TS % TallyInterior % Value, &
        TBV => TS % TallyBoundary % Value, &
        TTV => TS % TallyTotal % Value, &
        TCV => TS % TallyChange % Value, &
        nS  => TS % TallyTotal % nSelected, &
        iaS => TS % TallyTotal % iaSelected )
    do iS = 1, nS
      SIV ( iV, iS ) = TIV ( iaS ( iS ) )
      SBV ( iV, iS ) = TBV ( iaS ( iS ) )
      STV ( iV, iS ) = TTV ( iaS ( iS ) )
      SCV ( iV, iS ) = TCV ( iaS ( iS ) )
    end do !-- iS
    end associate !-- STV, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( TS )

    type ( TimeSeries_C_Form ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % SeriesChange ) ) &
      deallocate ( TS % SeriesChange )
    if ( allocated ( TS % SeriesTotal ) ) &
      deallocate ( TS % SeriesTotal )
    if ( allocated ( TS % SeriesBoundary ) ) &
      deallocate ( TS % SeriesBoundary )
    if ( allocated ( TS % SeriesInterior ) ) &
      deallocate ( TS % SeriesInterior )

    nullify ( TS % TallyChange )
    nullify ( TS % TallyTotal )
    nullify ( TS % TallyBoundary )
    nullify ( TS % TallyInterior )

  end subroutine Finalize


end module TimeSeries_C__Form
