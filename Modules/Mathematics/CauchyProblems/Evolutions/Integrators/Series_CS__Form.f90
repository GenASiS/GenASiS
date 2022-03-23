module Series_CS__Form

  !-- Series_CurrentSet__Form

  use Basics
  use Fields
  use Series_B__Form

  implicit none
  private

  type, public, extends ( Series_B_Form ) :: Series_CS_Form
    type ( StorageForm ), allocatable :: &
      Interior, &
      Boundary, &  
      Total, &
      Change
    class ( Tally_CS_Form ), pointer :: &
      TallyInterior => null ( ), &
      TallyBoundary => null ( ), &
      TallyTotal    => null ( ), &
      TallyChange   => null ( )
  contains
    procedure, private, pass :: &
      Initialize_CS
    generic, public :: &
      Initialize => Initialize_CS
    procedure, public, pass :: &
      Record
!     procedure, public, pass :: &
!       Restore
    final :: &
      Finalize
  end type Series_CS_Form


contains


  subroutine Initialize_CS &
               ( S, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    character ( * ), dimension ( : ), intent ( in ) :: &
      dT_Label
    type ( MeasuredValueForm ), intent ( in ) :: &
      Unit_T
    real ( KDR ), dimension ( : ), intent ( in ), target :: &
      dT_Candidate
    real ( KDR ), intent ( in ), target :: &
      T
    integer ( KDI ), intent ( in ) :: &
      CommunicatorRank, &
      nWrite
    integer ( KDI ), intent ( in ), target :: &
      iCycle

    integer ( KDI ) :: &
      iS  !-- iSelected
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_CS' 

    call S % Initialize &
           ( GIS, dT_Label, Unit_T, dT_Candidate, T, CommunicatorRank, &
             nWrite, iCycle  )

    S % TallyInterior  =>  CS % TallyInterior
    S % TallyBoundary  =>  CS % TallyBoundary ( 1 ) % Element
    S % TallyTotal     =>  CS % TallyTotal
    S % TallyChange    =>  CS % TallyChange

    associate &
      ( TT  => S % TallyTotal, &
        iaS => S % TallyTotal % iaSelected )
    allocate ( SeriesName ( TT % nSelected ) )
    allocate ( SeriesUnit ( TT % nSelected ) )
    do iS  =  1,  TT % nSelected
      SeriesName ( iS )  =  TT % Variable ( iaS ( iS ) )
      SeriesUnit ( iS )  =  TT % Unit ( iaS ( iS ) )
    end do !-- iS

    allocate ( S % Interior )
    allocate ( S % Boundary )
    allocate ( S % Total )
    allocate ( S % Change )
    associate &
      ( I   =>  S % Interior, &
        By  =>  S % Boundary, &
        T   =>  S % Total, &
        C   =>  S % Change, &
        Bc  =>  S % Basic )
    call I % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption =  trim ( CS % Name ) // '_Interior', &
             ClearOption = .true. )
    call By % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS % Name ) // '_Boundary', &
             ClearOption = .true. )
    call T % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS % Name ) // '_Total', &
             ClearOption = .true. )
    call C % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS % Name ) // '_Change', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( I )
      call CI % AddStorage ( By )
      call CI % AddStorage ( T )
      call CI % AddStorage ( C )
      end associate !-- CI
    end if
    end associate !-- I, etc.

    end associate !-- TT, etc.

  end subroutine Initialize_CS


  subroutine Record ( S )

    class ( Series_CS_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSelected

    call S % Series_B_Form % Record ( )

    associate &
      ( IV   =>  S % Interior % Value, &
        BV   =>  S % Boundary % Value, &
        TV   =>  S % Total % Value, &
        CV   =>  S % Change % Value, &
        iT   =>  S % iTime, &
        TIV  =>  S % TallyInterior % Value, &
        TBV  =>  S % TallyBoundary % Value, &
        TTV  =>  S % TallyTotal % Value, &
        TCV  =>  S % TallyChange % Value, &
        nS   =>  S % TallyTotal % nSelected, &
        iaS  =>  S % TallyTotal % iaSelected )
    do iS  =  1,  nS
      IV ( iT, iS )  =  TIV ( iaS ( iS ) )
      BV ( iT, iS )  =  TBV ( iaS ( iS ) )
      TV ( iT, iS )  =  TTV ( iaS ( iS ) )
      CV ( iT, iS )  =  TCV ( iaS ( iS ) )
    end do !-- iS
    end associate !-- TV, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( S )

    type ( Series_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Change ) ) &
      deallocate ( S % Change )
    if ( allocated ( S % Total ) ) &
      deallocate ( S % Total )
    if ( allocated ( S % Boundary ) ) &
      deallocate ( S % Boundary )
    if ( allocated ( S % Interior ) ) &
      deallocate ( S % Interior )

    nullify ( S % TallyChange )
    nullify ( S % TallyTotal )
    nullify ( S % TallyBoundary )
    nullify ( S % TallyInterior )

  end subroutine Finalize


end module Series_CS__Form
