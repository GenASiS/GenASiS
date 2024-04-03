module Series_CS_1D_CS__Form

  !-- Series_CurrentSet_1D_CurrentSet__Form

  use Basics
  use Fields
  use Series_CS__Form

  implicit none
  private

  type, public, extends ( Series_CS_Form ) :: Series_CS_1D_CS_Form
    type ( StorageForm ), allocatable :: &
      Interior_1D, &
      Boundary_1D, &  
      Total_1D, &
      Change_1D
    class ( Tally_CS_Form ), pointer :: &
      TallyInterior_1D => null ( ), &
      TallyBoundary_1D => null ( ), &
      TallyTotal_1D    => null ( ), &
      TallyChange_1D   => null ( )
  contains
    procedure, private, pass :: &
      Initialize_CS_1D_CS
    generic, public :: &
      Initialize => Initialize_CS_1D_CS
    procedure, public, pass :: &
      Record
    procedure, public, pass :: &
      Restore
    final :: &
      Finalize
  end type Series_CS_1D_CS_Form


contains


  subroutine Initialize_CS_1D_CS &
               ( S, CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_CS_1D_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    type ( GridImageStreamForm ), intent ( in ) :: &
      GIS
    character ( * ), dimension ( : ), intent ( in ) :: &
      dT_Label
    type ( QuantityForm ), intent ( in ) :: &
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
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_CS_1D_CS' 

    call S % Series_CS_Form % Initialize &
           ( CS, GIS, dT_Label, Unit_T, dT_Candidate, T, CommunicatorRank, &
             nWrite, iCycle )

    S % TallyInterior_1D  =>  CS_1D % TallyInterior
    S % TallyBoundary_1D  =>  CS_1D % TallyBoundary ( 1 ) % Element
    S % TallyTotal_1D     =>  CS_1D % TallyTotal
    S % TallyChange_1D    =>  CS_1D % TallyChange

    associate &
      ( TT  => S % TallyTotal_1D, &
        iaS => S % TallyTotal_1D % iaSelected )
    allocate ( SeriesName ( TT % nSelected ) )
    allocate ( SeriesUnit ( TT % nSelected ) )
    do iS  =  1,  TT % nSelected
      SeriesName ( iS )  =  TT % Variable ( iaS ( iS ) )
      SeriesUnit ( iS )  =  TT % Unit ( iaS ( iS ) )
    end do !-- iS

    allocate ( S % Interior_1D )
    allocate ( S % Boundary_1D )
    allocate ( S % Total_1D )
    allocate ( S % Change_1D )
    associate &
      ( I   =>  S % Interior_1D, &
        By  =>  S % Boundary_1D, &
        Tl  =>  S % Total_1D, &
        C   =>  S % Change_1D, &
        Bc  =>  S % Basic )
    call I % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption =  trim ( CS_1D % Name ) // '_Interior', &
             ClearOption = .true. )
    call By % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS_1D % Name ) // '_Boundary', &
             ClearOption = .true. )
    call Tl % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS_1D % Name ) // '_Total', &
             ClearOption = .true. )
    call C % Initialize &
           ( [ Bc % nValues, TT % nSelected ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption = trim ( CS_1D % Name ) // '_Change', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( I )
      call CI % AddStorage ( By )
      call CI % AddStorage ( Tl )
      call CI % AddStorage ( C )
      end associate !-- CI
    end if
    end associate !-- I, etc.

    end associate !-- TT, etc.

  end subroutine Initialize_CS_1D_CS


  subroutine Record ( S )

    class ( Series_CS_1D_CS_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iSelected

    call S % Series_CS_Form % Record ( )

    associate &
      ( IV   =>  S % Interior_1D % Value, &
        BV   =>  S % Boundary_1D % Value, &
        TV   =>  S % Total_1D % Value, &
        CV   =>  S % Change_1D % Value, &
        iR   =>  S % iRecord, &
        TIV  =>  S % TallyInterior_1D % Value, &
        TBV  =>  S % TallyBoundary_1D % Value, &
        TTV  =>  S % TallyTotal_1D % Value, &
        TCV  =>  S % TallyChange_1D % Value, &
        nS   =>  S % TallyTotal_1D % nSelected, &
        iaS  =>  S % TallyTotal_1D % iaSelected )
    do iS  =  1,  nS
      IV ( iR, iS )  =  TIV ( iaS ( iS ) )
      BV ( iR, iS )  =  TBV ( iaS ( iS ) )
      TV ( iR, iS )  =  TTV ( iaS ( iS ) )
      CV ( iR, iS )  =  TCV ( iaS ( iS ) )
    end do !-- iS
    end associate !-- IV, etc.

  end subroutine Record


  subroutine Restore ( S, iCycleRestart )

    class ( Series_CS_1D_CS_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iCycleRestart

    integer ( KDI ) :: &
      iS  !-- iSelected

    call S % Series_CS_Form % Restore ( iCycleRestart )

    associate &
      ( IV   =>  S % Interior_1D % Value, &
        BV   =>  S % Boundary_1D % Value, &
        TV   =>  S % Total_1D % Value, &
        CV   =>  S % Change_1D % Value, &
       iR    =>  S % iRecord, &
        TIV  =>  S % TallyInterior_1D % Value, &
        TBV  =>  S % TallyBoundary_1D % Value, &
        TTV  =>  S % TallyTotal_1D % Value, &
        TCV  =>  S % TallyChange_1D % Value, &
        nS   =>  S % TallyTotal_1D % nSelected, &
        iaS  =>  S % TallyTotal_1D % iaSelected )
    do iS = 1, nS
      TIV ( iaS ( iS ) )  =  IV ( iR, iS ) 
      TBV ( iaS ( iS ) )  =  BV ( iR, iS )
      TTV ( iaS ( iS ) )  =  TV ( iR, iS )
      TCV ( iaS ( iS ) )  =  CV ( iR, iS )
    end do !-- iS
    end associate !-- IV, etc.

    associate &
      (  I  =>  S % Interior_1D, &
         B  =>  S % Boundary_1D, &
         T  =>  S % Total_1D, &
         C  =>  S % Change_1D, &
        TI  =>  S % TallyInterior_1D, &
        TB  =>  S % TallyBoundary_1D, &
        TT  =>  S % TallyTotal_1D, &
        TC  =>  S % TallyChange_1D )
    call TI % Show ( I % Name, CONSOLE % INFO_1 )
    call TB % Show ( B % Name, CONSOLE % INFO_1 )
    call TT % Show ( T % Name, CONSOLE % INFO_1 )
    call TC % Show ( C % Name, CONSOLE % INFO_1 )
    end associate !-- I, etc.

  end subroutine Restore


  impure elemental subroutine Finalize ( S )

    type ( Series_CS_1D_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Change_1D ) ) &
      deallocate ( S % Change_1D )
    if ( allocated ( S % Total_1D ) ) &
      deallocate ( S % Total_1D )
    if ( allocated ( S % Boundary_1D ) ) &
      deallocate ( S % Boundary_1D )
    if ( allocated ( S % Interior_1D ) ) &
      deallocate ( S % Interior_1D )

    nullify ( S % TallyChange_1D )
    nullify ( S % TallyTotal_1D )
    nullify ( S % TallyBoundary_1D )
    nullify ( S % TallyInterior_1D )

  end subroutine Finalize


end module Series_CS_1D_CS__Form
