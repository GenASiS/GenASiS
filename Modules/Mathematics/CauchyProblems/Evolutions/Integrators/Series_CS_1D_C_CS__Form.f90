module Series_CS_1D_C_CS__Form

  !-- Series_CurrentSet_1D_Collected_CurrentSet__Form

  use Basics
  use Fields
  use Series_CS__Form

  implicit none
  private

  type, public, extends ( Series_CS_Form ) :: Series_CS_1D_C_CS_Form
    integer ( KDI ) :: &
      nCurrentSets_1D
    type ( StorageForm ), dimension ( : ), allocatable :: &
      Interior_1D, &
      Boundary_1D, &  
      Total_1D, &
      Change_1D
    class ( Tally_CS_Pointer ), dimension ( : ), allocatable :: &
      TallyInterior_1D, &
      TallyBoundary_1D, &
      TallyTotal_1D, &
      TallyChange_1D
  contains
    procedure, private, pass :: &
      Initialize_CS_1D_C_CS
    generic, public :: &
      Initialize => Initialize_CS_1D_C_CS
    procedure, public, pass :: &
      Record
    procedure, public, pass :: &
      Restore
    final :: &
      Finalize
  end type Series_CS_1D_C_CS_Form


contains


  subroutine Initialize_CS_1D_C_CS &
               ( S, CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_CS_1D_C_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), dimension ( : ), intent ( in ), target :: &
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
      iCS, &  !-- iCurrentSet
      iS      !-- iSelected
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_CS_1D_C_CS' 

    call S % Series_CS_Form % Initialize &
           ( CS, GIS, dT_Label, Unit_T, dT_Candidate, T, CommunicatorRank, &
             nWrite, iCycle )

    S % nCurrentSets_1D  =  size ( CS_1D )

    allocate ( S % TallyInterior_1D ( S % nCurrentSets_1D ) )
    allocate ( S % TallyBoundary_1D ( S % nCurrentSets_1D ) )
    allocate ( S % TallyTotal_1D    ( S % nCurrentSets_1D ) )
    allocate ( S % TallyChange_1D   ( S % nCurrentSets_1D ) )
    allocate ( S % Interior_1D      ( S % nCurrentSets_1D ) )
    allocate ( S % Boundary_1D      ( S % nCurrentSets_1D ) )
    allocate ( S % Total_1D         ( S % nCurrentSets_1D ) )
    allocate ( S % Change_1D        ( S % nCurrentSets_1D ) )

    do iCS  =  1,  S % nCurrentSets_1D

      S % TallyInterior_1D ( iCS ) % Pointer  &
        =>  CS_1D ( iCS ) % TallyInterior
      S % TallyBoundary_1D ( iCS ) % Pointer  &
        =>  CS_1D ( iCS ) % TallyBoundary ( 1 ) % Element
      S % TallyTotal_1D    ( iCS ) % Pointer  &
        =>  CS_1D ( iCS ) % TallyTotal
      S % TallyChange_1D   ( iCS ) % Pointer  &
        =>  CS_1D ( iCS ) % TallyChange

      associate &
        ( TT  => CS_1D ( iCS ) % TallyTotal, &
          iaS => CS_1D ( iCS ) % TallyTotal % iaSelected )
      allocate ( SeriesName ( TT % nSelected ) )
      allocate ( SeriesUnit ( TT % nSelected ) )
      do iS  =  1,  TT % nSelected
        SeriesName ( iS )  =  TT % Variable ( iaS ( iS ) )
        SeriesUnit ( iS )  =  TT % Unit ( iaS ( iS ) )
      end do !-- iS

      associate &
        ( I   =>  S % Interior_1D ( iCS ), &
          By  =>  S % Boundary_1D ( iCS ), &
          Tl  =>  S % Total_1D ( iCS ), &
          C   =>  S % Change_1D ( iCS ), &
          Bc  =>  S % Basic )
      call I % Initialize &
             ( [ Bc % nValues, TT % nSelected ], &
               VariableOption = SeriesName, UnitOption = SeriesUnit, &
               NameOption =  trim ( CS_1D ( iCS ) % Name ) // '_Interior', &
               ClearOption = .true. )
      call By % Initialize &
             ( [ Bc % nValues, TT % nSelected ], &
               VariableOption = SeriesName, UnitOption = SeriesUnit, &
               NameOption = trim ( CS_1D ( iCS ) % Name ) // '_Boundary', &
               ClearOption = .true. )
      call Tl % Initialize &
             ( [ Bc % nValues, TT % nSelected ], &
               VariableOption = SeriesName, UnitOption = SeriesUnit, &
               NameOption = trim ( CS_1D ( iCS ) % Name ) // '_Total', &
               ClearOption = .true. )
      call C % Initialize &
             ( [ Bc % nValues, TT % nSelected ], &
               VariableOption = SeriesName, UnitOption = SeriesUnit, &
               NameOption = trim ( CS_1D ( iCS ) % Name ) // '_Change', &
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

    end do !-- iCS

  end subroutine Initialize_CS_1D_C_CS


  subroutine Record ( S )

    class ( Series_CS_1D_C_CS_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iCS, &  !-- iCurrentSet
      iS      !-- iSelected

    call S % Series_CS_Form % Record ( )

    do iCS  =  1,  S % nCurrentSets_1D
      associate &
        ( IV   =>  S % Interior_1D ( iCS ) % Value, &
          BV   =>  S % Boundary_1D ( iCS ) % Value, &
          TV   =>  S % Total_1D    ( iCS ) % Value, &
          CV   =>  S % Change_1D   ( iCS ) % Value, &
          iR   =>  S % iRecord, &
          TIV  =>  S % TallyInterior_1D ( iCS ) % Pointer % Value, &
          TBV  =>  S % TallyBoundary_1D ( iCS ) % Pointer % Value, &
          TTV  =>  S % TallyTotal_1D    ( iCS ) % Pointer % Value, &
          TCV  =>  S % TallyChange_1D   ( iCS ) % Pointer % Value, &
          nS   =>  S % TallyTotal_1D    ( iCS ) % Pointer % nSelected, &
          iaS  =>  S % TallyTotal_1D    ( iCS ) % Pointer % iaSelected )
      do iS  =  1,  nS
        IV ( iR, iS )  =  TIV ( iaS ( iS ) )
        BV ( iR, iS )  =  TBV ( iaS ( iS ) )
        TV ( iR, iS )  =  TTV ( iaS ( iS ) )
        CV ( iR, iS )  =  TCV ( iaS ( iS ) )
      end do !-- iS
      end associate !-- IV, etc.
    end do !-- iCS

  end subroutine Record


  subroutine Restore ( S, iCycleRestart )

    class ( Series_CS_1D_C_CS_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iCycleRestart

    integer ( KDI ) :: &
      iCS, &  !-- iCurrentSet
      iS      !-- iSelected

    call S % Series_CS_Form % Restore ( iCycleRestart )

    do iCS  =  1,  S % nCurrentSets_1D

      associate &
        ( IV   =>  S % Interior_1D ( iCS ) % Value, &
          BV   =>  S % Boundary_1D ( iCS ) % Value, &
          TV   =>  S % Total_1D    ( iCS ) % Value, &
          CV   =>  S % Change_1D   ( iCS ) % Value, &
         iR    =>  S % iRecord, &
          TIV  =>  S % TallyInterior_1D ( iCS ) % Pointer % Value, &
          TBV  =>  S % TallyBoundary_1D ( iCS ) % Pointer % Value, &
          TTV  =>  S % TallyTotal_1D    ( iCS ) % Pointer % Value, &
          TCV  =>  S % TallyChange_1D   ( iCS ) % Pointer % Value, &
          nS   =>  S % TallyTotal_1D    ( iCS ) % Pointer % nSelected, &
          iaS  =>  S % TallyTotal_1D    ( iCS ) % Pointer % iaSelected )
      do iS = 1, nS
        TIV ( iaS ( iS ) )  =  IV ( iR, iS ) 
        TBV ( iaS ( iS ) )  =  BV ( iR, iS )
        TTV ( iaS ( iS ) )  =  TV ( iR, iS )
        TCV ( iaS ( iS ) )  =  CV ( iR, iS )
      end do !-- iS
      end associate !-- IV, etc.

      associate &
        (  I  =>  S % Interior_1D ( iCS ), &
           B  =>  S % Boundary_1D ( iCS ), &
           T  =>  S % Total_1D    ( iCS ), &
           C  =>  S % Change_1D   ( iCS ), &
          TI  =>  S % TallyInterior_1D ( iCS ) % Pointer, &
          TB  =>  S % TallyBoundary_1D ( iCS ) % Pointer, &
          TT  =>  S % TallyTotal_1D    ( iCS ) % Pointer, &
          TC  =>  S % TallyChange_1D   ( iCS ) % Pointer )
      call TI % Show ( I % Name, CONSOLE % INFO_1 )
      call TB % Show ( B % Name, CONSOLE % INFO_1 )
      call TT % Show ( T % Name, CONSOLE % INFO_1 )
      call TC % Show ( C % Name, CONSOLE % INFO_1 )
      end associate !-- I, etc.

     end do !-- iCS

  end subroutine Restore


  impure elemental subroutine Finalize ( S )

    type ( Series_CS_1D_C_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Change_1D ) ) &
      deallocate ( S % Change_1D )
    if ( allocated ( S % Total_1D ) ) &
      deallocate ( S % Total_1D )
    if ( allocated ( S % Boundary_1D ) ) &
      deallocate ( S % Boundary_1D )
    if ( allocated ( S % Interior_1D ) ) &
      deallocate ( S % Interior_1D )

    if ( allocated ( S % TallyChange_1D ) ) &
      deallocate ( S % TallyChange_1D )
    if ( allocated ( S % TallyTotal_1D ) ) &
      deallocate ( S % TallyTotal_1D )
    if ( allocated ( S % TallyBoundary_1D ) ) &
      deallocate ( S % TallyBoundary_1D )
    if ( allocated ( S % TallyInterior_1D ) ) &
      deallocate ( S % TallyInterior_1D )

  end subroutine Finalize


end module Series_CS_1D_C_CS__Form
