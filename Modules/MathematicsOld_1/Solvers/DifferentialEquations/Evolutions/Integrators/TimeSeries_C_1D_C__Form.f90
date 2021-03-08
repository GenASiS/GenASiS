!-- TimeSeries_C_1D_C adds to a time series integrals of interest in the 
!   evolution of multiple conserved currents, and an additional conserved 
!   current.

module TimeSeries_C_1D_C__Form

  !-- TimeSeries_Current_1D_Current_Form

  use Basics
  use Fields
  use Integrator_C_PS__Form
  use TimeSeries_C__Form

  implicit none
  private

  type, public, extends ( TimeSeries_C_Form ) :: TimeSeries_C_1D_C_Form
    type ( StorageForm ), dimension ( : ), allocatable :: &
      SeriesInterior_1D, &
      SeriesBoundary_1D, &  
      SeriesTotal_1D, &
      SeriesChange_1D
    class ( Tally_C_PointerForm ), dimension ( : ), allocatable :: &
      TallyInterior_1D, &
      TallyBoundary_1D, &
      TallyTotal_1D, &
      TallyChange_1D
  contains
    procedure, private, pass :: &
      Initialize_C_1D_C
    generic, public :: &
      Initialize => Initialize_C_1D_C
    procedure, public, pass :: &
      Record
    final :: &
      Finalize
  end type TimeSeries_C_1D_C_Form

contains


  subroutine Initialize_C_1D_C ( TS, I, CA_1D )

    class ( TimeSeries_C_1D_C_Form ), intent ( inout ) :: &
      TS
    class ( Integrator_C_PS_Form ), intent ( in ) :: &
      I
    class ( * ), dimension ( : ), intent ( in ), target :: &
      CA_1D

    integer ( KDI ) :: &
      iS,  &  !-- iSelected
      iCA, &  !-- iCurrentAtlas
      nCA
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName
    class ( Current_ASC_Template ), pointer :: &
      CA

    if ( TS % Type == '' ) &
      TS % Type = 'a TimeSeries_C_1D_C' 

    if ( allocated ( I % Current_ASC ) ) then
      associate ( CA => I % Current_ASC )
      call TS % Initialize ( I, CA )
      end associate !-- CA
    else
      call TS % Initialize ( I )
    end if

    nCA  =  size ( CA_1D )

    allocate ( TS % TallyInterior_1D ( nCA ) )
    allocate ( TS % TallyBoundary_1D ( nCA ) )
    allocate ( TS % TallyTotal_1D ( nCA ) )
    allocate ( TS % TallyChange_1D ( nCA ) )
    allocate ( TS % SeriesInterior_1D ( nCA ) )
    allocate ( TS % SeriesBoundary_1D ( nCA ) )
    allocate ( TS % SeriesTotal_1D ( nCA ) )
    allocate ( TS % SeriesChange_1D ( nCA ) )

    do iCA = 1, nCA

      select type ( CA_1D )
      class is ( Current_ASC_ElementForm )
        CA  =>  CA_1D ( iCA ) % Element
      class is ( Current_BSLL_ASC_CSLD_ElementForm )
        select type ( BI => CA_1D ( iCA ) % Element % BundleIntegral )
        class is ( Current_ASC_Template )
          CA  =>  BI
        end select !-- BI
      end select !-- CA_1D

      TS % TallyInterior_1D ( iCA ) % Pointer  &
        =>  CA % TallyInterior 
      TS % TallyBoundary_1D ( iCA ) % Pointer  &
        =>  CA % TallyBoundaryGlobal ( 1 ) % Element 
      TS % TallyTotal_1D ( iCA ) % Pointer  &
        =>  CA % TallyTotal
      TS % TallyChange_1D ( iCA ) % Pointer  &
        =>  CA % TallyChange  

      associate &
        ( TT  => CA % TallyTotal, &
          iaS => CA % TallyTotal % iaSelected )

      allocate ( SeriesName ( TT % nSelected ) )
      allocate ( SeriesUnit ( TT % nSelected ) )
      do iS = 1, TT % nSelected
        SeriesName ( iS ) = TT % Variable ( iaS ( iS ) )
        SeriesUnit ( iS ) = TT % Unit ( iaS ( iS ) )
      end do !-- iS

      associate &
        ( SI => TS % SeriesInterior_1D ( iCA ), &
          SB => TS % SeriesBoundary_1D ( iCA ), &
          ST => TS % SeriesTotal_1D ( iCA ), &
          SC => TS % SeriesChange_1D ( iCA ), &
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

      deallocate ( SeriesUnit )
      deallocate ( SeriesName )

      end associate !-- TT, etc.

    end do !-- iCA

    nullify ( CA )

  end subroutine Initialize_C_1D_C


  subroutine Record ( TS, MaxTime, MinTime, MeanTime )

    class ( TimeSeries_C_1D_C_Form ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iCA

    call TS % TimeSeries_C_Form % Record ( MaxTime, MinTime, MeanTime )

    do iCA = 1, size ( TS % SeriesInterior_1D )
      associate &
        ( SIV => TS % SeriesInterior_1D ( iCA ) % Value, &
          SBV => TS % SeriesBoundary_1D ( iCA ) % Value, &
          STV => TS % SeriesTotal_1D ( iCA ) % Value, &
          SCV => TS % SeriesChange_1D ( iCA ) % Value, &
          iV  => TS % iTime, &
          TIV => TS % TallyInterior_1D ( iCA ) % Pointer % Value, &
          TBV => TS % TallyBoundary_1D ( iCA ) % Pointer % Value, &
          TTV => TS % TallyTotal_1D ( iCA ) % Pointer % Value, &
          TCV => TS % TallyChange_1D ( iCA ) % Pointer % Value, &
          nS  => TS % TallyTotal_1D ( iCA ) % Pointer % nSelected, &
          iaS => TS % TallyTotal_1D ( iCA ) % Pointer % iaSelected )
      do iS = 1, nS
        SIV ( iV, iS ) = TIV ( iaS ( iS ) )
        SBV ( iV, iS ) = TBV ( iaS ( iS ) )
        STV ( iV, iS ) = TTV ( iaS ( iS ) )
        SCV ( iV, iS ) = TCV ( iaS ( iS ) )
      end do !-- iS
      end associate !-- STV, etc.
    end do !-- iCA

  end subroutine Record


  impure elemental subroutine Finalize ( TS )

    type ( TimeSeries_C_1D_C_Form ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % SeriesChange_1D ) ) &
      deallocate ( TS % SeriesChange_1D )
    if ( allocated ( TS % SeriesTotal_1D ) ) &
      deallocate ( TS % SeriesTotal_1D )
    if ( allocated ( TS % SeriesBoundary_1D ) ) &
      deallocate ( TS % SeriesBoundary_1D )
    if ( allocated ( TS % SeriesInterior_1D ) ) &
      deallocate ( TS % SeriesInterior_1D )

    if ( allocated ( TS % TallyChange_1D ) ) &
      deallocate ( TS % TallyChange_1D )
    if ( allocated ( TS % TallyTotal_1D ) ) &
      deallocate ( TS % TallyTotal_1D )
    if ( allocated ( TS % TallyBoundary_1D ) ) &
      deallocate ( TS % TallyBoundary_1D )
    if ( allocated ( TS % TallyInterior_1D ) ) &
      deallocate ( TS % TallyInterior_1D )

  end subroutine Finalize


end module TimeSeries_C_1D_C__Form
