module Series_F_CC__Form

  !-- Series_Fluid_CentralCore__Form

  use Basics
  use Mathematics
  use Measures_F_CC__Form

  implicit none
  private

  type, public, extends ( Series_CS_Form ) :: Series_F_CC_Form
    type ( StorageForm ), allocatable :: &
      Measures_F_CC
    class ( Measures_F_CC_Form ), pointer :: &
      Measures => null ( )
  contains
    procedure, private, pass :: &
      Initialize_F_CC
    generic, public :: &
      Initialize => Initialize_F_CC
    procedure, public, pass :: &
      Record
!     procedure, public, pass :: &
!       Restore
    final :: &
      Finalize
  end type Series_F_CC_Form


contains


  subroutine Initialize_F_CC &
               ( S, M, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_F_CC_Form ), intent ( inout ) :: &
      S
    class ( Measures_F_CC_Form ), intent ( in ), target :: &
      M
    class ( CurrentSetForm ), intent ( in ) :: &
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
      iM  !-- iMeasure
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_F_CC' 

    call S % Series_CS_Form % Initialize &
           ( CS, GIS, dT_Label, Unit_T, dT_Candidate, T, CommunicatorRank, &
             nWrite, iCycle  )

    S % Measures  =>  M

    allocate ( SeriesName ( M % nMeasures ) )
    allocate ( SeriesUnit ( M % nMeasures ) )
    do iM  =  1,  M % nMeasures
      SeriesName ( iM )  =  M % Name ( iM )
      SeriesUnit ( iM )  =  M % Unit ( iM )
    end do !-- iS

    allocate ( S % Measures_F_CC )
    associate &
      ( SM  =>  S % Measures_F_CC, &
        SB  =>  S % Basic )
    call SM % Initialize &
           ( [ SB % nValues, M % nMeasures ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption =  'Measures_F_CC', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( SM )
      end associate !-- CI
    end if
    end associate !-- SM, etc.

  end subroutine Initialize_F_CC


  subroutine Record ( S )

    class ( Series_F_CC_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iM  !-- iMeasure

    call S % Series_CS_Form % Record ( )

    associate &
      (  SMV  =>  S % Measures_F_CC % Value, &
        iR    =>  S % iRecord, &
          MV  =>  S % Measures % Value, &
        nM    =>  S % Measures % nMeasures )
    do iM  =  1,  nM
      SMV ( iR, iM )  =  MV ( iM )
    end do !-- iM
    end associate !-- SMV, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( S )

    type ( Series_F_CC_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Measures_F_CC ) ) &
      deallocate ( S % Measures_F_CC )

    nullify ( S % Measures )

  end subroutine Finalize


end module Series_F_CC__Form
