module Series_R_CC_C__Form

  !-- Series_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Measures_R_CC_C__Form

  implicit none
  private

  type, public, extends ( Series_CS_1D_C_CS_Form ) :: Series_R_CC_C_Form
    type ( StorageForm ), allocatable :: &
      Measures_CC
    class ( Measures_R_CC_C_Form ), pointer :: &
      Measures => null ( )
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    procedure, public, pass :: &
      Record
!     procedure, public, pass :: &
!       Restore
    final :: &
      Finalize
  end type Series_R_CC_C_Form


contains


  subroutine Initialize_R_CC_C &
               ( S, M, CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
                 CommunicatorRank, nWrite, iCycle )

    class ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S
    class ( Measures_R_CC_C_Form ), intent ( in ), target :: &
      M
    class ( CurrentSetForm ), dimension ( : ), intent ( in ) :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ) :: &
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
      iM  !-- iMeasure
    type ( QuantityForm ), dimension ( : ), allocatable :: &
      SeriesUnit
    character ( LDL ), dimension ( : ), allocatable :: &
      SeriesName

    if ( S % Type == '' ) &
      S % Type = 'a Series_R_CC_C' 

    call S % Series_CS_1D_C_CS_Form % Initialize &
           ( CS_1D, CS, GIS, dT_Label, Unit_T, dT_Candidate, T, &
             CommunicatorRank, nWrite, iCycle  )

    S % Measures  =>  M

    allocate ( SeriesName ( M % nMeasures ) )
    allocate ( SeriesUnit ( M % nMeasures ) )
    do iM  =  1,  M % nMeasures
      SeriesName ( iM )  =  M % Name ( iM )
      SeriesUnit ( iM )  =  M % Unit ( iM )
    end do !-- iS

    allocate ( S % Measures_CC )
    associate &
      ( SM  =>  S % Measures_CC, &
        SB  =>  S % Basic )
    call SM % Initialize &
           ( [ SB % nValues, M % nMeasures ], &
             VariableOption = SeriesName, UnitOption = SeriesUnit, &
             NameOption =  'Measures_CC', &
             ClearOption = .true. )
    if ( allocated ( S % CurveImage ) ) then
      associate ( CI => S % CurveImage )
      call CI % AddStorage ( SM )
      end associate !-- CI
    end if
    end associate !-- SM, etc.

  end subroutine Initialize_R_CC_C


  subroutine Record ( S )

    class ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iM  !-- iMeasure

    call S % Series_CS_1D_C_CS_Form % Record ( )

    associate &
      (  SMV  =>  S % Measures_CC % Value, &
        iR    =>  S % iRecord, &
          MV  =>  S % Measures % Value, &
        nM    =>  S % Measures % nMeasures )
    do iM  =  1,  nM
      SMV ( iR, iM )  =  MV ( iM )
    end do !-- iM
    end associate !-- SMV, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( S )

    type ( Series_R_CC_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Measures_CC ) ) &
      deallocate ( S % Measures_CC )

    nullify ( S % Measures )

  end subroutine Finalize


end module Series_R_CC_C__Form
