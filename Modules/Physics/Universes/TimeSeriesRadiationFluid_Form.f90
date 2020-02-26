module TimeSeriesRadiationFluid_Form

  use Basics
  use Mathematics
  use Universe_Template

  implicit none
  private

  type, public, extends ( TimeSeries_C_1D_C_Form ) :: &
    TimeSeriesRadiationFluidForm
      type ( StorageForm ), allocatable :: &
        SeriesChangeRadiationFluid
  contains
    procedure, private, pass :: &
      Initialize_RF
    generic, public :: &
      Initialize => Initialize_RF
    final :: &
      Finalize
  end type TimeSeriesRadiationFluidForm

contains


  subroutine Initialize_RF ( TS, U, iNumberPlusOption, iNumberMinusOption )

    class ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS
    class ( UniverseTemplate ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      iNumberPlusOption, &
      iNumberMinusOption

    integer ( KDI ) :: &
      nVariables
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    character ( LDL ), dimension ( : ), allocatable :: &
      Variable

    if ( TS % Type == '' ) &
      TS % Type = 'a TimeSeriesRadiationFluid' 

    select type ( I => U % Integrator )
    class is ( Integrator_C_1D_PS_C_PS_Form )
      call TS % Initialize ( I, I % Current_ASC_1D )
    class is ( Integrator_C_1D_MS_C_PS_Form )
      call TS % Initialize ( I, I % Current_BSLL_ASC_CSLD_1D )
    end select !-- I

    nVariables = 4
    if ( present ( iNumberPlusOption ) .or. present ( iNumberMinusOption ) ) &
    then
      nVariables = 5
    end if

    allocate ( Variable ( nVariables ) )
    allocate ( Unit ( nVariables ) )

    Variable ( : 4 ) &
      =  [ 'Energy    ', 'Momentum_1', 'Momentum_2', 'Momentum_3' ]
    Unit ( : 4 ) &
      =  [ U % Units % Energy, U % Units % Momentum, &
           U % Units % Momentum, U % Units % Momentum ]
    if ( nVariables == 5 ) then
      Variable ( 5 )  =  'Number'
          Unit ( 5 )  =  U % Units % Number    
    end if

    allocate ( TS % SeriesChangeRadiationFluid )
    associate &
      ( SCRF => TS % SeriesChangeRadiationFluid, &
        SBsc => TS % SeriesBasic )
    call SCRF % Initialize &
           ( [ SBsc % nValues, nVariables ], VariableOption = Variable, &
             UnitOption = Unit, NameOption = 'Change_RadiationFluid', &
             ClearOption = .true. )
    if ( allocated ( TS % CurveImage ) ) then
      associate ( CI => TS % CurveImage )
      call CI % AddStorage ( SCRF )
      end associate !-- CI
    end if
    end associate !-- STC, etc.

  end subroutine Initialize_RF


  impure elemental subroutine Finalize ( TS ) 

    type ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % SeriesChangeRadiationFluid ) ) &
      deallocate ( TS % SeriesChangeRadiationFluid )

  end subroutine Finalize


end module TimeSeriesRadiationFluid_Form
