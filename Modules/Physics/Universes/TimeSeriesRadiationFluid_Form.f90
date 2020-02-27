module TimeSeriesRadiationFluid_Form

  use Basics
  use Mathematics
  use Universe_Template

  implicit none
  private

  type, public, extends ( TimeSeries_C_1D_C_Form ) :: &
    TimeSeriesRadiationFluidForm
      integer ( KDI ) :: &
        iEnergy_F, &
        iEnergy_R, &
        iMomentum_1_F, iMomentum_2_F, iMomentum_3_F, &
        iMomentum_1_R, iMomentum_2_R, iMomentum_3_R, &
        iNumber_F = 0, &
        iNumber_R = 0, &
        iSpeciesNumberPlus = 0, &
        iSpeciesNumberMinus = 0
      type ( StorageForm ), allocatable :: &
        SeriesChangeRadiationFluid
  contains
    procedure, private, pass :: &
      Initialize_RF
    generic, public :: &
      Initialize => Initialize_RF
    procedure, public, pass :: &
      Record
    final :: &
      Finalize
  end type TimeSeriesRadiationFluidForm

contains


  subroutine Initialize_RF &
               ( TS, U, iSpeciesNumberPlusOption, iSpeciesNumberMinusOption )

    class ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS
    class ( UniverseTemplate ), intent ( in ) :: &
      U
    integer ( KDI ), intent ( in ), optional :: &
      iSpeciesNumberPlusOption, &
      iSpeciesNumberMinusOption

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      nVariables
    type ( MeasuredValueForm ), dimension ( : ), allocatable :: &
      Unit
    logical ( KDL ) :: &
      TotalEnergyFound
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
    if ( present ( iSpeciesNumberPlusOption ) &
         .and. present ( iSpeciesNumberMinusOption ) ) &
    then
      nVariables = 5
      TS % iSpeciesNumberPlus   =  iSpeciesNumberPlusOption
      TS % iSpeciesNumberMinus  =  iSpeciesNumberMinusOption
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

    !-- Fluid indices

    associate ( TC => TS % TallyChange )
    associate ( iaS => TC % iaSelected )

    !-- Check for Fluid TotalEnergy, Momentum, and ElectronNumber
    TotalEnergyFound = .false.
    do iS = 1, TC % nSelected
      if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'TotalEnergy' ) then
        TS % iEnergy_F  =  iaS ( iS )
        TotalEnergyFound  =  .true.
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_1' ) then
        TS % iMomentum_1_F  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_2' ) then
        TS % iMomentum_2_F  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_3' ) then
        TS % iMomentum_3_F  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'ElectronNumber' ) &
      then
        TS % iNumber_F  =  iaS ( iS )
      end if
    end do !-- iS

    !-- If TotalEnergy not there, use FluidEnergy
    if ( .not. TotalEnergyFound ) then
      do iS = 1, TC % nSelected
        if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'FluidEnergy' ) then
          TS % iEnergy_F  =  iaS ( iS )
          exit
        end if
      end do !-- iS
    end if

    end associate !-- iaS
    end associate !-- TC

    !-- Radiation indices

    associate ( TC => TS % TallyChange_1D ( 1 ) % Pointer )
    associate ( iaS => TC % iaSelected )
    do iS = 1, TC % nSelected
      if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Energy' ) then
        TS % iEnergy_R  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_1' ) then
        TS % iMomentum_1_R  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_2' ) then
        TS % iMomentum_2_R  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Momentum_3' ) then
        TS % iMomentum_3_R  =  iaS ( iS )
      else if ( trim ( TC % Variable ( iaS ( iS ) ) )  ==  'Number' ) &
      then
        TS % iNumber_R  =  iaS ( iS )
      end if
    end do !-- iS
    end associate !-- iaS
    end associate

  end subroutine Initialize_RF


  subroutine Record ( TS, MaxTime, MinTime, MeanTime )

    class ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iR  !-- iRadiation

    call TS % TimeSeries_C_1D_C_Form % Record ( MaxTime, MinTime, MeanTime )

    associate &
      ( SCRFV => TS % SeriesChangeRadiationFluid % Value, &
           iV => TS % iTime )
    !-- See SCRFV indices in Initialize_RF above

    !-- Fluid contribution
    associate ( TCV => TS % TallyChange % Value )
    SCRFV ( iV, 1 )  =  TCV ( TS % iEnergy_F )
    SCRFV ( iV, 2 )  =  TCV ( TS % iMomentum_1_F )
    SCRFV ( iV, 3 )  =  TCV ( TS % iMomentum_2_F )
    SCRFV ( iV, 4 )  =  TCV ( TS % iMomentum_3_F )
    if ( TS % iNumber_F > 0 ) &
      SCRFV ( iV, 5 )  =  TCV ( TS % iNumber_F )
    end associate !-- TCV

    !-- Radiation contribution

    do iR = 1, size ( TS % TallyChange_1D )

      associate ( TCV => TS % TallyChange_1D ( iR ) % Pointer % Value )
      SCRFV ( iV, 1 )  =  SCRFV ( iV, 1 )  +  TCV ( TS % iEnergy_R )
      SCRFV ( iV, 2 )  =  SCRFV ( iV, 2 )  +  TCV ( TS % iMomentum_1_R )
      SCRFV ( iV, 3 )  =  SCRFV ( iV, 3 )  +  TCV ( TS % iMomentum_2_R )
      SCRFV ( iV, 4 )  =  SCRFV ( iV, 4 )  +  TCV ( TS % iMomentum_3_R )
      if ( TS % iNumber_R > 0 ) then
        if ( iR == TS % iSpeciesNumberPlus ) &
          SCRFV ( iV, 5 )  =  SCRFV ( iV, 5 )  +  TCV ( TS % iNumber_R )
        if ( iR == TS % iSpeciesNumberMinus ) &
          SCRFV ( iV, 5 )  =  SCRFV ( iV, 5 )  -  TCV ( TS % iNumber_R )
      end if
      end associate !-- TCV

    end do !-- iR

    end associate !-- SCRF, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( TS ) 

    type ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % SeriesChangeRadiationFluid ) ) &
      deallocate ( TS % SeriesChangeRadiationFluid )

  end subroutine Finalize


end module TimeSeriesRadiationFluid_Form
