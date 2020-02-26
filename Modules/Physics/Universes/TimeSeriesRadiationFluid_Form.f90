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
    procedure, public, pass :: &
      Record
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


  subroutine Record ( TS, MaxTime, MinTime, MeanTime )

    class ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      MaxTime, &
      MinTime, &
      MeanTime

    integer ( KDI ) :: &
      iS, &  !-- iSelected
      iEnergy, &
      iMomentum_1, iMomentum_2, iMomentum_3
    logical ( KDL ) :: &
      TotalEnergyFound

    call TS % TimeSeries_C_1D_C_Form % Record ( MaxTime, MinTime, MeanTime )

    associate &
      ( SCRF => TS % SeriesChangeRadiationFluid % Value, &
          iV => TS % iTime )

    !-- Fluid contribution

    associate ( TC_F => TS % TallyChange )
    associate ( iaS => TC_F % iaSelected )

    !-- Check for Fluid TotalEnergy and Momentum
    TotalEnergyFound = .false.
    do iS = 1, TC_F % nSelected
      if ( trim ( TC_F % Variable ( iaS ( iS ) ) )  ==  'TotalEnergy' ) then
        iEnergy  =  iaS ( iS )
        TotalEnergyFound  =  .true.
      else if ( trim ( TC_F % Variable ( iaS ( iS ) ) )  ==  'Momentum_1' ) then
        iMomentum_1  =  iaS ( iS )
      else if ( trim ( TC_F % Variable ( iaS ( iS ) ) )  ==  'Momentum_2' ) then
        iMomentum_2  =  iaS ( iS )
      else if ( trim ( TC_F % Variable ( iaS ( iS ) ) )  ==  'Momentum_3' ) then
        iMomentum_3  =  iaS ( iS )
      end if
    end do !-- iS

    !-- If not there, use FluidEnergy
    if ( .not. TotalEnergyFound ) then
      do iS = 1, TC_F % nSelected
        if ( trim ( TC_F % Variable ( iaS ( iS ) ) )  ==  'FluidEnergy' ) then
          iEnergy  =  iaS ( iS )
          exit
        end if
      end do !-- iS
    end if

    !-- See SCRF indices in Initialize_RF above
    SCRF ( iV, 1 )  =  TC_F % Value ( iEnergy )
    SCRF ( iV, 2 )  =  TC_F % Value ( iMomentum_1 )
    SCRF ( iV, 3 )  =  TC_F % Value ( iMomentum_2 )
    SCRF ( iV, 4 )  =  TC_F % Value ( iMomentum_3 )

call Show ( '>>> Recording' )
call Show ( TC_F % Variable, '>>> Variable_F' )
call Show ( TotalEnergyFound, '>>> TotalEnergyFound' )
call Show ( iEnergy, '>>> iEnergy' )
call Show ( iMomentum_1, '>>> iMomentum_1' )
call Show ( iMomentum_2, '>>> iMomentum_2' )
call Show ( iMomentum_3, '>>> iMomentum_3' )
call Show ( SCRF ( iV, 1 ), TS % SeriesChangeRadiationFluid % Unit ( 1 ), '>>> Energy fluid change' )
call Show ( SCRF ( iV, 2 ), TS % SeriesChangeRadiationFluid % Unit ( 2 ), '>>> Momentum_1 fluid change' )
call Show ( SCRF ( iV, 3 ), TS % SeriesChangeRadiationFluid % Unit ( 3 ), '>>> Momentum_2 fluid change' )
call Show ( SCRF ( iV, 4 ), TS % SeriesChangeRadiationFluid % Unit ( 4 ), '>>> Momentum_3 fluid change' )

    end associate !-- iaS
    end associate !-- TC_F
    end associate !-- SCRF, etc.

  end subroutine Record


  impure elemental subroutine Finalize ( TS ) 

    type ( TimeSeriesRadiationFluidForm ), intent ( inout ) :: &
      TS

    if ( allocated ( TS % SeriesChangeRadiationFluid ) ) &
      deallocate ( TS % SeriesChangeRadiationFluid )

  end subroutine Finalize


end module TimeSeriesRadiationFluid_Form
