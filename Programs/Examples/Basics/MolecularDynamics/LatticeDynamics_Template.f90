module LatticeDynamics_Template

  use Basics
  use ParticleDynamics_Template
  use LatticeParticles_Form

  implicit none
  private

  type, public, extends ( ParticleDynamicsTemplate ), abstract :: &
    LatticeDynamicsTemplate
      integer ( KDI ) :: &
        nEquilibrateCycles, &
        TargetTemperatureInterval
      type ( MeasuredValueForm ) :: &
        EnergyUnit
    type ( StorageForm ) :: &
      Correlation
    type ( GridImageStreamForm ) :: &
      CorrelationImageStream
    type ( CurveImageForm ) :: &
      CorrelationFunction
  contains
    procedure, private, pass :: &
      Initialize_LD
    generic, public :: &
      Initialize => Initialize_LD
    procedure, public, pass :: &
      RecordObservables
    procedure, public, pass :: &
      SetExtensiveIntensive
  end type LatticeDynamicsTemplate

    private :: &
      WriteCorrelationFunction

contains


  subroutine Initialize_LD ( LD, EnergyUnit, TimeUnit )

    class ( LatticeDynamicsTemplate ), intent ( inout ) :: &
      LD
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyUnit, &
      TimeUnit

    allocate ( LatticeParticlesForm :: LD % DistributedParticles )
    select type ( LP => LD % DistributedParticles )
    type is ( LatticeParticlesForm )
    call LP % Initialize &
           ( PROGRAM_HEADER % Communicator, TimeUnitOption = TimeUnit )
    end select !-- LP

    call LD % Initialize ( TimeUnit = TimeUnit )

    LD % nEquilibrateCycles = LD % nCycles / 3
    call PROGRAM_HEADER % GetParameter &
           ( LD % nEquilibrateCycles, 'nEquilibrateCycles' )

    LD % TargetTemperatureInterval = LD % nEquilibrateCycles / 10
    call PROGRAM_HEADER % GetParameter &
           ( LD % TargetTemperatureInterval, 'TargetTemperatureInterval' )

    LD % EnergyUnit = EnergyUnit

  end subroutine Initialize_LD


  subroutine RecordObservables ( PD, iCycle )

    class ( LatticeDynamicsTemplate ), intent ( inout ) :: &
      PD
    integer ( KDI ), intent ( in ) :: &
      iCycle

    integer ( KDI ) :: &
      iVrbl, & !-- iVariable
      iB  !-- iBin
    real ( KDR ) :: &
      TimeScale_KE, &
      TimeScale_PE, &
      T_Average, &
      PE_Average, &
      V_Average
    real ( KDR ), dimension ( 3 ) :: &
      CenterOfMass
    type ( CollectiveOperation_R_Form ) :: &
      CO, &
      CO_Correlation

    select type ( LP => PD % DistributedParticles )
    type is ( LatticeParticlesForm )

    associate ( iVl => iCycle + 1 ) !-- iValue offset, iCycle starts at 0

    !-- Extensive

    !-- FIXME: At least with NAG, seems to be necessary to spell out indices of
    !          MOMENTUM and ANGULAR_MOMENTUM, at least in the associate statement
    associate &
      ( KE => PD % Extensive % Value ( iVl, LP % KINETIC_ENERGY ), &
        PE => PD % Extensive % Value ( iVl, LP % POTENTIAL_ENERGY ), &
        TE => PD % Extensive % Value ( iVl, LP % TOTAL_ENERGY ), &
        V  => PD % Extensive % Value ( iVl, LP % VIRIAL ), &
        MM => PD % Extensive % Value &
                 ( iVl, LP % MASS_MOMENT ( 1 ) : LP % MASS_MOMENT ( 3 ) ), &
        P  => PD % Extensive % Value &
                 ( iVl, LP % MOMENTUM ( 1 ) : LP % MOMENTUM ( 3 ) ), &
        L  => PD % Extensive % Value &
                 ( iVl, &
                   LP % ANGULAR_MOMENTUM ( 1 ) : LP % ANGULAR_MOMENTUM ( 3 ) ), &
        SD => PD % Extensive % Value ( iVl, LP % SQUARE_DISPLACEMENT ) )

    !-- CenterOfMass used to compute angular momentum lags by 1 step
    if ( iVl == 1 ) then
      CenterOfMass = 0.0_KDR
    else
      CenterOfMass = PD % Intensive % Value ( iVl - 1, LP % CENTER_OF_MASS )
    end if

    call CO % Initialize &
           ( LP % Communicator, [ LP % N_EXTENSIVE ], [ LP % N_EXTENSIVE ] )
    CO % Outgoing % Value ( LP % KINETIC_ENERGY ) = LP % MyKineticEnergy ( )
    CO % Outgoing % Value ( LP % POTENTIAL_ENERGY ) &
      = LP % MyPotentialEnergy ( PD % PotentialValue )
    CO % Outgoing % Value ( LP % TOTAL_ENERGY ) = 0.0_KDR
    CO % Outgoing % Value ( LP % VIRIAL ) &
      = LP % MyVirial ( PD % VirialValue )
    CO % Outgoing % Value ( LP % MASS_MOMENT ) = LP % MyMassMoment ( )
    CO % Outgoing % Value ( LP % MOMENTUM ) = LP % MyMomentum ( )
    CO % Outgoing % Value ( LP % ANGULAR_MOMENTUM ) &
      = LP % MyAngularMomentum ( CenterOfMass )
    CO % Outgoing % Value ( LP % SQUARE_DISPLACEMENT ) &
      = LP % MySquareDisplacement ( )
    call CO % Reduce ( REDUCTION % SUM )

    KE = CO % Incoming % Value ( LP % KINETIC_ENERGY )
    PE = CO % Incoming % Value ( LP % POTENTIAL_ENERGY )
    TE = KE + PE
     V = CO % Incoming % Value ( LP % VIRIAL )
    MM = CO % Incoming % Value ( LP % MASS_MOMENT )
     P = CO % Incoming % Value ( LP % MOMENTUM )
     L = CO % Incoming % Value ( LP % ANGULAR_MOMENTUM )
    SD = CO % Incoming % Value ( LP % SQUARE_DISPLACEMENT )

    if ( mod ( iCycle, PD % WriteCycleInterval ) == 0 ) then
      call Show ( 'Extensive variables', CONSOLE % INFO_3 )
      do iVrbl = 1, LP % N_EXTENSIVE
        call Show ( PD % Extensive % Value ( iVl, iVrbl ), &
                    PD % Extensive % Unit ( iVrbl ), &
                    PD % Extensive % Variable ( iVrbl ), CONSOLE % INFO_3 )
      end do
    end if

    !-- Intensive

    associate &
      ( CM  => PD % Intensive % Value &
                 ( iVl, &
                   LP % CENTER_OF_MASS ( 1 ) : LP % CENTER_OF_MASS ( 3 ) ), &
        T   => PD % Intensive % Value ( iVl, LP % TEMPERATURE ), &
        CF  => PD % Intensive % Value ( iVl, LP % COMPRESSIBILITY_FACTOR ), &
        MSD => PD % Intensive % Value ( iVl, LP % MEAN_SQUARE_DISPLACEMENT ) )

    CM  = MM / ( LP % ParticleMass  *  LP % nParticles )
    T   = LP % TemperatureExpression ( KE )
    CF  = LP % CompressibilityFactorExpression ( V, T )
    MSD = SD / LP % nParticles

    if ( mod ( iCycle, PD % WriteCycleInterval ) == 0 ) then
      call Show ( 'Intensive variables', CONSOLE % INFO_3 )
      do iVrbl = 1, LP % N_INTENSIVE
        call Show ( PD % Intensive % Value ( iVl, iVrbl ), &
                    PD % Intensive % Unit ( iVrbl ), &
                    PD % Intensive % Variable ( iVrbl ), CONSOLE % INFO_3 )
      end do
    end if

    !-- TimeStep

!     associate &
!       ( M => LP % ParticleMass * LP % nParticles, &
!         R => LP % BoxLength )
!     TimeScale_KE &
!       = PD % TimeScaleExpression &
!           ( EnergyScale = KE, LengthScale = R, MassScale = M )
!     TimeScale_PE &
!       = PD % TimeScaleExpression &
!           ( EnergyScale = abs ( PE ), LengthScale = R, MassScale = M )
!     call Show ( 'Time scales', CONSOLE % INFO_3 )
!     call Show ( TimeScale_KE, PD % TimeUnit, 'TimeScale_KE', CONSOLE % INFO_3 )
!     call Show ( TimeScale_PE, PD % TimeUnit, 'TimeScale_PE', CONSOLE % INFO_3 )
!     PD % TimeStep = PD % TimeStepFactor * min ( TimeScale_KE, TimeScale_PE )
!     end associate !-- M, R

    !-- Equilibration

    select type ( PD )
    class is ( LatticeDynamicsTemplate )
    associate ( TTI => PD % TargetTemperatureInterval )
    if ( iCycle > 0 .and. iCycle <= PD % nEquilibrateCycles &
         .and. mod ( iCycle, TTI ) == 0 ) &
    then

      T_Average &
        = sum ( PD % Intensive % Value ( iVl - ( TTI - 1 ) : iVl, &
                LP % TEMPERATURE ) ) &
           / TTI
      call LP % EnforceTargetTemperature ( T_Average )

    end if
    end associate !-- TTI
    end select !-- PD

    !-- Correlation function

    if ( iCycle == 0 ) call Clear ( LP % PairCount )

    if ( iCycle > PD % nEquilibrateCycles ) then
      call CO_Correlation % Initialize &
             ( LP % Communicator, [ LP % nCorrelationBins ], &
               [ LP % nCorrelationBins ] )
      CO_Correlation % Outgoing % Value = LP % MyPairCount
      call CO_Correlation % Reduce ( REDUCTION % SUM )
      LP % PairCount &
        = LP % PairCount + 0.5_KDR * CO_Correlation % Incoming % Value
    end if

    !-- Thermodynamic averages

    if ( iCycle == PD % nCycles ) then

      associate ( nAverage => PD % nCycles - PD % nEquilibrateCycles )

      T_Average &
        = sum ( PD % Intensive % Value ( iVl - ( nAverage - 1 ) : iVl, &
                LP % TEMPERATURE ) ) &
          / nAverage

      PE_Average &
        = sum ( PD % Extensive % Value ( iVl - ( nAverage - 1 ) : iVl, &
                LP % POTENTIAL_ENERGY ) ) &
          / nAverage

      V_Average &
        = sum ( PD % Extensive % Value ( iVl - ( nAverage - 1 ) : iVl, &
                LP % VIRIAL ) ) &
          / nAverage
      
      LP % PairCount = LP % PairCount / nAverage

      end associate !-- nAverage

      call PD % Correlation % Initialize &
             ( [ LP % nCorrelationBins, 1 ], &
               VariableOption = [ 'CorrelationFunction            ' ], &
               NameOption = 'CorrelationFunction' )

      do iB = 1, LP % nCorrelationBins
        associate &
          ( r => 0.5_KDR * ( LP % CorrelationBinEdge ( iB ) &
                             + LP % CorrelationBinEdge ( iB + 1 ) ), &
            dr => LP % CorrelationBinEdge ( iB + 1 ) &
                  - LP % CorrelationBinEdge ( iB ), &
            Volume => LP % BoxLength ** 3, &
            N      => LP % nParticles, &
            FourPi => 4.0_KDR * CONSTANT % PI )

        PD % Correlation % Value ( iB, 1 ) &
          = ( 2.0_KDR * Volume / ( N * ( N - 1.0_KDR ) ) ) &
            * LP % PairCount ( iB ) / ( FourPi * ( r ** 2 ) * dr )

        end associate !-- r
      end do

      call WriteCorrelationFunction ( PD )

      call Show ( '*** THERMODYNAMIC AVERAGES', CONSOLE % INFO_1 )
      call Show ( T_Average, LP % TemperatureUnit, 'TEMPERATURE', &
                  CONSOLE % INFO_1 )
      call Show ( PE_Average / LP % nParticles, PD % EnergyUnit, &
                  'POTENTIAL ENERGY (per particle)', CONSOLE % INFO_1 )
      call Show ( LP % CompressibilityFactorExpression &
                    ( V_Average, T_Average ), &
                  'COMPRESSIBILITY FACTOR', CONSOLE % INFO_1 )
    end if

    end associate !-- CM, etc.
    end associate !-- KE, etc.
    end associate !-- iVl
    end select !-- LP

  end subroutine RecordObservables


  subroutine SetExtensiveIntensive ( PD )

    class ( LatticeDynamicsTemplate ), intent ( inout ) :: &
      PD

    select type ( LP => PD % DistributedParticles )
    type is ( LatticeParticlesForm )

    LP % N_EXTENSIVE         = LP % N_EXTENSIVE_DP + LP % N_EXTENSIVE_LP
    LP % KINETIC_ENERGY      = 1
    LP % POTENTIAL_ENERGY    = 2
    LP % TOTAL_ENERGY        = 3
    LP % VIRIAL              = 4
    LP % MASS_MOMENT         = [  5,  6,  7 ]
    LP % MOMENTUM            = [  8,  9, 10 ]
    LP % ANGULAR_MOMENTUM    = [ 11, 12, 13 ]
    LP % SQUARE_DISPLACEMENT = 14
    call PD % Extensive % Initialize &
           ( [ PD % nCycles + 1, LP % N_EXTENSIVE ], &
             VariableOption = [ 'KineticEnergy                  ', &
                                'PotentialEnergy                ', &
                                'TotalEnergy                    ', &
                                'Virial                         ', &
                                'MassMoment_1                   ', & 
                                'MassMoment_2                   ', & 
                                'MassMoment_3                   ', & 
                                'Momentum_1                     ', & 
                                'Momentum_2                     ', & 
                                'Momentum_3                     ', & 
                                'AngularMomentum_1              ', & 
                                'AngularMomentum_2              ', & 
                                'AngularMomentum_3              ', &
                                'SquareDisplacement             ' ], & 
             NameOption = 'Extensive', &
             UnitOption &
               = [ PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   LP % MassUnit * LP % LengthUnit, &
                   LP % MassUnit * LP % LengthUnit, &
                   LP % MassUnit * LP % LengthUnit, &
                   LP % MassUnit * LP % LengthUnit / PD % TimeUnit, &
                   LP % MassUnit * LP % LengthUnit / PD % TimeUnit, &
                   LP % MassUnit * LP % LengthUnit / PD % TimeUnit, &
                   LP % MassUnit * LP % LengthUnit ** 2  /  PD % TimeUnit, &
                   LP % MassUnit * LP % LengthUnit ** 2  /  PD % TimeUnit, &
                   LP % MassUnit * LP % LengthUnit ** 2  /  PD % TimeUnit, &
                   LP % LengthUnit ** 2 ] )

    LP % N_INTENSIVE = LP % N_INTENSIVE_DP + LP % N_INTENSIVE_LP
    LP % CENTER_OF_MASS           = [ 1, 2, 3 ]
    LP % TEMPERATURE              = 4
    LP % COMPRESSIBILITY_FACTOR   = 5
    LP % MEAN_SQUARE_DISPLACEMENT = 6
    call PD % Intensive % Initialize &
           ( [ PD % nCycles + 1, LP % N_INTENSIVE ], &
             VariableOption = [ 'CenterOfMass_1                 ', &
                                'CenterOfMass_2                 ', &
                                'CenterOfMass_3                 ', &
                                'Temperature                    ', &
                                'CompressibilityFactor          ', &
                                'MeanSquareDisplacement         ' ], &
             NameOption = 'Intensive', &
             UnitOption &
             = [ LP % LengthUnit, LP % LengthUnit, LP % LengthUnit, &
                 LP % TemperatureUnit, UNIT % IDENTITY, LP % LengthUnit ** 2 ] )

    end select !-- LP

  end subroutine SetExtensiveIntensive


  subroutine WriteCorrelationFunction ( LD )

    class ( LatticeDynamicsTemplate ), intent ( inout ) :: &
      LD

    character ( LDF ) :: &
      OutputDirectory

    if ( PROGRAM_HEADER % Communicator % Rank /= CONSOLE % DisplayRank ) return

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )    

    associate ( GIS => LD % CorrelationImageStream )
    call GIS % Initialize &
           ( trim ( PROGRAM_HEADER % Name ) // '_CorrelationFunction', &
             WorkingDirectoryOption = OutputDirectory )
    call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )

    associate &
      ( CF => LD % CorrelationFunction, &
        LP => LD % DistributedParticles )
    call CF % Initialize ( GIS ) 
    call CF % SetGrid  &
           ( Directory = 'CorrelationFunction', &
             NodeCoordinate = LP % CorrelationBinEdge, &
             nProperCells = LP % nCorrelationBins, oValue = 0, &
             CoordinateUnitOption = LP % LengthUnit, &
             CoordinateLabelOption = 'r' )
    call CF % AddStorage ( LD % Correlation )
    call CF % Write ( )
    end associate !-- CF

    call GIS % Close ( )
    end associate !-- GIS

  end subroutine WriteCorrelationFunction


end module LatticeDynamics_Template
