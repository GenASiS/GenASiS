module ClusterDynamics_Template

  use Basics
  use ParticleDynamics_Template
  use ClusterParticles_Form

  implicit none
  private

  type, public, extends ( ParticleDynamicsTemplate ), abstract :: &
    ClusterDynamicsTemplate
      type ( MeasuredValueForm ) :: &
        EnergyUnit
  contains
    procedure, private, pass :: &
      Initialize_CD
    generic, public :: &
      Initialize => Initialize_CD
    procedure, public, pass :: &
      RecordObservables
    procedure, public, pass :: &
      SetExtensiveIntensive
    procedure ( CLS_Interface ), public, nopass, deferred :: &
      ClusterLengthScale
  end type ClusterDynamicsTemplate

  abstract interface

    function CLS_Interface ( CP, PotentialEnergy ) result ( CLS )
      use Basics
      use ClusterParticles_Form
      class ( ClusterParticlesForm ) :: CP
      real ( KDR ), intent ( in ) :: PotentialEnergy
      real ( KDR ) :: CLS
    end function CLS_Interface

  end interface

contains


  subroutine Initialize_CD ( CD, EnergyUnit, TimeUnit, PotentialEnergyScale )

    class ( ClusterDynamicsTemplate ), intent ( inout ) :: &
      CD
    type ( MeasuredValueForm ), intent ( in ) :: &
      EnergyUnit, &
      TimeUnit
    procedure ( PES_Interface ), pointer, intent ( in ) :: &
      PotentialEnergyScale

    allocate ( ClusterParticlesForm :: CD % DistributedParticles )
    select type ( CP => CD % DistributedParticles )
    type is ( ClusterParticlesForm )
    call CP % Initialize &
           ( PROGRAM_HEADER % Communicator, PotentialEnergyScale, &
             TimeUnitOption = TimeUnit )
    end select !-- CP

    call CD % Initialize ( TimeUnit = TimeUnit )

    CD % EnergyUnit = EnergyUnit

  end subroutine Initialize_CD


  subroutine RecordObservables ( PD, iCycle )

    class ( ClusterDynamicsTemplate ), intent ( inout ) :: &
      PD
    integer ( KDI ), intent ( in ) :: &
      iCycle

    integer ( KDI ) :: &
      iVrbl  !-- iVariable
    real ( KDR ) :: &
      TimeScale_KE, &
      TimeScale_PE
    real ( KDR ), dimension ( 3 ) :: &
      CenterOfMass
    type ( CollectiveOperation_R_Form ) :: &
      CO

    select type ( CP => PD % DistributedParticles )
    type is ( ClusterParticlesForm )

    associate ( iVl => iCycle + 1 ) !-- iValue offset, iCycle starts at 0

    !-- Extensive

    !-- FIXME: At least with NAG, seems to be necessary to spell out indices of
    !          MOMENTUM and ANGULAR_MOMENTUM, at least in the associate statement
    associate &
      ( KE => PD % Extensive % Value ( iVl, CP % KINETIC_ENERGY ), &
        PE => PD % Extensive % Value ( iVl, CP % POTENTIAL_ENERGY ), &
        TE => PD % Extensive % Value ( iVl, CP % TOTAL_ENERGY ), &
        V  => PD % Extensive % Value ( iVl, CP % VIRIAL ), &
        MM => PD % Extensive % Value &
                 ( iVl, CP % MASS_MOMENT ( 1 ) : CP % MASS_MOMENT ( 3 ) ), &
        P  => PD % Extensive % Value &
                 ( iVl, CP % MOMENTUM ( 1 ) : CP % MOMENTUM ( 3 ) ), &
        L  => PD % Extensive % Value &
                 ( iVl, &
                   CP % ANGULAR_MOMENTUM ( 1 ) : CP % ANGULAR_MOMENTUM ( 3 ) ) )

    !-- CenterOfMass used to compute angular momentum lags by 1 step
    if ( iVl == 1 ) then
      CenterOfMass = 0.0_KDR
    else
      CenterOfMass = PD % Intensive % Value ( iVl - 1, CP % CENTER_OF_MASS )
    end if

    call CO % Initialize &
           ( CP % Communicator, [ CP % N_EXTENSIVE ], [ CP % N_EXTENSIVE ] )
    CO % Outgoing % Value ( CP % KINETIC_ENERGY ) = CP % MyKineticEnergy ( )
    CO % Outgoing % Value ( CP % POTENTIAL_ENERGY ) &
      = CP % MyPotentialEnergy ( PD % PotentialValue )
    CO % Outgoing % Value ( CP % TOTAL_ENERGY ) = 0.0_KDR
    CO % Outgoing % Value ( CP % VIRIAL ) &
      = CP % MyVirial ( PD % VirialValue )
    CO % Outgoing % Value ( CP % MASS_MOMENT ) = CP % MyMassMoment ( )
    CO % Outgoing % Value ( CP % MOMENTUM ) = CP % MyMomentum ( )
    CO % Outgoing % Value ( CP % ANGULAR_MOMENTUM ) &
      = CP % MyAngularMomentum ( CenterOfMass )
    call CO % Reduce ( REDUCTION % SUM )

    KE = CO % Incoming % Value ( CP % KINETIC_ENERGY )
    PE = CO % Incoming % Value ( CP % POTENTIAL_ENERGY )
    TE = KE + PE
     V = CO % Incoming % Value ( CP % VIRIAL )
    MM = CO % Incoming % Value ( CP % MASS_MOMENT )
     P = CO % Incoming % Value ( CP % MOMENTUM )
     L = CO % Incoming % Value ( CP % ANGULAR_MOMENTUM )

    if ( mod ( iCycle, PD % WriteCycleInterval ) == 0 ) then
      call Show ( 'Extensive variables', CONSOLE % INFO_3 )
      do iVrbl = 1, CP % N_EXTENSIVE
        call Show ( PD % Extensive % Value ( iVl, iVrbl ), &
                    PD % Extensive % Unit ( iVrbl ), &
                    PD % Extensive % Variable ( iVrbl ), CONSOLE % INFO_3 )
      end do
    end if

    !-- Intensive

    associate &
      ( CM  => PD % Intensive % Value &
                 ( iVl, &
                   CP % CENTER_OF_MASS ( 1 ) : CP % CENTER_OF_MASS ( 3 ) ), &
        CLS => PD % Intensive % Value ( iVl, CP % CLUSTER_LENGTH_SCALE ), &
        GTS => PD % Intensive % Value ( iVl, CP % GLOBAL_TIME_SCALE ), &
        LTS => PD % Intensive % Value ( iVl, CP % LOCAL_TIME_SCALE ) )

    CM  = MM / ( CP % ParticleMass  *  CP % nParticles )
    CLS = PD % ClusterLengthScale ( CP, PE )
!    LTS = PD % TimeScaleMin

!     !-- TimeStep

!     associate ( M => CP % ParticleMass * CP % nParticles )
!     TimeScale_KE &
!       = PD % TimeScaleExpression &
!           ( EnergyScale = KE, LengthScale = CLS, MassScale = M )
!     TimeScale_PE &
!       = PD % TimeScaleExpression &
!           ( EnergyScale = abs ( PE ), LengthScale = CLS, MassScale = M )
!     call Show ( 'Time scales', CONSOLE % INFO_3 )
!     call Show ( TimeScale_KE, PD % TimeUnit, 'TimeScale_KE', CONSOLE % INFO_3 )
!     call Show ( TimeScale_PE, PD % TimeUnit, 'TimeScale_PE', CONSOLE % INFO_3 )
!     call Show ( LTS, PD % TimeUnit, 'TimeScale_PE_Local', CONSOLE % INFO_3 )
!     PD % TimeStep = PD % TimeStepFactor * min ( TimeScale_KE, TimeScale_PE )
! !    PD % TimeStep &
! !      = PD % TimeStepFactor * min ( TimeScale_KE, TimeScale_PE, LTS )
!     end associate !-- M

!     GTS = TimeScale_PE

    if ( mod ( iCycle, PD % WriteCycleInterval ) == 0 ) then
      call Show ( 'Intensive variables', CONSOLE % INFO_3 )
      do iVrbl = 1, CP % N_INTENSIVE
        call Show ( PD % Intensive % Value ( iVl, iVrbl ), &
                    PD % Intensive % Unit ( iVrbl ), &
                    PD % Intensive % Variable ( iVrbl ), CONSOLE % INFO_3 )
      end do
    end if

    end associate !-- CM, etc.
    end associate !-- KE, etc.
    end associate !-- iVl
    end select !-- CP

  end subroutine RecordObservables


  subroutine SetExtensiveIntensive ( PD )

    class ( ClusterDynamicsTemplate ), intent ( inout ) :: &
      PD

    select type ( CP => PD % DistributedParticles )
    type is ( ClusterParticlesForm )

    CP % N_EXTENSIVE      = CP % N_EXTENSIVE_DP + CP % N_EXTENSIVE_CP
    CP % KINETIC_ENERGY   = 1
    CP % POTENTIAL_ENERGY = 2
    CP % TOTAL_ENERGY     = 3
    CP % VIRIAL           = 4
    CP % MASS_MOMENT      = [  5,  6,  7 ]
    CP % MOMENTUM         = [  8,  9, 10 ]
    CP % ANGULAR_MOMENTUM = [ 11, 12, 13 ]
    call PD % Extensive % Initialize &
           ( [ PD % nCycles + 1, CP % N_EXTENSIVE ], &
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
                                'AngularMomentum_3              ' ], & 
             NameOption = 'Extensive', &
             UnitOption &
               = [ PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   PD % EnergyUnit, &
                   CP % MassUnit * CP % LengthUnit, &
                   CP % MassUnit * CP % LengthUnit, &
                   CP % MassUnit * CP % LengthUnit, &
                   CP % MassUnit * CP % LengthUnit / PD % TimeUnit, &
                   CP % MassUnit * CP % LengthUnit / PD % TimeUnit, &
                   CP % MassUnit * CP % LengthUnit / PD % TimeUnit, &
                   CP % MassUnit * CP % LengthUnit ** 2  /  PD % TimeUnit, &
                   CP % MassUnit * CP % LengthUnit ** 2  /  PD % TimeUnit, &
                   CP % MassUnit * CP % LengthUnit ** 2  /  PD % TimeUnit ] )

    CP % N_INTENSIVE = CP % N_INTENSIVE_DP + CP % N_INTENSIVE_CP
    CP % CENTER_OF_MASS       = [ 1, 2, 3 ]
    CP % CLUSTER_LENGTH_SCALE = 4
    CP % GLOBAL_TIME_SCALE = 5
    CP % LOCAL_TIME_SCALE  = 6
    call PD % Intensive % Initialize &
           ( [ PD % nCycles + 1, CP % N_INTENSIVE ], &
             VariableOption = [ 'CenterOfMass_1                 ', &
                                'CenterOfMass_2                 ', &
                                'CenterOfMass_3                 ', &
                                'ClusterLengthScale             ', &
                                'GlobalTimeScale                ', &
                                'LocalTimeScale                 ' ], &
             NameOption = 'Intensive', &
             UnitOption &
             = [ CP % LengthUnit, CP % LengthUnit, CP % LengthUnit, &
                 CP % LengthUnit, PD % TimeUnit, PD % TimeUnit ] )

    end select !-- CP

  end subroutine SetExtensiveIntensive


end module ClusterDynamics_Template
