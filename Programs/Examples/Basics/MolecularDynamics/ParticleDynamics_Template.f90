module ParticleDynamics_Template

  use Basics
  use DistributedParticles_Form

  implicit none
  private

  type, public, abstract :: ParticleDynamicsTemplate
    integer ( KDI ) :: &
      iCycle, &
      nCycles, &
      nWrite, &
      WriteCycleInterval
    real ( KDR ) :: &
      Time, &
      TimeStep, &
      TimeStepFactor, &
      TimeScaleMin
    real ( KDR ), dimension ( : ), allocatable :: &
      TimeValue
    type ( MeasuredValueForm ) :: &
      TimeUnit
    real ( KDR ), dimension ( : ), allocatable :: &
      PotentialValue, &
      VirialValue
    real ( KDR ), dimension ( :, : ), allocatable :: &
      ForceValue
    character ( LDL ) :: &
      Type = ''
    type ( VariableGroupForm ) :: &
      Extensive, &
      Intensive
    type ( GridImageStreamForm ) :: &
      GridImageStream
    type ( CurveImageForm ) :: &
      TimeSeries
    class ( DistributedParticlesForm ), allocatable :: &
      DistributedParticles
    procedure ( StepInterface ), pointer :: &
      Step => null ( )
  contains
    procedure, private, pass :: &
      Initialize_PD
    generic, public :: &
      Initialize => Initialize_PD
    procedure, public, pass :: &
      Evolve
    procedure, public, nopass :: &
      TimeScaleExpression
    procedure ( PotentialInterface ), public, pass, deferred :: &
      Potential
    procedure ( PotentialInterface ), public, pass, deferred :: &
      Virial
    procedure ( ForceInterface ), public, pass, deferred :: &
      Force
    procedure ( RecordObservablesInterface ), public, pass, deferred :: &
      RecordObservables
    procedure ( SetExtensiveIntensiveInterface ), public, pass, deferred :: &
      SetExtensiveIntensive
  end type ParticleDynamicsTemplate

  abstract interface

    subroutine StepInterface ( PD )
      import ParticleDynamicsTemplate
      class ( ParticleDynamicsTemplate ), intent ( inout ) :: PD
    end subroutine StepInterface

    function PotentialInterface ( PD, Distance ) result ( P )
      use Basics
      import ParticleDynamicsTemplate
      class ( ParticleDynamicsTemplate ), intent ( inout ) :: PD
      real ( KDR ), dimension ( : ), intent ( in ) :: Distance
      real ( KDR ) :: P
    end function PotentialInterface

    function ForceInterface ( PD, Displacement, Distance ) result ( F )
      use Basics
      import ParticleDynamicsTemplate
      class ( ParticleDynamicsTemplate ), intent ( inout ) :: PD
      real ( KDR ), dimension ( :, : ), intent ( in ) :: Displacement
      real ( KDR ), dimension ( : ), intent ( in ) :: Distance
      real ( KDR ), dimension ( 3 ) :: F
    end function ForceInterface

    subroutine RecordObservablesInterface ( PD, iCycle )
      use Basics
      import ParticleDynamicsTemplate
      class ( ParticleDynamicsTemplate ), intent ( inout ) :: PD
      integer ( KDI ), intent ( in ) :: iCycle
    end subroutine RecordObservablesInterface

    subroutine SetExtensiveIntensiveInterface ( PD )
      import ParticleDynamicsTemplate
      class ( ParticleDynamicsTemplate ), intent ( inout ) :: PD
    end subroutine SetExtensiveIntensiveInterface

  end interface

    private :: &
      WriteTimeSeries, &
      StepLeapfrog, &
      StepVelocityVerlet

      private :: &
        ComputeForce, &
        StepVelocity, &
        StepPosition

        private :: &
          BinParticleDistances

contains


  subroutine Initialize_PD ( PD, TimeUnit )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD
    type ( MeasuredValueForm ), intent ( in ) :: &
      TimeUnit

    character ( LDL ) :: &
      Step

    call Show ( 'Initializing ' // trim ( PD % Type ), CONSOLE % INFO_3 )

    PD % iCycle  = 0
    PD % nCycles = 1000
    PD % nWrite  = 100
    call PROGRAM_HEADER % GetParameter ( PD % nCycles, 'nCycles' )
    call PROGRAM_HEADER % GetParameter ( PD % nWrite, 'nWrite' )

    PD % WriteCycleInterval = PD % nCycles / PD % nWrite
    call Show ( PD % WriteCycleInterval, 'WriteCycleInterval' )

    PD % TimeUnit = TimeUnit

    PD % TimeStepFactor = 0.004_KDR
    call PROGRAM_HEADER % GetParameter &
           ( PD % TimeStepFactor, 'TimeStepFactor' )

    Step = 'LEAPFROG'
    call PROGRAM_HEADER % GetParameter ( Step, 'Step' )
    select case ( trim ( Step ) )
    case ( 'LEAPFROG' )
      PD % Step => StepLeapfrog
    case ( 'VELOCITY_VERLET' )
      PD % Step => StepVelocityVerlet
    case default
      call Show ( 'Step not implemented', CONSOLE % ERROR )
      call Show ( Step, 'Step', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select

    !-- DistributedParticles must be initialized by the extension

    associate &
      ( DP => PD % DistributedParticles, &
        MP => PD % DistributedParticles % MyParticles )
    allocate ( PD % PotentialValue ( DP % nMyParticles ) )
    allocate ( PD % VirialValue ( DP % nMyParticles ) )
    allocate ( PD % ForceValue ( DP % nMyParticles, 3 ) )
    allocate ( PD % TimeValue ( 0 : PD % nCycles ) )
    end associate !-- DP, MP

  end subroutine Initialize_PD


  subroutine Evolve ( PD )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD

    integer ( KDI ) :: &
      iC, &  !-- iCycle
      iTimerComputation
    type ( VariableGroupForm ), dimension ( 1 ) :: &
      VGP

    call PROGRAM_HEADER % AddTimer ( 'Computational', iTimerComputation )

    call Show ( 'Evolving particles', CONSOLE % INFO_3 )

    associate &
      ( DP => PD % DistributedParticles, &
        MP => PD % DistributedParticles % MyParticles, &
        T => PROGRAM_HEADER % Timer ( iTimerComputation ) )

    PD % Time = 0.0_KDR
    call Show ( PD % Time, PD % TimeUnit, 'Time', CONSOLE % INFO_3 )

    call PD % SetExtensiveIntensive ( )
    call VGP ( 1 ) % Initialize ( MP )
    call DP % SetImage ( VGP, VGP, PROGRAM_HEADER % Name )

    PD % TimeValue ( 0 ) = PD % Time
    call ComputeForce ( PD )  !-- needed to get initial potential energy
    call PD % RecordObservables ( iCycle = 0 )
    call DP % Write &
           ( TimeOption = PD % Time / PD % TimeUnit, &
             CycleNumberOption = PD % iCycle )

    call T % Start ( ) 

    do iC = 1, PD % nCycles

      call Show ( 'Stepping particles', CONSOLE % INFO_3 )
      PD % iCycle = PD % iCycle + 1
      call Show ( PD % iCycle, 'iCycle', CONSOLE % INFO_3 )
      call Show ( PD % TimeStep, PD % TimeUnit, 'TimeStep', CONSOLE % INFO_3 )

      call PD % Step ( )

      PD % Time = PD % Time + PD % TimeStep
      call Show ( PD % Time, PD % TimeUnit, 'Time', CONSOLE % INFO_3 )

      PD % TimeValue ( iC ) = PD % Time
      call PD % RecordObservables ( iC )

      if ( mod ( iC, PD % WriteCycleInterval ) == 0 ) then

        call T % Stop ( )

        call DP % Write &
               ( TimeOption = PD % Time / PD % TimeUnit, &
                 CycleNumberOption = PD % iCycle )

        call T % Start ( )
        
      end if

    end do

    call T % Stop ( )
    
    call WriteTimeSeries ( PD )

    end associate !-- DP, etc.

  end subroutine Evolve


  function TimeScaleExpression &
             ( EnergyScale, LengthScale, MassScale ) result ( TS )

    real ( KDR ), intent ( in ) :: &
      EnergyScale, &
      LengthScale, &
      MassScale
    real ( KDR ) :: &
      TS

    if ( EnergyScale > 0.0_KDR ) then
      TS = sqrt ( MassScale / EnergyScale ) * LengthScale
    else
      TS = sqrt ( huge ( 1.0_KDR ) )
    end if
  
  end function TimeScaleExpression


  subroutine WriteTimeSeries ( PD )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD

    character ( LDF ) :: &
      OutputDirectory

    if ( PROGRAM_HEADER % Communicator % Rank /= CONSOLE % DisplayRank ) return

    OutputDirectory = '../Output/'
    call PROGRAM_HEADER % GetParameter ( OutputDirectory, 'OutputDirectory' )    

    associate ( GIS => PD % GridImageStream )
    call GIS % Initialize &
           ( trim ( PROGRAM_HEADER % Name ) // '_TimeSeries', &
             WorkingDirectoryOption = OutputDirectory )
    call GIS % Open ( GIS % ACCESS_CREATE, SeriesOption = .false. )

    associate ( TS => PD % TimeSeries )
    call TS % Initialize ( GIS ) 
    call TS % SetGrid  &
           ( Directory = 'TimeSeries', &
             NodeCoordinate = PD % TimeValue, &
             nProperCells = size ( PD % TimeValue ), oValue = 0, &
             CoordinateUnitOption = PD % TimeUnit, &
             CoordinateLabelOption = 't' )
    call TS % AddVariableGroup ( PD % Extensive )
    call TS % AddVariableGroup ( PD % Intensive )
    call TS % Write ( )
    end associate !-- TS

    call GIS % Close ( )
    end associate !-- GIS

  end subroutine WriteTimeSeries


  subroutine StepLeapfrog ( PD )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD

    associate ( DP => PD % DistributedParticles )
    associate ( MP => DP % MyParticles )
    !-- At least with NAG, R => MP % Value ( :, MP % POSITION ) doesn't leave
    !   values altered after return
    associate &
      ( R => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) ), &
        R_Box => MP % Value &
                   ( :, MP % POSITION_BOX ( 1 ) : MP % POSITION_BOX ( 3 ) ), &
        V => MP % Value ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )

    !-- Assume force at t already computed

    !-- Advance velocity from t - dt/2 to t + dt/2
    call StepVelocity &
           ( V, PD % ForceValue, DP % ParticleMass, PD % TimeStep, &
             UpdateCoefficient = 1.0_KDR )

    !-- Advance position from t to t + dt
    call StepPosition &
           ( R, R_Box, DP % IsPeriodic, V, DP % BoxLength, PD % TimeStep )

    !-- Compute force at t + dt
    call ComputeForce ( PD )

    end associate !-- R, V
    end associate !-- MP
    end associate !-- DP

  end subroutine StepLeapfrog


  subroutine StepVelocityVerlet ( PD )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD

    associate ( DP => PD % DistributedParticles )
    associate ( MP => DP % MyParticles )
    !-- At least with NAG, R => MP % Value ( :, MP % POSITION ) doesn't leave
    !   values altered after return
    associate &
      ( R => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) ), &
        R_Box => MP % Value &
                   ( :, MP % POSITION_BOX ( 1 ) : MP % POSITION_BOX ( 3 ) ), &
        V => MP % Value ( :, MP % VELOCITY ( 1 ) : MP % VELOCITY ( 3 ) ) )

    !-- Assume force at t already computed

    !-- First half of velocity update to t + dt
    call StepVelocity &
           ( V, PD % ForceValue, DP % ParticleMass, PD % TimeStep, &
             UpdateCoefficient = 0.5_KDR )

    !-- Advance position from t to t + dt, using partial velocity update
    call StepPosition &
           ( R, R_Box, DP % IsPeriodic, V, DP % BoxLength, PD % TimeStep )

    !-- Compute force at t + dt
    call ComputeForce ( PD )

    !-- Second half of velocity update to t + dt
    call StepVelocity &
           ( V, PD % ForceValue, DP % ParticleMass, PD % TimeStep, &
             UpdateCoefficient = 0.5_KDR )

    end associate !-- R, V
    end associate !-- MP
    end associate !-- DP

  end subroutine StepVelocityVerlet


  subroutine ComputeForce ( PD )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD

    integer ( KDI ) :: &
      iR, &  !-- iRank
      iP, &  !-- iParticle
      iD     !-- iDimension
    real ( KDR ), dimension ( : ), allocatable :: &
      Distance
    real ( KDR ), dimension ( :, : ), allocatable :: &
      Displacement
    real ( KDR ), dimension ( :, : ), pointer :: &
      R_Local

    associate ( DP => PD % DistributedParticles )
    associate &
      ( C   => DP % Communicator, &
        MP  => DP % MyParticles, &
        nMP => DP % nMyParticles, &
        L   => DP % BoxLength )
    associate ( R_Guest => DP % GuestPosition )

    if ( DP % IsPeriodic ) then
      R_Local &
        => MP % Value ( :, MP % POSITION_BOX ( 1 ) : MP % POSITION_BOX ( 3 ) )
    else
      R_Local => MP % Value ( :, MP % POSITION ( 1 ) : MP % POSITION ( 3 ) )
    end if

    allocate ( Distance ( nMP ) )
    allocate ( Displacement ( nMP, 3 ) )

    call Clear ( PD % PotentialValue )
    call Clear ( PD % ForceValue )
    call Clear ( PD % VirialValue )
    call Clear ( DP % MyPairCount )
    PD % TimeScaleMin = huge ( 1.0_KDR )

    do iR = 0, C % Size - 1

      if ( iR == 0 ) then
        call Copy ( R_Local, R_Guest )
      else
        call DP % Incoming % Wait ( )
        call Copy ( reshape ( DP % Incoming % Value, [ nMP, 3 ] ), R_Guest )
        call DP % Outgoing % Wait ( )
      end if
      
      if ( iR < C % Size - 1 ) then
        call DP % Incoming % Receive ( )
        call Copy ( reshape ( R_guest, [ 3 * nMP ] ), DP % Outgoing % Value )
        call DP % Outgoing % Send ( )
      end if

      do iP = 1, nMP
        call Clear ( Distance )
        do iD = 1, 3
          Displacement ( :, iD ) = R_Local ( iP, iD ) - R_Guest ( :, iD )
          if ( DP % IsPeriodic ) then
            !-- Minimum image convention
            Displacement ( :, iD ) &
              = Displacement ( :, iD ) - nint ( Displacement ( :, iD ) / L ) * L 
          end if
          Distance = Distance + Displacement ( :, iD ) ** 2
        end do
        Distance = sqrt ( Distance )
        PD % PotentialValue ( iP ) &
          = PD % PotentialValue ( iP ) + PD % Potential ( Distance )
        PD % VirialValue ( iP ) &
          = PD % VirialValue ( iP ) + PD % Virial ( Distance )
        PD % ForceValue ( iP, : ) &
          = PD % ForceValue ( iP, : ) + PD % Force ( Displacement, Distance )
        call BinParticleDistances ( PD, Distance )
        associate &
          ( G => CONSTANT % GRAVITATIONAL, &
            M => DP % ParticleMass, &
            R => minval ( Distance, mask = Distance > 0.0_KDR ) )
        PD % TimeScaleMin &
          = min ( PD % TimeScaleMin, &
                  PD % TimeScaleExpression &
                    ( EnergyScale = abs ( PD % PotentialValue ( iP ) ), &
                      LengthScale = R, MassScale = M ) )
        end associate !-- M, R
      end do !-- iP

    end do !-- iR

    nullify ( R_Local )

    end associate !-- R_Guest
    end associate !-- C, etc.
    end associate !-- DP

  end subroutine ComputeForce


  subroutine StepVelocity ( V, F, M, dt, UpdateCoefficient )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      V
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      F
    real ( KDR ), intent ( in ) :: &
      M, &
      dt, &
      UpdateCoefficient

    V = V + UpdateCoefficient * dt * F / M

  end subroutine StepVelocity


  subroutine StepPosition ( R, R_Box, IsPeriodic, V, L, dt )

    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      R, &
      R_Box
    logical ( KDL ), intent ( in ) :: &
      IsPeriodic
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      V
    real ( KDR ), intent ( in ) :: &
      L, &
      dt

    R = R + dt * V

    if ( IsPeriodic ) then
      !-- Position inside the periodic box
      associate ( Half_L => 0.5_KDR * L )
      R_Box = sign ( 1.0_KDR, R ) * ( mod ( abs ( R ) + Half_L, L ) - Half_L )
      end associate !-- Half_L
    else
      R_Box = R
    end if

  end subroutine StepPosition


  subroutine BinParticleDistances ( PD, Distance )

    class ( ParticleDynamicsTemplate ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance

    integer ( KDI ) :: &
      iP, &  !-- iParticle
      iB     !-- iBin

    associate ( DP => PD % DistributedParticles )
    associate ( L  => DP % BoxLength )

    do iP = 1, size ( Distance )
      if ( Distance ( iP ) > 0.0_KDR .and. Distance ( iP ) < 0.5_KDL * L ) then
        call Search ( DP % CorrelationBinEdge, Distance ( iP ), iB )
        DP % MyPairCount ( iB ) = DP % MyPairCount ( iB ) + 1
      end if
    end do

    end associate !-- L
    end associate !-- DP

  end subroutine BinParticleDistances


end module ParticleDynamics_Template
