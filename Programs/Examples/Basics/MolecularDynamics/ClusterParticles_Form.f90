module ClusterParticles_Form

  use Basics
  use NormalRandomNumber_Function
  use DistributedParticles_Form

  implicit none
  private

  type, public, extends ( DistributedParticlesForm ) :: ClusterParticlesForm
    integer ( KDI ) :: &
      N_EXTENSIVE_CP = 0
    integer ( KDI ) :: &
      CLUSTER_LENGTH_SCALE, &
      GLOBAL_TIME_SCALE, &
      LOCAL_TIME_SCALE
    integer ( KDI ) :: &
      N_INTENSIVE_CP = 3
    real ( KDR ) :: &
      SofteningLength
  contains
    procedure, private, pass :: &
      Initialize_CP
    generic, public :: &
      Initialize => Initialize_CP
  end type ClusterParticlesForm
  
  interface
    function PES_Interface ( CP, ClusterSize ) result ( PES )
      use Basics
      import ClusterParticlesForm
      class ( ClusterParticlesForm ), intent ( inout ) :: CP
      real ( KDR ), intent ( in ) :: ClusterSize
      real ( KDR ) :: PES
    end function PES_Interface
  end interface

  public &
    :: PES_Interface

contains


  subroutine Initialize_CP ( CP, C, PotentialEnergyScale, TimeUnitOption )

    class ( ClusterParticlesForm ), intent ( inout ) :: &
      CP
    type ( CommunicatorForm ), intent ( in ), target :: &
      C
    procedure ( PES_Interface ), pointer, intent ( in ) :: &
      PotentialEnergyScale
    type ( MeasuredValueForm ), intent ( in ), optional :: &
      TimeUnitOption

    integer ( KDI ) :: &
      nParticles
    real ( KDR ) :: &
      ClusterSize, &
      ClusterMass, &
      KineticEnergyRatio, &
      SofteningLengthFactor, &
      BoxLength, &
      ParticleMass
    type ( MeasuredValueForm ) :: &
      LengthUnit

    nParticles = 256
    call PROGRAM_HEADER % GetParameter ( nParticles, 'nParticles' )

    ClusterSize = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( ClusterSize, 'ClusterSize', InputUnitOption = LengthUnit )

    ClusterMass = 1.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( ClusterSize, 'ClusterMass', InputUnitOption = CP % MassUnit )

    KineticEnergyRatio = 0.01_KDR
    call PROGRAM_HEADER % GetParameter &
           ( KineticEnergyRatio, 'KineticEnergyRatio', &
             InputUnitOption = LengthUnit )

    SofteningLengthFactor = 0.01_KDR
    call PROGRAM_HEADER % GetParameter &
           ( SofteningLengthFactor, 'SofteningLengthFactor' )

    CP % SofteningLength = SofteningLengthFactor * ClusterSize

!    BoxLength = 20.0_KDR * ClusterSize
    BoxLength = huge ( 1.0_KDR )

    ParticleMass = ClusterMass / nParticles

    call CP % DistributedParticlesForm % Initialize &
           ( C, BoxLength, ParticleMass, nParticles, &
             LengthUnitOption = LengthUnit, TimeUnitOption = TimeUnitOption )
    CP % IsPeriodic = .false.

    call InitializeRandomSeed ( C )
    call SetParticles &
           ( CP, ClusterSize, KineticEnergyRatio, PotentialEnergyScale )

  end subroutine Initialize_CP


  subroutine SetParticles &
               ( CP, ClusterSize, KineticEnergyRatio, PotentialEnergyScale )

    class ( ClusterParticlesForm ), intent ( inout ) :: &
      CP
    real ( KDR ), intent ( in ) :: &
      ClusterSize, &
      KineticEnergyRatio
    procedure ( PES_Interface ), pointer, intent ( in ) :: &
      PotentialEnergyScale

    integer ( KDI ) :: &
      iP, &  !-- iParticle
      iD     !-- iDimension
    real ( KDR ) :: &
      KineticEnergyPerParticle

    call Show ( 'Setting particle initial conditions', CONSOLE % INFO_3 )

    associate ( MP => CP % MyParticles )

    !-- Initial positions
  
    associate ( PositionVariance => ( 0.5_KDR * ClusterSize ) ** 2 )
    do iP = 1, CP % nMyParticles
      do iD = 1, 3

        MP % Value ( iP, MP % POSITION ( iD ) ) &
          = NormalRandomNumber ( PositionVariance )

      end do
    end do
    end associate !-- PositionVariance

    !-- Positions in box

    MP % Value ( :, MP % POSITION_BOX ) = MP % Value ( :, MP % POSITION )

    !-- Initial velocities

    KineticEnergyPerParticle &
      = KineticEnergyRatio * PotentialEnergyScale ( CP, ClusterSize ) &
        / CP % nParticles

    associate &
      ( VelocityVariance => KineticEnergyPerParticle / CP % ParticleMass )
    do iP = 1, CP % nMyParticles
      do iD = 1, 3

        MP % Value ( iP, MP % VELOCITY ( iD ) ) &
          = NormalRandomNumber ( VelocityVariance )

      end do
    end do
    end associate !-- VelocityVariance

    end associate !-- MP

  end subroutine SetParticles


end module ClusterParticles_Form
