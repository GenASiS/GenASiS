module GravitationalDynamics_Form

  use Basics
  use ClusterParticles_Form
  use ClusterDynamics_Template

  implicit none
  private

  type, public, extends ( ClusterDynamicsTemplate ) :: GravitationalDynamicsForm
  contains
    procedure, private, pass :: &
      Initialize_GD
    generic, public :: &
      Initialize => Initialize_GD
    procedure, public, pass :: &
      Potential
    procedure, public, pass :: &
      Virial
    procedure, public, pass :: &
      Force  
    procedure, public, nopass :: &
      ClusterLengthScale
  end type GravitationalDynamicsForm
    
    private :: &
      PotentialEnergyScale

contains


  subroutine Initialize_GD ( GD )

    class ( GravitationalDynamicsForm ), intent ( inout ) :: &
      GD

    type ( MeasuredValueForm ) :: &
      EnergyUnit, &
      TimeUnit
    procedure ( PES_Interface ), pointer :: &
      PES

    GD % Type = 'a gravitational cluster' 

    EnergyUnit = UNIT % IDENTITY
    TimeUnit   = UNIT % IDENTITY

    PES => PotentialEnergyScale

    call GD % Initialize &
           ( EnergyUnit = EnergyUnit, TimeUnit = TimeUnit, &
             PotentialEnergyScale = PES )

    select type ( CP => GD % DistributedParticles )
    class is ( ClusterParticlesForm )
    associate &
      ( G => 1.0_KDR, &
        M => CP % ParticleMass  *  CP % nParticles, &
        L => CP % SofteningLength )
    GD % TimeStep &
      = GD % TimeScaleExpression &
          ( EnergyScale = G * M * M / L, &
            LengthScale = L, &
            MassScale = M )
    end associate !-- G, etc.
    end select !-- CP

  end subroutine Initialize_GD


  function Potential ( PD, Distance ) result ( P )

    class ( GravitationalDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ) :: &
      P

    real ( KDR ), dimension ( size ( Distance ) ) :: &
      One_R

    associate ( DP => PD % DistributedParticles )
    associate &
      ( G => 1.0_KDR, &
        M => DP % ParticleMass ) 

    where ( Distance > 0.0_KDR )
      One_R = 1.0_KDR / Distance
    elsewhere
      One_R = 0.0_KDR
    end where

    P = - G * M * M * sum ( One_R )

    end associate !-- G, M
    end associate !-- DP

  end function Potential


  function Virial ( PD, Distance ) result ( V )

    class ( GravitationalDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ) :: &
      V

    real ( KDR ), dimension ( size ( Distance ) ) :: &
      One_R

    associate ( DP => PD % DistributedParticles )
    associate &
      ( G => 1.0_KDR, &
        M => DP % ParticleMass ) 

    where ( Distance > 0.0_KDR )
      One_R = 1.0_KDR / Distance
    elsewhere
      One_R = 0.0_KDR
    end where

    V = - G * M * M * sum ( One_R )

    end associate !-- G, M
    end associate !-- DP

  end function Virial


  function Force ( PD, Displacement, Distance ) result ( F )

    class ( GravitationalDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Displacement
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ), dimension ( 3 ) :: &
      F

    integer ( KDR ) :: &
      iD
    real ( KDR ), dimension ( size ( Distance ), 3 ) :: &
      R_R3

    select type ( CP => PD % DistributedParticles )
    class is ( ClusterParticlesForm )

    associate &
      ( G => 1.0_KDR, &
        M => CP % ParticleMass, &
        Epsilon => CP % SofteningLength ) 

    do iD = 1, 3

      where ( Distance > 0.0_KDR )
        R_R3 ( :, iD ) &
          = Displacement ( :, iD ) &
            / ( Epsilon ** 2  +  Distance ** 2 ) ** 1.5_KDR
      elsewhere
        R_R3 ( :, iD ) = 0.0_KDR
      end where

      F ( iD ) = - G * M * M * sum ( R_R3 ( :, iD ) )

    end do

    end associate !-- G, M
    end select !-- CP

  end function Force


  function ClusterLengthScale ( CP, PotentialEnergy ) result ( CLS )

    class ( ClusterParticlesForm ) :: &
      CP
    real ( KDR ), intent ( in ) :: &
      PotentialEnergy
    real ( KDR ) :: &
      CLS

    associate &
      ( G => 1.0_KDR, &
        M_Total => CP % ParticleMass * CP % nParticles ) 

    CLS = - G * M_Total * M_Total / PotentialEnergy

    end associate !-- G, M_Total

  end function ClusterLengthScale


  function PotentialEnergyScale ( CP, ClusterSize ) result ( PES )

    class ( ClusterParticlesForm ), intent ( inout ) :: &
      CP
    real ( KDR ), intent ( in ) :: &
      ClusterSize
    real ( KDR ) :: &
      PES

    associate &
      ( G => 1.0_KDR, &
        M_Total => CP % ParticleMass * CP % nParticles ) 

    PES = G * M_Total * M_Total / ClusterSize

    end associate !-- G, M_Total

  end function PotentialEnergyScale


end module GravitationalDynamics_Form
