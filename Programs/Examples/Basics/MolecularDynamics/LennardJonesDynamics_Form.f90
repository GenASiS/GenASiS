module LennardJonesDynamics_Form

  use Basics
  use LatticeDynamics_Template

  implicit none
  private

  type, public, extends ( LatticeDynamicsTemplate ) :: LennardJonesDynamicsForm
    real ( KDR ) :: &
      EnergyParameter, &
      LengthParameter
  contains
    procedure, private, pass :: &
      Initialize_LJD
    generic, public :: &
      Initialize => Initialize_LJD
    procedure, public, pass :: &
      Potential
    procedure, public, pass :: &
      Virial
    procedure, public, pass :: &
      Force  
  end type LennardJonesDynamicsForm

contains


  subroutine Initialize_LJD ( LJD )

    class ( LennardJonesDynamicsForm ), intent ( inout ) :: &
      LJD

    type ( QuantityForm ) :: &
      EnergyUnit, &
      LengthUnit, &
      TimeUnit

    LJD % Type = 'a Lennard-Jones lattice' 

    LJD % EnergyParameter = 0.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( LJD % EnergyParameter, 'EnergyParameter', &
             InputUnitOption = EnergyUnit )

    LJD % LengthParameter = 0.0_KDR
    call PROGRAM_HEADER % GetParameter &
           ( LJD % LengthParameter, 'LengthParameter', &
             InputUnitOption = LengthUnit )

    TimeUnit = UNIT % IDENTITY
    if ( EnergyUnit == UNIT % ELECTRON_VOLT &
         .or. LengthUnit == UNIT % ANGSTROM ) TimeUnit = UNIT % FEMTOSECOND

    call LJD % Initialize &
           ( EnergyUnit = EnergyUnit, TimeUnit = TimeUnit )

    LJD % TimeStep &
      = LJD % TimeStepFactor &
          * LJD % TimeScaleExpression &
              ( EnergyScale = LJD % EnergyParameter, &
                LengthScale = LJD % LengthParameter, &
                MassScale = LJD % DistributedParticles % ParticleMass )


  end subroutine Initialize_LJD


  function Potential ( PD, Distance ) result ( P )

    class ( LennardJonesDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ) :: &
      P

    real ( KDR ), dimension ( size ( Distance ) ) :: &
      S_R

    associate &
      ( Epsilon => PD % EnergyParameter, &
        Sigma   => PD % LengthParameter )

    where ( Distance > 0.0_KDR )
      S_R = Sigma / Distance
    elsewhere
      S_R = 0.0_KDR
    end where

    P = 4.0_KDR * Epsilon * sum ( S_R ** 12 - S_R ** 6 )

    end associate !-- Epsilon, Sigma

  end function Potential


  function Virial ( PD, Distance ) result ( V )

    class ( LennardJonesDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ) :: &
      V

    real ( KDR ), dimension ( size ( Distance ) ) :: &
      S_R

    associate &
      ( Epsilon => PD % EnergyParameter, &
        Sigma   => PD % LengthParameter )

    where ( Distance > 0.0_KDR )
      S_R = Sigma / Distance
    elsewhere
      S_R = 0.0_KDR
    end where

    V = 24.0_KDR * Epsilon * sum ( 2.0_KDR * S_R ** 12 - S_R ** 6 )

    end associate !-- Epsilon, Sigma

  end function Virial


  function Force ( PD, Displacement, Distance ) result ( F )

    class ( LennardJonesDynamicsForm ), intent ( inout ) :: &
      PD
    real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      Displacement
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Distance
    real ( KDR ), dimension ( 3 ) :: &
      F

    integer ( KDR ) :: &
      iD
    real ( KDR ), dimension ( size ( Distance ) ) :: &
      S_R
    real ( KDR ), dimension ( size ( Distance ), 3 ) :: &
      R_R2

    associate &
      ( Epsilon => PD % EnergyParameter, &
        Sigma   => PD % LengthParameter )

    where ( Distance > 0.0_KDR )
      S_R = Sigma / Distance
    elsewhere
      S_R = 0.0_KDR
    end where

    do iD = 1, 3

      where ( Distance > 0.0_KDR )
        R_R2 ( :, iD ) = Displacement ( :, iD ) / Distance ** 2
      elsewhere
        R_R2 ( :, iD ) = 0.0_KDR
      end where

      F ( iD ) = 24.0_KDR * Epsilon &
                 * sum ( R_R2 ( :, iD ) * ( 2.0_KDR * S_R ** 12 - S_R ** 6 ) )

    end do

    end associate !-- Epsilon, Sigma

  end function Force


end module LennardJonesDynamics_Form
