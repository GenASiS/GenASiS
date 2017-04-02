module SetPlanckSpectrum_Command

  use Basics

  implicit none
  private

  public :: &
    SetPlanckSpectrum

contains


  subroutine SetPlanckSpectrum ( Energy, Temperature, EnergyDensity )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      Temperature
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EnergyDensity

    real ( KDR ) :: &
      hc_Minus_3, &
      LogHuge
    real ( KDR ), dimension ( size ( Energy ) ) :: &
      ExpArgument

    associate &
      ( J      => EnergyDensity, &
        E      => Energy, &
        T      => Temperature, &
        kB     => CONSTANT % BOLTZMANN, &
        g      => 2.0_KDR, &  !-- Degeneracy
        hBar_c => CONSTANT % PLANCK_REDUCED  *  CONSTANT % SPEED_OF_LIGHT, &
        TwoPi  => 2.0_KDR * CONSTANT % PI )

    hc_Minus_3 =  ( TwoPi * hBar_c ) ** (-3)

    LogHuge = log ( huge ( 1.0_KDR ) )
    ExpArgument = E / ( kB * T )
    where ( ExpArgument > LogHuge )
      ExpArgument = LogHuge
    end where

    J  =  g * hc_Minus_3  *  E / ( exp ( ExpArgument ) - 1.0_KDR )

    end associate !-- J, etc.

  end subroutine SetPlanckSpectrum


end module SetPlanckSpectrum_Command
