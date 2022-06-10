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
      TwoPi, &
      hBar_c, &
      hc_M_3, &
      kB, &
      LogSqrtHuge
    real ( KDR ), dimension ( size ( Energy ) ) :: &
      ExpArgument

    associate &
      ( J  =>  EnergyDensity, &
        E  =>  Energy, &
        T  =>  Temperature, &
        g  =>  2.0_KDR )  !-- Degeneracy

    !-- FIXME: IBM XL compiler does not like associate with CONSTANT members
    TwoPi   =  2.0_KDR * CONSTANT % PI 
    hBar_c  =  CONSTANT % PLANCK_REDUCED  *  CONSTANT % SPEED_OF_LIGHT
    hc_M_3  =  ( TwoPi * hBar_c ) ** (-3)
    kB      =  CONSTANT % BOLTZMANN

    LogSqrtHuge = log ( sqrt ( huge ( 1.0_KDR ) ) )
    ExpArgument = E / ( kB * T )
    where ( ExpArgument > LogSqrtHuge )
      ExpArgument = LogSqrtHuge
    end where

    J  =  g * hc_M_3  *  E / ( exp ( ExpArgument ) - 1.0_KDR )

    end associate !-- J, etc.

  end subroutine SetPlanckSpectrum


end module SetPlanckSpectrum_Command
