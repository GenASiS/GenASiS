module SetFermiDiracSpectrum_Command

  use Basics

  implicit none
  private

  public :: &
    SetFermiDiracSpectrum

contains


  subroutine SetFermiDiracSpectrum &
               ( Energy, Temperature, ChemicalPotential, EnergyDensity, &
                 NumberDensity )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      Temperature, &
      ChemicalPotential
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EnergyDensity, &
      NumberDensity

    real ( KDR ) :: &
      TwoPi, &
      hBar_c, &
      hc_M_3, &
      kB, &
      LogSqrtHuge
    real ( KDR ), dimension ( size ( Energy ) ) :: &
      ExpArgument

    associate &
      ( N   =>  NumberDensity, &
        J   =>  EnergyDensity, &
        E   =>  Energy, &
        T   =>  Temperature, &
        Mu  =>  ChemicalPotential, &
        g   =>  1.0_KDR )  !-- Degeneracy

    !-- FIXME: IBM XL compiler does not like associate with CONSTANT members
    TwoPi   =  2.0_KDR * CONSTANT % PI 
    hBar_c  =  CONSTANT % PLANCK_REDUCED  *  CONSTANT % SPEED_OF_LIGHT
    hc_M_3  =  ( TwoPi * hBar_c ) ** (-3)
    kB      =  CONSTANT % BOLTZMANN

    LogSqrtHuge = log ( sqrt ( huge ( 1.0_KDR ) ) )
    ExpArgument = ( E - Mu ) / ( kB * T )
    where ( ExpArgument > LogSqrtHuge )
      ExpArgument = LogSqrtHuge
    end where

    N  =  g * hc_M_3  /  ( exp ( ExpArgument ) + 1.0_KDR )

    J  =  E  *  N

    end associate !-- J, etc.

  end subroutine SetFermiDiracSpectrum


end module SetFermiDiracSpectrum_Command
