module SetFermiDiracSpectrum_Command

  use Basics

  implicit none
  private

  public :: SetFermiDiracSpectrum

contains


  subroutine SetFermiDiracSpectrum &
               ( Energy, Temperature, ChemicalPotential, EnergyDensity )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Energy
    real ( KDR ), intent ( in ) :: &
      Temperature, &
      ChemicalPotential
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      EnergyDensity

    real ( KDR ) :: &
      LogHuge
    real ( KDR ), dimension ( size ( Energy ) ) :: &
      ExpArgument

    associate &
      ( J  => EnergyDensity, &
        E  => Energy, &
        T  => Temperature, &
        Mu => ChemicalPotential, &
        kB => CONSTANT % BOLTZMANN )

    LogHuge = log ( huge ( 1.0_KDR ) )
    ExpArgument = ( E - Mu ) / ( kB * T )
    where ( ExpArgument > LogHuge )
      ExpArgument = LogHuge
    end where

    J = E / ( exp ( ExpArgument ) + 1.0_KDR )

    end associate !-- J, etc.

  end subroutine SetFermiDiracSpectrum



end module SetFermiDiracSpectrum_Command
