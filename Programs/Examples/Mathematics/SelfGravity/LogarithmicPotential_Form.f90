!! Parameters for LogarithmicPotential_1 taken from 
!!   Muller Stienmetz 1995 for comparison.
!! Parameters for LogarithmicPotential_2 & LogarithmicPotential_3 
!!   taken from Binney and Tremaine "Galactic Dynamics"

!! Analytic solutions taken from Binney and Tremaine "Galactic Dynamics" second edition 
!!   Chapter 2 Section 2.3.2 pg 75-76, Eqns. 2.71a & 2.71c.

! Constructs DD % N_equations density profiles that corespond 
!   to a logarithmic potential.

module LogarithmicPotential_Form

  use Basics
  use Mathematics
  use DensityDistribution_Template
  
  implicit none

  type, public, extends ( DensityDistributionTemplate ) :: &
    LogarithmicForm
  contains
    procedure, public, pass :: &
      SetLogarithmic
  end type LogarithmicForm

    private :: &
      SetLogarithmicKernal

contains


  subroutine SetLogarithmic ( L, v0, rho_c, q_phi )
    class ( LogarithmicForm ), intent ( inout ) :: &
      L
    real ( KDR ), dimension ( L % N_Equations ) :: &
      v0, &
      rho_c, &
      q_phi

    integer ( KDI ) :: &
      iE !-- iEquation

    associate ( A  => L % Atlas )
    associate ( SA => L % Source )
    associate ( RA => L % Reference )

    do iE = 1, L % N_Equations
      call SetLogarithmicKernal &
             ( SA, RA, A, v0 ( iE ), rho_c ( iE ), q_phi ( iE ), iE )
    end do

    end associate !-- RA
    end associate !-- SA
    end associate !-- A

  end subroutine SetLogarithmic


  subroutine SetLogarithmicKernal &
               ( Source_ASC, Reference_ASC, A, v0, rho_c, q_phi, iVariable )
    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      v0, &
      rho_c, &
      q_phi
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( StorageForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      rho_sq, & !-- X^2 + Y^2 
      Z_sq
    integer ( KDI ) :: &
      i

    !-- Geometry

    G => A % Geometry ( )

    associate &
      (  R    => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ) )

    allocate &
      ( rho_sq ( size ( R ) ) , &
        Z_sq   ( size ( R ) ) )
         
    rho_sq = ( R * sin ( Theta ) ) ** 2
    Z_sq   = ( R * cos ( Theta ) ) ** 2

    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D     => Source % Value ( :, iVariable ) )

    D = ( v0 / q_phi ) ** 2 & 
          * ( (  2 * q_phi ** 2 + 1.0_KDR ) * rho_c ** 2 &
                + rho_sq + ( 2.0_KDR - q_phi ** ( - 2 ) ) * Z_sq ) &
          / ( rho_c ** 2 + rho_sq + Z_sq / q_phi ** 2 ) ** 2

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI )

    Phi = v0 ** 2 * log ( rho_c ** 2 + rho_sq + Z_sq / ( q_phi ** 2 ) ) / 2

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq )
    nullify ( G, Source, Reference )

  end subroutine SetLogarithmicKernal

end module LogarithmicPotential_Form
