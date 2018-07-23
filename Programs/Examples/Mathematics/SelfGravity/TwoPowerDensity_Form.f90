!! Analytic solutions taken from Binney and Tremaine "Galactic Dynamics" second edition 
!!   Chapter 2 Section 2.2.2 pgs 70-71.

!! Analytic solutions only included for :
!! Jaffe     model                     - alpha = 2, beta = 4
!! Hernquist model                     - alpha = 1, beta = 4
!! Navarro, Frenk, & White (NFW) model - alpha = 1, beta = 3

! Constructs DD % N_equations two-power density profiles

module TwoPowerDensity_Form

  use Basics
  use Mathematics
  use DensityDistribution_Template
  
  implicit none

  type, public, extends ( DensityDistributionTemplate ) :: &
    TwoPowerDensityForm
  contains
    procedure, public, pass :: &
      SetTwoPowerDensity
  end type TwoPowerDensityForm

    private :: &
      SetTwoPowerDensityKernal

contains


  subroutine SetTwoPowerDensity ( TP, Density, rr, alpha, beta )
    class ( TwoPowerDensityForm ), intent ( inout ) :: &
      TP
    real ( KDR ), intent ( in ) :: &
      Density
    real ( KDR ), dimension ( TP % N_Equations ) :: &
      rr, &     !-- rr = reference radius = a
      alpha, &
      beta

    integer ( KDI ) :: &
      iE !-- iEquation

    associate ( A  => TP % Atlas )
    associate ( SA => TP % Source )
    associate ( RA => TP % Reference )

    do iE = 1, TP % N_Equations
      call SetTwoPowerDensityKernal &
             ( SA, RA, A, Density, rr ( iE ), alpha ( iE), beta ( iE ), iE )
    end do

    end associate !-- RA
    end associate !-- SA
    end associate !-- A

  end subroutine SetTwoPowerDensity


  subroutine SetTwoPowerDensityKernal &
               ( Source_ASC, Reference_ASC, A, Density, rr, alpha, beta, &
                 iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, &
      rr, &      
      alpha, &
      beta
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( StorageForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Geometry

    G => A % Geometry ( )

    associate (  R => G % Value ( :, G % CENTER_U ( 1 ) ) )


    !-- Source

    Source => Source_ASC % Storage ( )

    associate ( D => Source % Value ( :, iVariable ) )

      D = Density &
          / ( ( R / rr ) ** alpha &
              * ( 1.0_KDR + ( R / rr ) ) ** ( beta - alpha ) )

    end associate !-- D


    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ) )

    if ( beta == 4.0_KDR ) then
      if ( alpha == 1.0_KDR ) then
        Phi = - Density * rr ** 2 * 1.0_KDR / ( 2 * ( 1 + R / rr ) )
      else if ( alpha == 2.0_KDR ) then
        Phi = - Density * rr ** 2 * log ( 1 + rr / max ( R, tiny ( 0.0_KDR ) ) )
      else
        call Show ( '**** Analytic Solution Not Implemented ****' )
      end if
    else if ( beta == 3.0_KDR ) then
      Phi = - Density * rr ** 2 * log ( 1 + R / rr ) * rr / max ( R, tiny ( 0.0_KDR ) )
    else
      call Show ( '**** Analytic Solution Not Implemented ****' )
    end if

    end associate !-- Phi


    !-- Cleanup

    end associate !-- R, etc.

    nullify ( G, Source, Reference )

  end subroutine SetTwoPowerDensityKernal


end module TwoPowerDensity_Form
