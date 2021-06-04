!! Analytic solutions taken from Stone and Norman 1992
!!  Test 2 under section 5.6

! Constructs DD % N_equations condensed sphere density profiles

module CondensedSphere_Form

  use Basics
  use Mathematics
  use DensityDistribution_Template
  
  implicit none

  type, public, extends ( DensityDistributionTemplate ) :: &
    CondensedSphereForm
  contains
    procedure, public, pass :: &
      SetCondensedSphere
  end type CondensedSphereForm

    private :: &
      SetCondensedSphereKernal

contains


  subroutine SetCondensedSphere ( CS, Density, SR, rr )
    class ( CondensedSphereForm ), intent ( inout ) :: &
      CS
    real ( KDR ), intent ( in ) :: &
      Density
    real ( KDR ), dimension ( CS % N_Equations ) :: &
      SR, &       !-- sphere radius
      rr          !-- reference radius

    integer ( KDI ) :: &
      iE !-- iEquation

    associate ( A  => CS % Atlas )
    associate ( SA => CS % Source )
    associate ( RA => CS % Reference )

    do iE = 1, CS % N_Equations
      call SetCondensedSphereKernal &
             ( SA, RA, A, Density, SR ( iE ), rr ( iE ), iE )
    end do

    end associate !-- RA
    end associate !-- SA
    end associate !-- A

  end subroutine SetCondensedSphere


  subroutine SetCondensedSphereKernal &
               ( Source_ASC, Reference_ASC, A, Density, SR, rr, iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, &
      SR, &
      rr
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

    associate ( D => Source % Value ( :, iVariable ), &
                FourPi => 4 * CONSTANT % PI )

    where ( R <= SR )
      D = Density / ( 1.0_KDR + ( R / rr ) ** 2 ) * FourPi
    elsewhere 
      D = 0.0_KDR
    end where

    end associate !-- D


    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        FourPi => 4 * CONSTANT % PI )

    where ( R <= SR ) 
      Phi = &
            Density * FourPi * rr ** 2 * ( atan ( R / rr ) / ( R / rr ) &
              + log ( ( 1.0_KDR + ( R / rr ) ** 2 ) &
                        / ( 1.0_KDR + ( SR / rr ) ** 2 ) ) / 2 - 1.0_KDR )
    elsewhere
      Phi = - Density * FourPi * rr ** 3 / R * ( SR / rr - atan ( SR / rr ) )
    end where 

    end associate !-- Phi

    !-- Cleanup

    end associate !-- R

    nullify ( G, Source, Reference )

  end subroutine SetCondensedSphereKernal


end module CondensedSphere_Form
