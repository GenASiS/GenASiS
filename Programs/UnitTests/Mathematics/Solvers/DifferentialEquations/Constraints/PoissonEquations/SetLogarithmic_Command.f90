module SetLogarithmic_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetLogarithmic

contains

  subroutine SetLogarithmic &
               ( Source_ASC, Reference_ASC, A, v0, rho_c, q, iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      v0, &
      rho_c, &
      q
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

    D = ( v0 / q ) ** 2 & 
          * ( (  2 * q ** 2 + 1.0_KDR ) * rho_c ** 2 &
                + rho_sq + ( 2.0_KDR - q ** ( - 2 ) ) * Z_sq ) &
          / ( rho_c ** 2 + rho_sq + Z_sq / q ** 2 ) ** 2

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI )
    
    !Phi = 1.0_KDR
    
   ! where ( D >= 0.0_KDR )
      Phi = v0 ** 2 * log ( rho_c ** 2 + rho_sq + Z_sq / ( q ** 2 ) ) / 2
   ! end where

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq )
    nullify ( G, Source, Reference )

  end subroutine SetLogarithmic


end module SetLogarithmic_Command
