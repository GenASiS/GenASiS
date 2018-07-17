module SetSatohDisk_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetSatohDisk

contains

  subroutine SetSatohDisk &
               ( Source_ASC, Reference_ASC, A, M, C_a, C_b, iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      M, &
      C_a, &
      C_b
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( StorageForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      rho_sq, & !-- X^2 + Y^2 
      Z_sq, &
      S
    integer ( KDI ) :: &
      i

    !-- Geometry

    G => A % Geometry ( )

    associate &
      (  R    => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ) )

    allocate &
      ( rho_sq ( size ( R ) ) , &
        Z_sq   ( size ( R ) ), &
        S      ( size ( R ) ) )
         

    rho_sq = ( R * sin ( Theta ) ) ** 2
    Z_sq   = ( R * cos ( Theta ) ) ** 2

    S = sqrt ( rho_sq + Z_sq + C_a * ( C_a + 2 * sqrt ( Z_sq + C_b ** 2 ) ) )

    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D     => Source % Value ( :, iVariable ) )

    D = C_a * C_b ** 2 * M / ( S ** 3 * ( Z_sq + C_b ** 2 ) ) & 
          * ( 1.0_KDR / sqrt ( Z_sq + C_b ** 2 ) + 3.0_KDR / C_a  &
              * ( 1.0_KDR - ( rho_sq + Z_sq ) / S ** 2  ) )

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI )

    Phi = - M  / S

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq )
    nullify ( G, Source, Reference )

  end subroutine SetSatohDisk


end module SetSatohDisk_Command
