module SetKuzminPlummerDisk_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetKuzminPlummerDisk

contains

  subroutine SetKuzminPlummerDisk &
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

    type ( VariableGroupForm ), pointer :: &
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
      (  R    => G % Value ( :, G % CENTER ( 1 ) ), &
        dR    => G % Value ( :, G % WIDTH ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER ( 2 ) ) )

    allocate &
      ( rho_sq ( size ( R ) ) , &
        Z_sq   ( size ( R ) ) )
         

    rho_sq = ( R * sin ( Theta ) ) ** 2
    Z_sq   = ( R * cos ( Theta ) ) ** 2

    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D     => Source % Value ( :, iVariable ) )

    D = C_b ** 2 * M & 
          * ( C_a * rho_sq + ( C_a + 3 * sqrt ( Z_sq + C_b ** 2 ) )  &
                           * ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 ) &
          / ( ( rho_sq + ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 ) ** 2.5_KDR &
              * ( Z_sq + C_b ** 2 ) ** 1.5_KDR )

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI, &
        G_C  =>  1.0_KDR  ) !CONSTANT % GRAVITATIONAL )

    Phi = - G_C * M  / sqrt ( rho_sq + ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 )

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq )
    nullify ( G, Source, Reference )

  end subroutine SetKuzminPlummerDisk


end module SetKuzminPlummerDisk_Command
