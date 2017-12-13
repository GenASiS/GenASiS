module SetHomogeneousSphere_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetHomogeneousSphere

contains


  subroutine SetHomogeneousSphere &
               ( Source_ASC, Reference_ASC, A, Density, RadiusDensity, &
                 iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, RadiusDensity
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( VariableGroupForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Geometry

    G => A % Geometry ( )

    associate &
      (    R => G % Value ( :, G % CENTER ( 1 ) ), &
        dR_L => G % Value ( :, G % WIDTH_LEFT ( 1 ) ), &
        dR_R => G % Value ( :, G % WIDTH_RIGHT ( 1 ) ) )


    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D => Source % Value ( :, iVariable ), &
        R_In  => R - dR_L, &
        R_Out => R + dR_R, &
        RD    => RadiusDensity )
    D = 0.0_KDR
    where ( R_Out <= RD )
      D = Density
    end where
    where ( R_In < RD .and. R_Out > RD )
      D = Density * ( RD - R_In ) / ( dR_L + dR_R )
    end where
    end associate !-- R_In, etc.


    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI )
    where ( R < RadiusDensity )
      Phi  =  1.0_KDR / 6.0_KDR  *  Density  *  R ** 2  &
              -  1.0_KDR / 2.0_KDR  *  Density  *  RadiusDensity ** 2
    elsewhere
      Phi  =  - 1.0_KDR / 3.0_KDR  *  Density  *  RadiusDensity ** 3  /  R
    end where
    end associate !-- R


    !-- Cleanup

    end associate !-- R, etc.

    nullify ( G, Source, Reference )

  end subroutine SetHomogeneousSphere


end module SetHomogeneousSphere_Command
