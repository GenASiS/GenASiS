module SetCondensedSphere_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetCondensedSphere

contains


  subroutine SetCondensedSphere &
               ( Source_ASC, Reference_ASC, A, Density, RadiusDensity, &
                 RadiusCondensed, iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, &
      RadiusDensity, &
      RadiusCondensed
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
      (  R => G % Value ( :, G % CENTER ( 1 ) ), &
        dR => G % Value ( :, G % WIDTH ( 1 ) ) )


    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D => Source % Value ( :, iVariable ), &
        R_In  => R - 0.5_KDR * dR, &
        R_Out => R + 0.5_KDR * dR, &
        RD    => RadiusDensity, &
        RC    => RadiusCondensed )
    D = 0.0_KDR
    where ( R <= RD )
      D = Density / ( 1.0_KDR + ( R / RC ) ** 2 )
    end where
    end associate !-- R_In, etc.


    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        RD    => RadiusDensity, &
        RC    => RadiusCondensed, &
        Pi   =>  CONSTANT % PI )

    where ( R < RadiusDensity )
      Phi  =  &
             Density * RC ** 2 &
             * ( atan ( R / RC ) / ( R / RC ) &
                 + log ( ( 1.0_KDR + ( R / RC ) ** 2 ) &
                          / ( 1.0_KDR + ( RD / RC ) ** 2) ) / 2 &
                 - 1.0_KDR )
                       
    elsewhere
      Phi  =  - Density * RC ** 3 / R * ( RD / RC - atan ( RD / RC ) )
    end where
    end associate !-- R


    !-- Cleanup

    end associate !-- R, etc.

    nullify ( G, Source, Reference )

  end subroutine SetCondensedSphere


end module SetCondensedSphere_Command
