module SetHomogeneousSpheroid_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetHomogeneousSpheroid

contains

  subroutine SetHomogeneousSpheroid &
               ( Source_ASC, Reference_ASC, A, Density, SemiMajor, SemiMinor, &
                 iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, &
      SemiMajor, &
      SemiMinor
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( VariableGroupForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G
    real ( KDR ), allocatable, dimension ( : ) :: &
      rho_sq, & !-- X^2 + Y^2 
      Z_sq, &
      l, &
      h
    real ( KDR ) :: &
      e, &
      C_A, &
      C_B
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
        Z_sq   ( size ( R ) ), &
        l   ( size ( R ) ), &
        h   ( size ( R ) ) )
         

    rho_sq = ( R * sin ( Theta ) ) ** 2
    Z_sq   = ( R * cos ( Theta ) ) ** 2

    !-- Source

    Source => Source_ASC % Storage ( )

    associate &
      ( D     => Source % Value ( :, iVariable ), &
        R_In  => R - 0.5_KDR * dR, &
        R_Out => R + 0.5_KDR * dR, &
        a     => SemiMajor, &
        b     => SemiMinor )

    D = 0.0_KDR

call Show ( 'Before density' )
    where ( rho_sq /a ** 2 + Z_sq / b ** 2 <= 1.0_KDR )
      D = Density
    end where

    ! where ( R_In * sin ( Theta ) < a .and. R_In * cos ( Theta ) < b )
    !   where ( R_Out * sin ( Theta ) > a .and. R_Out * cos ( Theta ) > b )
    !     D = Density * ( RD - R_In ) / dR
    !   end where
    ! end where


    !-- Reference

    Reference => Reference_ASC % Storage ( )

call Show ( b, 'b' )

    e = sqrt ( 1.0_KDR - ( b / a ) ** 2 )

call Show ( e, 'e' )
    C_A = sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e ) &
          - ( 1.0_KDR - e ** 2 ) / e ** 2

    C_B = 2.0_KDR / e ** 2 - 2 * sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e )

    where ( rho_sq /a ** 2 + Z_sq / b ** 2 > 1.0_KDR )
      l = 0.5_KDR * ( ( rho_sq + Z_sq - a ** 2 - b **2  ) &
            + sqrt ( 2.0_KDR * rho_sq * ( - a ** 2 + b ** 2 + Z_sq ) &
                     + ( a ** 2 - b ** 2 + Z_sq ) ** 2 + rho_sq ** 2 ) )
         ! + sqrt ( ( a ** 2 - b ** 2 ) ** 2 + ( rho + Z ) ** 2 &
         !          + 2.0_KDR * a ** 2 * ( Z - rho ) &
         !          - 2.0_KDR * b ** 2 * ( Z - rho ) ) )


      where ( l < 0.0_KDR )
        l = 0.5_KDR * ( ( rho_sq + Z_sq - a ** 2 - b **2  ) &
            - sqrt ( 2.0_KDR * rho_sq * ( - a ** 2 + b ** 2 + Z_sq ) &
                     + ( a ** 2 - b ** 2 + Z_sq ) ** 2 + rho_sq ** 2 ) )
         ! - sqrt ( ( a ** 2 - b ** 2 ) ** 2 + ( rho + Z ) ** 2 &
         !          + 2.0_KDR * a ** 2 * ( Z - rho ) &
         !          - 2.0_KDR * b ** 2 * ( Z - rho ) ) )
      end where

      h = a * e / sqrt ( b ** 2 + l )
    end where


call Show ( 'before phi' )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        Pi   =>  CONSTANT % PI )
    where ( rho_sq /a ** 2 + Z_sq / b ** 2 <= 1.0_KDR  )
      Phi  =   - Density / 4.0_KDR * &
                   ( C_A * ( 2.0_KDR * a ** 2 - rho_sq ) + C_B * ( b ** 2 - Z_sq ) )
    elsewhere
      Phi  =   -  Density / 2.0_KDR * a * b / e &
                  * ( atan ( h ) &
                        - 1 / ( 2 * a ** 2 * e ** 2 )&
                            * ( rho_sq * ( atan ( h ) - h / ( 1 + h ** 2 ) ) &
                                + 2 * Z_sq * ( h - atan ( h ) ) ) )
    end where

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq, l, h )
    nullify ( G, Source, Reference )

  end subroutine SetHomogeneousSpheroid


end module SetHomogeneousSpheroid_Command
