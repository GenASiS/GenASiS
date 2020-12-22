module SetHomogeneousSphere_Command

  use Basics
  use Manifolds

  implicit none
  private

  public :: &
    SetHomogeneousSphere

    private :: &
      SetDensityKernel

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

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    type ( StorageForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G


    !-- Geometry

    G => A % Geometry ( )

    associate &
      (    R => G % Value ( :, G % CENTER_U ( 1 ) ), &
        dR_L => G % Value ( :, G % WIDTH_LEFT_U ( 1 ) ), &
        dR_R => G % Value ( :, G % WIDTH_RIGHT_U ( 1 ) ) )

    nValues  =  size ( R )


    !-- Source

    Source => Source_ASC % Storage ( )

    associate ( D => Source % Value ( :, iVariable ) )
    call SetDensityKernel ( R, dR_L, dR_R, RadiusDensity, Density, nValues, D )
    end associate !-- D


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


  subroutine SetDensityKernel ( R, dR_L, dR_R, RD, Density, nValues, D )

    real ( KDR ), dimension ( : ), intent ( in ) :: &
      R, &
      dR_L, dR_R
    real ( KDR ), intent ( in ) :: &
      RD, &
      Density
    integer ( KDI ), intent ( in ) :: &
      nValues
    real ( KDR ), dimension ( : ), intent ( out ) :: &
      D

    integer ( KDI ) :: &
      iV
    real ( KDR ) :: &
      R_In, R_Out

    do iV  =  1, nValues
      R_In   =  R ( iV )  -  dR_L ( iV )
      R_Out  =  R ( iV )  +  dR_R ( iV )
      if ( R_Out  <=  RD ) then
        D ( iV )  =  Density
      else if ( R_In  <  RD .and. R_Out  >  RD ) then
        D ( iV )  =  Density * ( RD ** 3  -  R_In ** 3 ) &
                     / ( R_Out ** 3  -  R_In ** 3 )
      else
        D ( iV )  =  0.0_KDR
      end if
    end do !-- iV

  end subroutine SetDensityKernel


end module SetHomogeneousSphere_Command
