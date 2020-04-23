!! Parameters for spheroids taken from Muller Stienmetz 1995 for comparison.
!! Analytic solutions taken from Binney and Tremaine "Galactic Dynamics" second edition 
!!   Chapter 2 Section 2.5.2 Table 2.1 pg 90.

! Constructs DD % N_equations homogeneous spheroid density profiles, 
!   making the last profile prolate.

module HomogeneousSpheroid_Form

  use Basics
  use Mathematics
  use DensityDistribution_Template
  
  implicit none

  type, public, extends ( DensityDistributionTemplate ) :: &
    HomogeneousSpheroidForm
  contains
    procedure, public, pass :: &
      SetHomogeneousSpheroid
  end type HomogeneousSpheroidForm

    private :: &
      SetHomogeneousSpheroidKernal

contains


  subroutine SetHomogeneousSpheroid ( HS, Density, SemiMajor, Eccentricity )
    class ( HomogeneousSpheroidForm ), intent ( inout ) :: &
      HS
    real ( KDR ), intent ( in ) :: &
      Density
    real ( KDR ), dimension ( HS % N_Equations ) :: &
      Eccentricity, &
      SemiMajor

    integer ( KDI ) :: &
      iE !-- iEquation
    real ( KDR ) :: &
      ProlateMajor
    real ( KDR ), dimension ( HS % N_Equations ) :: &
      SemiMinor

    SemiMinor = sqrt ( 1.0_KDR - Eccentricity ** 2 ) * SemiMajor

    !-- Prolate Spheroid 
    ProlateMajor                   = SemiMajor ( HS % N_Equations )
    SemiMajor ( HS % N_Equations ) = SemiMinor ( HS % N_Equations )
    SemiMinor ( HS % N_Equations ) = ProlateMajor

    associate ( A => HS % Atlas )
    associate ( SA => HS % Source )
    associate ( RA => HS % Reference )

    do iE = 1, HS % N_Equations
      call SetHomogeneousSpheroidKernal &
             ( SA, RA, A, Density, SemiMajor ( iE), SemiMinor ( iE ), iE )
    end do

    end associate !-- RA
    end associate !-- SA
    end associate !-- A

  end subroutine SetHomogeneousSpheroid 


  subroutine SetHomogeneousSpheroidKernal &
               ( Source_ASC, Reference_ASC, A, Density, a_1, a_3, &
                 iVariable )

    class ( Storage_ASC_Form ), intent ( inout ) :: &
      Source_ASC, &
      Reference_ASC
    class ( Atlas_SC_Form ), intent ( in ) :: &
      A
    real ( KDR ), intent ( in ) :: &
      Density, &
      a_1, &
      a_3
    integer ( KDI ), intent ( in ) :: &
      iVariable

    type ( StorageForm ), pointer :: &
      Source, &
      Reference
    class ( GeometryFlatForm ), pointer :: &
      G
    integer ( KDI ) :: &
      i, &
      iC, &       !-- iCell
      iS, jS, kS  !-- iSubcell
    integer ( KDI ), dimension ( 3 ) :: &
      nSubcells
    real ( KDR ) :: &
      dVS, &  !-- dVolumeSubcell
      VS, &   !-- VolumeSubcell
      e, &
      a_1_sq, &
      a_3_sq, &
      C_I, &
      C_A, &
      C_B, &
      rho_sq_II, &
      rho_sq_IO, &
      rho_sq_OI, &
      rho_sq_OO, &
      Z_sq_II, &
      Z_sq_IO, &
      Z_sq_OI, &
      Z_sq_OO, &
      rho_S_sq, &
      Z_s_sq, &
      Pi, &
      FourPi
    real ( KDR ), dimension ( 3 ) :: &
      X_I, &
      X_O, &
      dXS, &  !-- dX_Subcell
      XS,  &  !--  X_Subcell
      X, &
      dX_L, &
      dX_R
    real ( KDR ), dimension ( : ), allocatable :: &
      BVF, &  !-- VolumeFraction
      rho_sq, & !-- X^2 + Y^2 
      Z_sq, &
      l, &
      e_vec, &
      C_I_vec, &
      C_A_vec, &
      C_B_vec

        Pi  =  CONSTANT % PI
    FourPi  =  4.0_KDR * Pi

    !-- Geometry

    G => A % Geometry ( )

    select type ( C => A % Chart )
    class is ( Chart_SLD_Form )
      
    associate ( dV => G % Value ( :, G % VOLUME ) )

    allocate ( BVF ( size ( dV ) ) )
    call Clear ( BVF )

    nSubCells = 1
    nSubcells ( : C % nDimensions ) = 20

    a_1_sq  =  a_1 ** 2
    a_3_sq  =  a_3 ** 2

    do iC = 1, size ( dV )

      if ( .not. C % IsProperCell ( iC ) ) cycle

         X  =  G % Value ( iC, G % CENTER_U ( 1 ) : G % CENTER_U ( 3 ) )
      dX_L  =  G % Value ( iC, G % WIDTH_LEFT_U ( 1 ) &
                               : G % WIDTH_LEFT_U ( 3 ) )
      dX_R  =  G % Value ( iC, G % WIDTH_RIGHT_U ( 1 ) &
                               : G % WIDTH_RIGHT_U ( 3 ) )

      X_I  =  X - dX_L
      X_O  =  X + dX_R

      rho_sq_II  =  ( X_I ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
      rho_sq_IO  =  ( X_O ( 1 ) * sin ( X_I ( 2 ) ) ) ** 2
      rho_sq_OI  =  ( X_I ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2
      rho_sq_OO  =  ( X_O ( 1 ) * sin ( X_O ( 2 ) ) ) ** 2

      Z_sq_II  =  ( X_I ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
      Z_sq_IO  =  ( X_O ( 1 ) * cos ( X_I ( 2 ) ) ) ** 2
      Z_sq_OI  =  ( X_I ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2
      Z_sq_OO  =  ( X_O ( 1 ) * cos ( X_O ( 2 ) ) ) ** 2

      if (       rho_sq_II / a_1_sq  +  Z_sq_II / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_IO / a_1_sq  +  Z_sq_IO / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_OI / a_1_sq  +  Z_sq_OI / a_3_sq  <=  1.0_KDR &
           .and. rho_sq_OO / a_1_sq  +  Z_sq_OO / a_3_sq  <=  1.0_KDR ) &
      then 
        BVF ( iC ) = 1.0_KDR
        cycle
      end if

      if (       rho_sq_II / a_1_sq  +  Z_sq_II / a_3_sq  >  1.0_KDR &
           .and. rho_sq_IO / a_1_sq  +  Z_sq_IO / a_3_sq  >  1.0_KDR &
           .and. rho_sq_OI / a_1_sq  +  Z_sq_OI / a_3_sq  >  1.0_KDR &
           .and. rho_sq_OO / a_1_sq  +  Z_sq_OO / a_3_sq  >  1.0_KDR ) &
      then 
        BVF ( iC ) = 0.0_KDR
        cycle
      end if

      dXS  =  ( X_O  -  X_I ) / nSubcells

      VS = 0.0_KDR
      do kS = 1, nSubcells ( 3 )
        do jS = 1, nSubcells ( 2 )
          do iS = 1, nSubcells ( 1 )
            XS  =  X_I  +  ( [ iS, jS, kS ] - 0.5_KDR ) * dXS
            rho_S_sq = ( XS ( 1 ) * sin ( XS ( 2 ) ) ) ** 2
            Z_S_sq = ( XS ( 1 ) * cos ( XS ( 2 ) ) ) ** 2
            select case ( C % nDimensions )
            case ( 2 )
              dVS = 2 * Pi * XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                      * dXS ( 1 ) * dXS ( 2 )
            case ( 3 )
              dVS = XS ( 1 ) ** 2  * sin ( XS ( 2 ) ) &
                     * dXS ( 1 ) * dXS ( 2 ) * dXS ( 3 )
            end select 
            VS = VS + dVS
            if ( rho_S_sq / a_1_sq + Z_S_sq / a_3_sq <= 1.0_KDR ) &
              BVF ( iC ) = BVF ( iC ) + dVS
          end do !-- iS
        end do !-- jS
      end do !-- kS
      BVF ( iC ) = BVF ( iC ) / VS

    end do !-- iC

    end associate !-- VJ
!    end select !-- C

    associate &
      (  R    => G % Value ( :, G % CENTER_U ( 1 ) ), &
        Theta => G % Value ( :, G % CENTER_U ( 2 ) ) )

    allocate &
      ( rho_sq  ( size ( R ) ) , &
        Z_sq    ( size ( R ) ), &
        l       ( size ( R ) ), &
        C_I_vec ( size ( R ) ) )
         

    rho_sq = ( R * sin ( Theta ) ) ** 2
    Z_sq   = ( R * cos ( Theta ) ) ** 2

    !-- Source

    Source => Source_ASC % Storage ( )

    associate ( D => Source % Value ( :, iVariable ) )

    D = Density * BVF * FourPi

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    e = sqrt ( 1.0_KDR &
               - min ( ( a_3 / a_1 ), ( a_1 / a_3 ) )  ** 2 )

    if ( a_3 < a_1 ) then
      C_I = 2 * sqrt ( 1.0_KDR - e ** 2 ) / e * asin ( e )
      C_A = sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e ) &
            - ( 1.0_KDR - e ** 2 ) / e ** 2

      C_B = 2.0_KDR / e ** 2 &
            - 2 * sqrt ( 1.0_KDR - e ** 2 ) / e ** 3 * asin ( e )
    else
      C_I  = 1.0_KDR / e * log ( ( 1 + e ) / ( 1 - e ) )
      C_A = 1.0_KDR / e ** 2 &
            - ( 1.0_KDR - e ** 2 ) / ( 2 * e ** 3 ) &
               * log ( ( 1 + e ) / ( 1 - e ) )

      C_B = ( 1.0_KDR - e ** 2 ) / e ** 3  * log ( ( 1 + e ) / ( 1 - e ) ) &
            - 2 * ( 1.0_KDR - e ** 2 ) / e ** 2
    end if

    associate ( Phi  =>  Reference % Value ( :, iVariable ) )

    where ( rho_sq / a_1_sq + Z_sq / a_3_sq <= 1.0_KDR )
      Phi  =   - Density * Pi &
                   * ( C_I * a_1_sq  - C_A * rho_sq - C_B * Z_sq )
    end where

    l = 0.0_KDR

    where ( rho_sq / a_1_sq + Z_sq / a_3_sq > 1.0_KDR )
      l = 0.5_KDR * ( ( rho_sq + Z_sq - a_1_sq - a_3 **2  ) &
            + sqrt ( ( a_1_sq + a_3_sq - Z_sq - rho_sq) ** 2 &
                     - 4.0_KDR * ( a_1_sq * a_3_sq - rho_sq * a_3_sq &
                                   - Z_sq * a_1_sq ) ) )

    end where

    where ( l < 0.0_KDR )
      l = 0.5_KDR * ( ( rho_sq + Z_sq - a_1_sq - a_3 **2  ) &
          - sqrt ( ( a_1_sq + a_3_sq - Z_sq - rho_sq) ** 2 &
                   - 4.0_KDR * ( a_1_sq * a_3_sq - rho_sq * a_3_sq &
                                 - Z_sq * a_1_sq ) ) )
    end where

    if ( a_3 < a_1 ) then
      C_I_vec = Pi / sqrt ( a_1_sq - a_3_sq ) &
            - 2.0_KDR / sqrt ( a_1_sq - a_3_sq ) &
            * atan ( sqrt ( ( a_3_sq + l ) / ( a_1_sq - a_3_sq ) ) )
    else
      C_I_vec = -1.0_KDR / sqrt ( a_3_sq - a_1_sq ) &
                * log ( ( sqrt ( a_3_sq + l ) &
                           - sqrt ( a_3_sq - a_1_sq ) ) ** 2 &
                             / ( a_1_sq + l ) )
    end if

    where (  rho_sq / a_1_sq + Z_sq / a_3_sq > 1.0_KDR )
      Phi  = - Density * Pi * a_1_sq * a_3 &
               * ( ( 1.0_KDR + rho_sq / ( 2 * ( a_3_sq - a_1_sq ) ) &
                   - Z_sq / ( a_3_sq - a_1_sq ) ) * C_I_vec &
             - rho_sq * sqrt ( a_3_sq + l ) &
               / ( ( a_3_sq - a_1_sq ) * ( a_1_sq + l ) ) &
             - Z_sq * ( 2.0_KDR &
                        / ( ( a_1_sq + l ) * sqrt ( a_3_sq + l ) ) &
                        - 2.0_KDR  * sqrt ( a_3_sq + l ) &
                          / ( ( a_3_sq - a_1_sq ) * ( a_1_sq + l ) ) ))
    end where

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq, l, C_I_vec )
    nullify ( G, Source, Reference )
    end select !-- C

  end subroutine SetHomogeneousSpheroidKernal


end module HomogeneousSpheroid_Form
