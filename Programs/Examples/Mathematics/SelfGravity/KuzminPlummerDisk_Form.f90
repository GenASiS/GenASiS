!! Parameters for KuzminPlummerDisk_1 taken from 
!!   Muller Stienmetz 1995 for comparison.
!! Parameters for KuzminPlummerDisk_2 & KuzminPlummerDisk_3 
!!   taken from Binney and Tremaine "Galactic Dynamics"
!! Parameters for KuzminPlummerDisk_4 made so that most of the 
!!  density profile is inside the domain

!! Analytic solutions taken from Binney and Tremaine "Galactic Dynamics" second edition 
!!   Chapter 2 Section 2.3.1 pg 72-73, Eqns. 2.69a & 2.69b.

!! Density distribution ( Eqn. 2.69b in BT ) from Miyamoto & Nagai 1975

! Constructs DD % N_equations Kuzmin-Plummer density profiles, 
!   making the last profile contained mostly within the computational domain.

module KuzminPlummerDisk_Form

  use Basics
  use Mathematics
  use DensityDistribution_Template
  
  implicit none

  type, public, extends ( DensityDistributionTemplate ) :: &
    KuzminPlummerForm
  contains
    procedure, public, pass :: &
      SetKuzminPlummer
  end type KuzminPlummerForm

    private :: &
      SetKuzminPlummerKernal

contains


  subroutine SetKuzminPlummer ( KP, TotalMass, KPa, ratio )
    class ( KuzminPlummerForm ), intent ( inout ) :: &
      KP
    real ( KDR ), intent ( in ) :: &
      TotalMass
    real ( KDR ), dimension ( KP % N_Equations ) :: &
      KPa, &
      ratio

    integer ( KDI ) :: &
      iE !-- iEquation
    real ( KDR ) :: &
      ProlateMajor
    real ( KDR ), dimension ( KP % N_Equations ) :: &
      b

    KPa ( KP % N_Equations ) = 0.4_KDR

    b = ratio * KPa
    b ( KP % N_Equations ) = 0.1_KDR

    associate ( A  => KP % Atlas )
    associate ( SA => KP % Source )
    associate ( RA => KP % Reference )

    do iE = 1, KP % N_Equations
      call SetKuzminPlummerKernal &
             ( SA, RA, A, TotalMass, KPa ( iE ), b ( iE ), iE )
    end do

    end associate !-- RA
    end associate !-- SA
    end associate !-- A

  end subroutine SetKuzminPlummer


  subroutine SetKuzminPlummerKernal &
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

    D = C_b ** 2 * M & 
          * ( C_a * rho_sq + ( C_a + 3 * sqrt ( Z_sq + C_b ** 2 ) )  &
                           * ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 ) &
          / ( ( rho_sq + ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 ) ** 2.5_KDR &
              * ( Z_sq + C_b ** 2 ) ** 1.5_KDR )

    !-- Reference

    Reference => Reference_ASC % Storage ( )

    associate &
      ( Phi  =>  Reference % Value ( :, iVariable ), &
        G_C  =>  1.0_KDR  ) 
    Phi = - G_C * M & 
            / sqrt ( rho_sq + ( C_a + sqrt ( Z_sq + C_b ** 2 ) ) ** 2 )

    end associate !-- Phi, etc
    end associate !-- R_In, etc.
    end associate !-- R, etc.

    deallocate ( rho_sq, Z_sq )
    nullify ( G, Source, Reference )

  end subroutine SetKuzminPlummerKernal


end module KuzminPlummerDisk_Form
