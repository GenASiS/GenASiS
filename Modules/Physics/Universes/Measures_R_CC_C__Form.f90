module Measures_R_CC_C__Form

  !-- Measures_Radiation_CentralCore_Collected__Form

  use Basics
  use Mathematics
  use Fluids
  use Radiations
  use Measures_F_CC__Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_MEASURES_NM = 2

  type, public, extends ( Measures_F_CC_Form ) :: Measures_R_CC_C_Form
    integer ( KDI ) :: &
      N_MEASURES_NM = N_MEASURES_NM
    integer ( KDI ) :: &
      nRadiations, &
      nMeasures_R
    real ( KDR ), dimension ( : ), allocatable :: &
      Luminosity, &
      EnergyAverage
    type ( FieldSet_BM_Pointer ), dimension ( : ), allocatable :: &
      Radiation_SA_1D
    class ( Units_R_Form ), pointer :: &
      Units_R => null ( )
  contains
    procedure, private, pass :: &
      Initialize_R_CC_C
    generic, public :: &
      Initialize => Initialize_R_CC_C
    final :: &
      Finalize
  end type Measures_R_CC_C_Form


contains


  subroutine Initialize_R_CC_C &
               ( M, R_SA_1D, F_SA, G_SA, A_SA, Units_R, Units_F )

    class ( Measures_R_CC_C_Form ), intent ( inout ) :: &
      M
    class ( FieldSet_BM_Pointer ), dimension ( : ), intent ( in ) :: &
      R_SA_1D
    class ( FieldSet_BM_Form ), intent ( in ), target :: &
      F_SA, &
      G_SA
    class ( Atlas_H_Form ), intent ( in ), target :: &
      A_SA
    class ( Units_R_Form ), intent ( in ), target :: &
      Units_R
    class ( Units_F_Form ), intent ( in ), target :: &
      Units_F

    integer ( KDI ) :: &
      oV, &
      iR

    M % nRadiations  =  size ( R_SA_1D )
    M % nMeasures_R  =  M % nRadiations  *  M % N_MEASURES_NM

    call M % Measures_F_CC_Form % Initialize &
           ( F_SA, G_SA, A_SA, Units_F, nMeasuresAddOption = M % nMeasures_R )

    M % Units_R  =>  Units_R

    associate &
      ( UR   =>  M % Units_R, &
        nMF  =>  M % nMeasures_F, &
        nMN  =>  M % N_MEASURES_NM, &
        nR   =>  M % nRadiations )

    allocate ( M % Radiation_SA_1D ( nR ) )
    do iR  =  1, nR

      M % Radiation_SA_1D ( iR ) % Pointer  =>  R_SA_1D ( iR ) % Pointer
      associate ( R  =>  M % Radiation_SA_1D ( iR ) % Pointer)

      oV  =  nMF  +  ( iR - 1 ) * nMN

      M % Name ( oV + 1 )  =  'Luminosity_' // trim ( R % Name )
      M % Name ( oV + 2 )  =  'EnergyAve_'  // trim ( R % Name )

      M % Unit ( oV + 1 )  =  UR % Luminosity
      M % Unit ( oV + 2 )  =  UR % EnergyAverage
  
      end associate !-- R

    end do !-- iR

    end associate !-- UR, etc.

  end subroutine Initialize_R_CC_C


  impure elemental subroutine Finalize ( M )

    type ( Measures_R_CC_C_Form ), intent ( inout ) :: &
      M

    if ( allocated ( M % Radiation_SA_1D ) ) &
      deallocate ( M % Radiation_SA_1D )
    if ( allocated ( M % EnergyAverage ) ) &
      deallocate ( M % EnergyAverage )
    if ( allocated ( M % Luminosity ) ) &
      deallocate ( M % Luminosity )

    nullify ( M % Units_R )

  end subroutine Finalize


end module Measures_R_CC_C__Form
