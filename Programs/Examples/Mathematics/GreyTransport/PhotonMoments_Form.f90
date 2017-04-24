module PhotonMoments_Form

  use Basics
  use RadiationMoments_Form

  type, public, extends ( RadiationMomentsForm ) :: PhotonMomentsForm
  contains
    procedure, public, pass ( RM ) :: &
      ComputeSpectralParameters      
  end type PhotonMomentsForm

contains


  subroutine ComputeSpectralParameters &
               ( T, Eta, E_Ave, F_Ave, J_EQ, RM, J, N, T_EQ, Eta_EQ )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T, &
      Eta, &
      E_Ave, &
      F_Ave, &
      J_EQ
    class ( PhotonMomentsForm ), intent ( in ) :: &
      RM
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      N, &
      T_EQ, &
      Eta_EQ

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      a

    nValues = size ( T )

    a = CONSTANT % RADIATION

    !$OMP parallel do private ( iV )
    do iV = 1, nValues
      T  ( iV )   =  ( J ( iV )  /  a ) ** ( 0.25_KDR )
      Eta ( iV )  =  0.0_KDR
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeSpectralParameters


end module PhotonMoments_Form
