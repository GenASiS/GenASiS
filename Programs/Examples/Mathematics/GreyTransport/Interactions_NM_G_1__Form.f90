module Interactions_NM_G_1__Form

  !-- Interactions_NeutrinoMoments_Grey_1__Form
  
  use Basics
  use Mathematics
  use Fluid_P_MHN__Form
  use NeutrinoMoments_Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_NM_G_1_Form
    real ( KDR ) :: &
      SpecificOpacity = 0.0_KDR
    class ( Fluid_P_MHN_Form ), pointer :: &
      Fluid => null ( )
    class ( NeutrinoMomentsForm ), pointer :: &
      Neutrino_E_Nu          => null ( ), &
      Neutrino_E_NuBar       => null ( ), &
      Neutrino_MuTau_NuNuBar => null ( )
  contains
    procedure, public, pass :: &
      Compute
    procedure, public, pass ( I ) :: &
      ComputeDegeneracyParameter_EQ
  end type Interactions_NM_G_1_Form

    private :: &
      ComputeDegeneracyParameter_EQ_Kernel

contains


  subroutine Compute ( I )

    class ( Interactions_NM_G_1_Form ), intent ( inout ) :: &
      I

  end subroutine Compute


  subroutine ComputeDegeneracyParameter_EQ ( Eta_EQ, I, C )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Eta_EQ
    class ( Interactions_NM_G_1_Form ), intent ( in ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      C

    associate ( F => I % Fluid )
    associate &
      (  Mu_E =>  F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
        Mu_NP =>  F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
            T =>  F % Value ( :, F % TEMPERATURE ) )

    select case ( trim ( C % Type ) )
    case ( 'NEUTRINO_E_NU' )
      call ComputeDegeneracyParameter_EQ_Kernel &
             ( Eta_EQ, Mu_E, Mu_NP, T, Sign = 1.0_KDR )
    case ( 'NEUTRINO_E_NU_BAR' )
      call ComputeDegeneracyParameter_EQ_Kernel &
             ( Eta_EQ, Mu_E, Mu_NP, T, Sign = -1.0_KDR )
    end select !-- NM % Type

    end associate !-- Mu_E, etc.
    end associate !-- F

  end subroutine ComputeDegeneracyParameter_EQ


  subroutine ComputeDegeneracyParameter_EQ_Kernel &
               ( Eta_EQ, Mu_E, Mu_NP, T, Sign )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Eta_EQ
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Mu_E, &
      Mu_NP, &
      T
    real ( KDR ), intent ( in ) :: &
      Sign

    integer ( KDI ) :: &
      iV, &
      nValues

    nValues = size ( Eta_EQ )

    !$OMP parallel do private ( iV ) 
    do iV = 1, nValues
      Eta_EQ ( iV )  =  Sign  *  ( Mu_E ( iV )  -  Mu_NP ( iV ) )  /  T ( iV )
    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeDegeneracyParameter_EQ_Kernel


end module Interactions_NM_G_1__Form
