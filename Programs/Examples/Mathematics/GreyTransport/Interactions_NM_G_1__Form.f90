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
    class ( Fluid_P_MHN_Form ), pointer :: &
      Fluid => null ( )
    class ( NeutrinoMomentsForm ), pointer :: &
      Neutrinos_E_Nu          => null ( ), &
      Neutrinos_E_NuBar       => null ( ), &
      Neutrinos_MuTau_NuNuBar => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_NM_G_1
    generic, public :: &
      Initialize => InitializeAllocate_NM_G_1
    procedure, private, pass :: &
      Set_NM_G_1
    generic, public :: &
      Set => Set_NM_G_1
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
    procedure, public, pass ( I ) :: &
      ComputeDegeneracyParameter_EQ
  end type Interactions_NM_G_1_Form

    private :: &
      ComputeDegeneracyParameter_EQ_Kernel

contains


  subroutine InitializeAllocate_NM_G_1 &
               ( I, LengthUnit, EnergyDensityUnit, nValues, NameOption, &
                 ClearOption, UnitOption )

    class ( Interactions_NM_G_1_Form ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_NM_G_1'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, nValues, NameOption, &
             ClearOption, UnitOption )

  end subroutine InitializeAllocate_NM_G_1


  subroutine Set_NM_G_1 &
               ( I, Neutrinos_E_Nu, Neutrinos_E_NuBar, &
                 Neutrinos_MuTau_NuNuBar, Fluid )

    class ( Interactions_NM_G_1_Form ), intent ( inout ) :: &
      I
    class ( NeutrinoMomentsForm ), intent ( in ), target :: &
      Neutrinos_E_Nu, &
      Neutrinos_E_NuBar, &
      Neutrinos_MuTau_NuNuBar
    class ( Fluid_P_MHN_Form ), intent ( in ), target :: &
      Fluid

    I % Fluid                    =>  Fluid
    I % Neutrinos_E_Nu           =>  Neutrinos_E_Nu
    I % Neutrinos_E_NuBar        =>  Neutrinos_E_NuBar
    I % Neutrinos_MuTau_NuNuBar  =>  Neutrinos_MuTau_NuNuBar 

  end subroutine Set_NM_G_1


  subroutine Compute ( I )

    class ( Interactions_NM_G_1_Form ), intent ( inout ) :: &
      I

    !-- FIXME: fill in

  end subroutine Compute


  impure elemental subroutine Finalize ( I )

    type ( Interactions_NM_G_1_Form ), intent ( inout ) :: &
      I

    nullify ( I % Neutrinos_MuTau_NuNuBar )
    nullify ( I % Neutrinos_E_NuBar )
    nullify ( I % Neutrinos_E_Nu )
    nullify ( I % Fluid )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


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
    case ( 'NEUTRINOS_E_NU' )
      call ComputeDegeneracyParameter_EQ_Kernel &
             ( Eta_EQ, Mu_E, Mu_NP, T, Sign = 1.0_KDR )
    case ( 'NEUTRINOS_E_NU_BAR' )
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
