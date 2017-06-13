module Interactions_NM_1_G__Form

  !-- Interactions_NeutrinoMoments_1_Grey__Form
  
  use Basics
  use Mathematics
  use Fluid_P_MHN__Form
  use NeutrinoMoments_G__Form
  use Interactions_Template

  implicit none
  private

  type, public, extends ( InteractionsTemplate ) :: Interactions_NM_1_G_Form
    real ( KDR ) :: &
      RegulationParameter = 0.01_KDR
    class ( Fluid_P_MHN_Form ), pointer :: &
      Fluid => null ( )
    class ( NeutrinoMoments_G_Form ), pointer :: &
      Neutrinos_E_Nu          => null ( ), &
      Neutrinos_E_NuBar       => null ( ), &
      Neutrinos_MuTau_NuNuBar => null ( )
  contains
    procedure, private, pass :: &
      InitializeAllocate_NM_1_G
    generic, public :: &
      Initialize => InitializeAllocate_NM_1_G
    procedure, private, pass :: &
      Set_NM_1_G
    generic, public :: &
      Set => Set_NM_1_G
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      Regulate
    final :: &
      Finalize
    procedure, private, pass ( I ) :: &
      ComputeEquilibrium_T_Eta
    procedure, public, nopass :: &
      Compute_NuE_N_EA
    procedure, public, nopass :: &
      Compute_NuBarE_N_EA
    procedure, public, nopass :: &
      Compute_Nu_N_A_S
    ! procedure, public, pass :: &
    !   Compute_E_Nu_Nuclei_Kernel
    ! procedure, public, pass :: &
    !   Compute_E_NuBar_Kernel
  end type Interactions_NM_1_G_Form

    private :: &
      ComputeEquilibrium_T_Eta_Kernel, &
      RegulateKernel

    real ( KDR ), private :: &
      Pi             =  CONSTANT % PI, &
      TwoPi          =  2.0_KDR * CONSTANT % PI, &
      FourPi         =  4.0_KDR * CONSTANT % PI, &
      AMU            =  CONSTANT % ATOMIC_MASS_UNIT, &
      m_n            =  CONSTANT % NEUTRON_MASS, &
      m_p            =  CONSTANT % PROTON_MASS, &
      G_F            =  CONSTANT % FERMI_COUPLING, &
      Sin_2_Theta_W  =  CONSTANT % SIN_2_WEINBERG, &
      g_A            =  CONSTANT % NEUTRON_AXIAL_COUPLING

contains


  subroutine InitializeAllocate_NM_1_G &
               ( I, LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
                 VariableOption, NameOption, ClearOption, UnitOption )

    class ( Interactions_NM_1_G_Form ), intent ( inout ) :: &
      I
    type ( MeasuredValueForm ), intent ( in ) :: &
      LengthUnit, &
      EnergyDensityUnit, &
      TemperatureUnit
    integer ( KDI ), intent ( in ) :: &
      nValues
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      VariableOption
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ClearOption
    type ( MeasuredValueForm ), dimension ( : ), intent ( in ), optional :: &
      UnitOption

    if ( I % Type == '' ) &
      I % Type = 'an Interactions_NM_1_G'

    call I % InitializeTemplate &
           ( LengthUnit, EnergyDensityUnit, TemperatureUnit, nValues, &
             VariableOption, NameOption, ClearOption, UnitOption )

  end subroutine InitializeAllocate_NM_1_G


  subroutine Set_NM_1_G &
               ( I, Neutrinos_E_Nu, Neutrinos_E_NuBar, &
                 Neutrinos_MuTau_NuNuBar, Fluid )

    class ( Interactions_NM_1_G_Form ), intent ( inout ) :: &
      I
    class ( NeutrinoMoments_G_Form ), intent ( in ), target :: &
      Neutrinos_E_Nu, &
      Neutrinos_E_NuBar, &
      Neutrinos_MuTau_NuNuBar
    class ( Fluid_P_MHN_Form ), intent ( in ), target :: &
      Fluid

    I % Fluid                    =>  Fluid
    I % Neutrinos_E_Nu           =>  Neutrinos_E_Nu
    I % Neutrinos_E_NuBar        =>  Neutrinos_E_NuBar
    I % Neutrinos_MuTau_NuNuBar  =>  Neutrinos_MuTau_NuNuBar 

  end subroutine Set_NM_1_G


  subroutine Compute ( I, Current )

    class ( Interactions_NM_1_G_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current

    call Show ( 'Computing Interactions', I % IGNORABILITY )
    call Show ( Current % Name, 'Species', I % IGNORABILITY )

    call Clear ( I % Value ( :, I % EMISSIVITY_J ) )
    call Clear ( I % Value ( :, I % EMISSIVITY_H ) )
    call Clear ( I % Value ( :, I % EMISSIVITY_N ) )
    call Clear ( I % Value ( :, I % OPACITY_J ) )
    call Clear ( I % Value ( :, I % OPACITY_H ) )
    call Clear ( I % Value ( :, I % OPACITY_N ) )

    select type ( NM => Current )
    class is ( NeutrinoMoments_G_Form )

    associate ( F  =>  I % Fluid )

    select case ( trim ( Current % Type ) )
    case ( 'NEUTRINOS_E_NU' )
       
      call I % Compute_NuE_N_EA &
             ( I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % EMISSIVITY_N ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               I % Value ( :, I % OPACITY_N ), &
               F % Value ( :, F % BARYON_MASS ), &
               F % Value ( :, F % COMOVING_DENSITY ), &
               F % Value ( :, F % TEMPERATURE ), &
               F % Value ( :, F % MASS_FRACTION_PROTON ), &
               F % Value ( :, F % MASS_FRACTION_NEUTRON ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
               NM % Value ( :, NM % TEMPERATURE_PARAMETER ), &
               NM % Value ( :, NM % DEGENERACY_PARAMETER ), &
               F % Value ( :, F % ELECTRON_FRACTION ) )

      call I % Compute_Nu_N_A_S &
             ( I % Value ( :, I % OPACITY_H ), &
               F % Value ( :, F % BARYON_MASS ), &
               F % Value ( :, F % COMOVING_DENSITY ), &
               F % Value ( :, F % MASS_FRACTION_PROTON ), &
               F % Value ( :, F % MASS_FRACTION_NEUTRON ), &
               F % Value ( :, F % MASS_FRACTION_HEAVY ), &
               F % Value ( :, F % HEAVY_ATOMIC_NUMBER ), &
               F % Value ( :, F % HEAVY_MASS_NUMBER ), &
               NM % Value ( :, NM % TEMPERATURE_PARAMETER ), &
               NM % Value ( :, NM % DEGENERACY_PARAMETER ) )

    case ( 'NEUTRINOS_E_NU_BAR' )

      call I % Compute_NuBarE_N_EA &
             ( I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % EMISSIVITY_N ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               I % Value ( :, I % OPACITY_N ), &
               F % Value ( :, F % BARYON_MASS ), &
               F % Value ( :, F % COMOVING_DENSITY ), &
               F % Value ( :, F % TEMPERATURE ), &
               F % Value ( :, F % MASS_FRACTION_PROTON ), &
               F % Value ( :, F % MASS_FRACTION_NEUTRON ), &
               F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
               NM % Value ( :, NM % TEMPERATURE_PARAMETER ), &
               NM % Value ( :, NM % DEGENERACY_PARAMETER ), &
               F % Value ( :, F % ELECTRON_FRACTION ) )

      call I % Compute_Nu_N_A_S &
             ( I % Value ( :, I % OPACITY_H ), &
               F % Value ( :, F % BARYON_MASS ), &
               F % Value ( :, F % COMOVING_DENSITY ), &
               F % Value ( :, F % MASS_FRACTION_PROTON ), &
               F % Value ( :, F % MASS_FRACTION_NEUTRON ), &
               F % Value ( :, F % MASS_FRACTION_HEAVY ), &
               F % Value ( :, F % HEAVY_ATOMIC_NUMBER ), &
               F % Value ( :, F % HEAVY_MASS_NUMBER ), &
               NM % Value ( :, NM % TEMPERATURE_PARAMETER ), &
               NM % Value ( :, NM % DEGENERACY_PARAMETER ) )

!    case default
!      call Show ( 'Radiation Type not recognized', CONSOLE % ERROR )
!      call Show ( 'Interactions_NM_G_1__Form', 'module', CONSOLE % ERROR )
!      call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
    end select !-- Radiation % Type

    end associate !-- F
    end select !-- NM
    
  end subroutine Compute


  subroutine Regulate ( I, Current, dT )

    class ( Interactions_NM_1_G_Form ), intent ( inout ) :: &
      I
    class ( CurrentTemplate ), intent ( in ) :: &
      Current
    real ( KDR ), intent ( in ) :: &
      dT

    call Show ( 'Regulating Interactions', I % IGNORABILITY )

    select type ( NM => Current )
    class is ( NeutrinoMoments_G_Form )

    associate ( F => I % Fluid )

    call RegulateKernel &
           ( I % Value ( :, I % EMISSIVITY_J ), &
             I % Value ( :, I % EMISSIVITY_H ), &
             I % Value ( :, I % EMISSIVITY_N ), &
             I % Value ( :, I % OPACITY_J ), &
             I % Value ( :, I % OPACITY_H ), &
             I % Value ( :, I % OPACITY_N ), &
             NM % Value ( :, NM % COMOVING_ENERGY ), &
             F % Value ( :, F % INTERNAL_ENERGY ), &
             I % RegulationParameter, dT )

    end associate !-- F
    end select !-- NM

  end subroutine Regulate


  impure elemental subroutine Finalize ( I )

    type ( Interactions_NM_1_G_Form ), intent ( inout ) :: &
      I

    nullify ( I % Neutrinos_MuTau_NuNuBar )
    nullify ( I % Neutrinos_E_NuBar )
    nullify ( I % Neutrinos_E_Nu )
    nullify ( I % Fluid )

    call I % FinalizeTemplate ( )

  end subroutine Finalize


  subroutine ComputeEquilibrium_T_Eta ( T_EQ, Eta_EQ, I, C )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ, &
      Eta_EQ
    class ( Interactions_NM_1_G_Form ), intent ( in ) :: &
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
      call ComputeEquilibrium_T_Eta_Kernel &
             ( T_EQ, Eta_EQ, Mu_E, Mu_NP, T, Sign = 1.0_KDR )
    case ( 'NEUTRINOS_E_NU_BAR' )
      call ComputeEquilibrium_T_Eta_Kernel &
             ( T_EQ, Eta_EQ, Mu_E, Mu_NP, T, Sign = -1.0_KDR )
    end select !-- NM % Type

    end associate !-- Mu_E, etc.
    end associate !-- F

  end subroutine ComputeEquilibrium_T_Eta


  subroutine Compute_NuE_N_EA &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, M, N, T, X_p, X_n, Mu_e, &
                 T_nu, Eta_nu, Y_e )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Xi_J, Xi_N, &
      Chi_J, Chi_H, Chi_N
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T, &
      X_p, X_n, &
      Mu_e, &
      T_nu, &
      Eta_nu, &
      Y_e

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Factor_Chi, Factor_Xi, &
      Q, &
      N_p, N_n, &
      Eta_e, Eta_e_Q, &
      Fermi_2_e_Q, Fermi_3_e_Q, Fermi_4_e_Q, Fermi_5_e_Q, &
      Fermi_2_nu, Fermi_3_nu, Fermi_4_nu, Fermi_5_nu, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta, &
      OneMinus_F_e_J, OneMinus_F_nu_J, &
      OneMinus_F_e_N, OneMinus_F_nu_N

    nValues  =  size ( Xi_J )
    
    Factor_Chi  =  G_F ** 2  /  Pi  *  ( 1  +  3 * g_A ** 2 )
     Factor_Xi  =  FourPi / TwoPi ** 3  *  Factor_Chi
             Q  =  m_n - m_p

    !$OMP parallel do &
    !$OMP   private ( iV, N_p, N_n, Eta_e, Eta_e_Q, &
    !$OMP             Fermi_2_e_Q, Fermi_3_e_Q, Fermi_4_e_Q, Fermi_5_e_Q, &
    !$OMP             Fermi_2_nu, Fermi_3_nu, Fermi_4_nu, Fermi_5_nu, &
    !$OMP             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta, &
    !$OMP             OneMinus_F_e_J, OneMinus_F_nu_J, & 
    !$OMP             OneMinus_F_e_N, OneMinus_F_nu_N ) 
    do iV = 1, nValues

      if ( T ( iV ) == 0.0_KDR ) &
        cycle

      N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  AMU
      N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  AMU

      Eta_e    =  Mu_e ( iV )  /  T ( iV )
      Eta_e_Q  =  ( Mu_e ( iV )  -  Q )  /  T ( iV )

      call DFERMI ( 2.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_2_e_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_3_e_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 4.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_4_e_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 5.0_KDR, Eta_e_Q, 0.0_KDR, Fermi_5_e_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      call DFERMI ( 2.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_2_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 4.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_4_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 5.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_5_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      OneMinus_F_e_J   =  1.0_KDR
      OneMinus_F_nu_J  =  1.0_KDR
      if ( Fermi_4_nu > 0.0_KDR ) &
        OneMinus_F_e_J  &
          =  1. / ( exp ( Eta_e  -  T_nu ( iV ) / T ( iV ) &
                                    * Fermi_5_nu / Fermi_4_nu )  &
                  +  1. )
      if ( Fermi_4_e_Q > 0.0_KDR ) &
        OneMinus_F_nu_J  &
          =  1. / ( exp ( Eta_nu ( iV )  -  T ( iV ) / T_nu ( iV ) &
                                            * Fermi_5_e_Q / Fermi_4_e_Q )  &
                    +  1. )

      OneMinus_F_e_N   =  1.0_KDR
      OneMinus_F_nu_N  =  1.0_KDR
      if ( Fermi_3_nu > 0.0_KDR ) &
        OneMinus_F_e_N  &
          =  1. / ( exp ( Eta_e  -  T_nu ( iV ) / T ( iV ) &
                                    * Fermi_4_nu / Fermi_3_nu )  &
                  +  1. )
      if ( Fermi_3_e_Q > 0.0_KDR ) &
        OneMinus_F_nu_N  &
          =  1. / ( exp ( Eta_nu ( iV )  -  T ( iV ) / T_nu ( iV ) &
                                            * Fermi_4_e_Q / Fermi_3_e_Q )  &
                    +  1. )

!if ( Y_e ( iV ) > 0.1_KDR ) then

      Xi_J ( iV )  &
        =  Xi_J ( iV )  &
           +  Factor_Xi  *  N_p  *  T ( iV ) ** 4 &
              *  (    T ( iV ) ** 2     *  Fermi_5_e_Q  &
                   +  2 * Q * T ( iV )  *  Fermi_4_e_Q  &
                   +  Q ** 2            *  Fermi_3_e_Q )  &
              *  OneMinus_F_nu_J

      Xi_N ( iV )  &
        =  Xi_N ( iV )  &
           +  Factor_Xi  *  N_p  *  T ( iV ) ** 3 &
              *  (    T ( iV ) ** 2     *  Fermi_4_e_Q  &
                   +  2 * Q * T ( iV )  *  Fermi_3_e_Q  &
                   +  Q ** 2            *  Fermi_2_e_Q )  &
              *  OneMinus_F_nu_N

! else
!   if ( iV > 2 ) then 
! call Show ( '>>> Preventing Y_e < 0.2', CONSOLE % ERROR )
! call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
! call Show ( iV, '>>> iV', CONSOLE % ERROR )
! call Show ( 'Compute_NuE_N', '>>> Subroutine', CONSOLE % ERROR )
!   end if
! end if

!if ( Y_e ( iV ) < 0.51_KDR ) then

      if ( Fermi_2_nu > 0.0_KDR .and. Fermi_3_nu > 0.0_KDR ) then

        Chi_J ( iV )  &
          =  Chi_J ( iV )  &
             +  Factor_Chi  *  N_n  /  Fermi_3_nu &
                *  (    T_nu ( iV ) ** 2     *  Fermi_5_nu  &
                     +  2 * Q * T_nu ( iV )  *  Fermi_4_nu  &
                     +  Q ** 2               *  Fermi_3_nu )  &
                *  OneMinus_F_e_J

        Chi_H ( iV )  &
          =  Chi_H ( iV )  &
             +  Factor_Chi  *  N_n  /  Fermi_3_nu &
                *  (    T_nu ( iV ) ** 2     *  Fermi_5_nu  &
                     +  2 * Q * T_nu ( iV )  *  Fermi_4_nu  &
                     +  Q ** 2               *  Fermi_3_nu )  &
                *  OneMinus_F_e_J

      end if
     
      if ( Fermi_2_nu > 0.0_KDR ) then
        Chi_N ( iV )  &
          =  Chi_N ( iV )  &
             +  Factor_Chi  *  N_n  /  Fermi_2_nu &
                *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu  &
                     +  2 * Q * T_nu ( iV )  *  Fermi_3_nu  &
                     +  Q ** 2               *  Fermi_2_nu )  &
                *  OneMinus_F_e_N
      end if

! else 
! call Show ( '>>> Preventing Y_e > 0.51', CONSOLE % ERROR )
! call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
! call Show ( iV, '>>> iV', CONSOLE % ERROR )
! call Show ( 'Compute_NuE_N', '>>> Subroutine', CONSOLE % ERROR )
! end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_NuE_N_EA


  subroutine Compute_NuBarE_N_EA &
               ( Xi_J, Xi_N, Chi_J, Chi_H, Chi_N, M, N, T, X_p, X_n, Mu_e, &
                 T_nu, Eta_nu, Y_e )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Xi_J, Xi_N, &
      Chi_J, Chi_H, Chi_N
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      T, &
      X_p, X_n, &
      Mu_e, &
      T_nu, &
      Eta_nu, &
      Y_e

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Factor_Chi, Factor_Xi, &
      Q, &
      N_p, N_n, &
      Eta_e, &
      Eta_nu_Q, &
      Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e, &
      Fermi_2_nu_Q, Fermi_3_nu_Q, Fermi_4_nu_Q, Fermi_5_nu_Q, &
      Fermi_2_nu, Fermi_3_nu, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta, &
      OneMinus_F_e_J, OneMinus_F_nu_J, &
      OneMinus_F_e_N, OneMinus_F_nu_N

    nValues  =  size ( Xi_J )
    
    Factor_Chi  =  G_F ** 2  /  Pi  *  ( 1  +  3 * g_A ** 2 )
     Factor_Xi  =  FourPi / TwoPi ** 3  *  Factor_Chi
             Q  =  m_n - m_p

    !$OMP parallel do &
    !$OMP   private ( iV, N_p, N_n, Eta_e, Eta_nu_Q, &
    !$OMP             Fermi_2_e, Fermi_3_e, Fermi_4_e, Fermi_5_e, &
    !$OMP             Fermi_2_nu_Q, Fermi_3_nu_Q, Fermi_4_nu_Q, Fermi_5_nu_Q, &
    !$OMP             Fermi_2_nu, Fermi_3_nu, &
    !$OMP             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta, & 
    !$OMP             OneMinus_F_e_J, OneMinus_F_nu_J, & 
    !$OMP             OneMinus_F_e_N, OneMinus_F_nu_N ) 
    do iV = 1, nValues

      if ( T ( iV ) == 0.0_KDR ) &
        cycle

      N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  AMU
      N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  AMU

      Eta_e  =  Mu_e ( iV )  /  T ( iV )
      call DFERMI ( 2.0_KDR, -Eta_e, 0.0_KDR, Fermi_2_e, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, -Eta_e, 0.0_KDR, Fermi_3_e, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 4.0_KDR, -Eta_e, 0.0_KDR, Fermi_4_e, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 5.0_KDR, -Eta_e, 0.0_KDR, Fermi_5_e, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      Eta_nu_Q  =  Eta_nu ( iV )  -  Q / T ( iV )
      call DFERMI ( 2.0_KDR, Eta_nu_Q, 0.0_KDR, Fermi_2_nu_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta_nu_Q, 0.0_KDR, Fermi_3_nu_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 4.0_KDR, Eta_nu_Q, 0.0_KDR, Fermi_4_nu_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 5.0_KDR, Eta_nu_Q, 0.0_KDR, Fermi_5_nu_Q, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      call DFERMI ( 2.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_2_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      OneMinus_F_e_J   =  1.0_KDR
      OneMinus_F_nu_J  =  1.0_KDR
      if ( Fermi_4_nu_Q > 0.0_KDR ) &
        OneMinus_F_e_J  &
        =  1. / ( exp ( -Eta_e  -  T_nu ( iV ) / T ( iV ) &
                                  * Fermi_5_nu_Q / Fermi_4_nu_Q )  &
                  +  1. )
      if ( Fermi_4_e > 0.0_KDR ) &
        OneMinus_F_nu_J  &
          =  1. / ( exp ( Eta_nu ( iV )  -  T ( iV ) / T_nu ( iV ) &
                                            * Fermi_5_e / Fermi_4_e )  &
                    +  1. )

      OneMinus_F_e_N   =  1.0_KDR
      OneMinus_F_nu_N  =  1.0_KDR
      if ( Fermi_3_nu_Q > 0.0_KDR ) &
        OneMinus_F_e_N  &
        =  1. / ( exp ( -Eta_e  -  T_nu ( iV ) / T ( iV ) &
                                  * Fermi_4_nu_Q / Fermi_3_nu_Q )  &
                  +  1. )
      if ( Fermi_3_e > 0.0_KDR ) &
        OneMinus_F_nu_N  &
          =  1. / ( exp ( Eta_nu ( iV )  -  T ( iV ) / T_nu ( iV ) &
                                            * Fermi_4_e / Fermi_3_e )  &
                    +  1. )

!if ( Y_e ( iV ) < 0.51_KDR ) then

      Xi_J ( iV )  &
        =  Xi_J ( iV )  &
           +  Factor_Xi  *  N_n  *  T ( iV ) ** 3 &
              *  (    T ( iV ) ** 3              *  Fermi_5_e  &
                   +  3  *  Q  *  T ( iV ) ** 2  *  Fermi_4_e  &
                   +  3  *  Q ** 2  *  T ( iV )  *  Fermi_3_e  &
                   +  Q ** 3                     *  Fermi_2_e )  &
              *  OneMinus_F_nu_J

      Xi_N ( iV )  &
        =  Xi_N ( iV )  &
           +  Factor_Xi  *  N_n  *  T ( iV ) ** 3 &
              *  (    T ( iV ) ** 2     *  Fermi_4_e  &
                   +  2 * Q * T ( iV )  *  Fermi_3_e  &
                   +  Q ** 2            *  Fermi_2_e )  &
              *  OneMinus_F_nu_N

! else 
! call Show ( '>>> Preventing Y_e > 0.51', CONSOLE % ERROR )
! call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
! call Show ( iV, '>>> iV', CONSOLE % ERROR )
! call Show ( 'Compute_NuBarE_N', '>>> Subroutine', CONSOLE % ERROR )
! end if

!if ( Y_e ( iV ) > 0.1_KDR ) then

      if ( T_nu ( iV ) * Fermi_3_nu > 0.0_KDR ) then

        Chi_J ( iV )  &
          =  Chi_J ( iV )  &
             +  Factor_Chi  *  N_p  /  ( T_nu ( iV )  *  Fermi_3_nu )  &
                *  (    T_nu ( iV ) ** 3              *  Fermi_5_nu_Q  &
                     +  3  *  Q  *  T_nu ( iV ) ** 2  *  Fermi_4_nu_Q  &
                     +  3  *  Q ** 2  *  T_nu ( iV )  *  Fermi_3_nu_Q  &
                     +  Q ** 3                        *  Fermi_2_nu_Q )  &
                *  OneMinus_F_e_J

        Chi_H ( iV )  &
          =  Chi_H ( iV )  &
             +  Factor_Chi  *  N_p  /  ( T_nu ( iV )  *  Fermi_3_nu )  &
                *  (    T_nu ( iV ) ** 3              *  Fermi_5_nu_Q  &
                     +  3  *  Q  *  T_nu ( iV ) ** 2  *  Fermi_4_nu_Q  &
                     +  3  *  Q ** 2  *  T_nu ( iV )  *  Fermi_3_nu_Q  &
                     +  Q ** 3                        *  Fermi_2_nu_Q )  &
                *  OneMinus_F_e_J

      end if
     
      if ( T_nu ( iV ) * Fermi_2_nu > 0.0_KDR ) then
        Chi_N ( iV )  &
          =  Chi_N ( iV )  &
             +  Factor_Chi  *  N_p  /  ( T_nu ( iV )  *  Fermi_2_nu )  &
                *  (    T_nu ( iV ) ** 2     *  Fermi_4_nu_Q  &
                     +  2 * Q * T_nu ( iV )  *  Fermi_3_nu_Q  &
                     +  Q ** 2               *  Fermi_2_nu_Q )  &
                *  OneMinus_F_e_N
      end if

! else 
!   if ( iV > 2 ) then
! call Show ( '>>> Preventing Y_e < 0.2', CONSOLE % ERROR )
! call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
! call Show ( iV, '>>> iV', CONSOLE % ERROR )
! call Show ( 'Compute_NuBarE_N', '>>> Subroutine', CONSOLE % ERROR )
!   end if
! end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_NuBarE_N_EA


  subroutine Compute_Nu_N_A_S ( Chi_H, M, N, X_p, X_n, X_A, Z, A, T_nu, Eta_nu )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Chi_H
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      X_p, X_n, X_A, &
      Z, A, &
      T_nu, &
      Eta_nu

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Factor_n, Factor_p, Factor_A, &
      N_p, N_n, N_A, &
      Fermi_3_nu, Fermi_5_nu, &
      fdeta, fdeta2, &
      fdtheta, fdtheta2, &
      fdetadtheta

    nValues  =  size ( Chi_H )
    
    Factor_p  =  2  *  G_F ** 2  /  3 * Pi  &
                 *  ( ( 2 * Sin_2_Theta_W  -  1. / 2. ) ** 2 &
                      +  5. / 4. * g_A ** 2 )

    Factor_n  =  2  *  G_F ** 2  /  3 * Pi  &
                 *  ( 1. / 4.  +  5. / 4. * g_A ** 2 )

    !$OMP parallel do &
    !$OMP   private ( iV, Factor_A, N_p, N_n, N_A, Fermi_3_nu, Fermi_5_nu, &
    !$OMP             fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta ) 
    do iV = 1, nValues

      if ( A ( iV ) == 0.0_KDR ) &
        cycle

      Factor_A  =  2  *  G_F ** 2  /  3 * Pi  &
                   *  ( Z ( iV ) / A ( iV )  *  ( 1  -  2 * Sin_2_Theta_W ) &
                        -  1. / 2. ) ** 2  &
                   *  A ( iV ) ** 2

      N_p  =  M ( iV )  *  N ( iV )  *  X_p ( iV )  /  AMU
      N_n  =  M ( iV )  *  N ( iV )  *  X_n ( iV )  /  AMU
      N_A  =  M ( iV )  *  N ( iV )  *  X_A ( iV )  /  ( A ( iV ) * AMU )

      call DFERMI ( 3.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_3_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )
      call DFERMI ( 5.0_KDR, Eta_nu ( iV ), 0.0_KDR, Fermi_5_nu, &
                    fdeta, fdtheta, fdeta2, fdtheta2, fdetadtheta )

      if ( Fermi_3_nu > 0.0_KDR ) then

        Chi_H ( iV )  &
          =  Chi_H ( iV )  &
             +  ( Factor_p * N_p  +  Factor_n * N_n  +  Factor_A * N_A )  &
                *  T_nu ( iV ) ** 2  *  Fermi_5_nu / Fermi_3_nu 

      end if
     
    end do !-- iV
    !$OMP end parallel do

  end subroutine Compute_Nu_N_A_S


  subroutine ComputeEquilibrium_T_Eta_Kernel &
               ( T_EQ, Eta_EQ, Mu_E, Mu_NP, T, Sign )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      T_EQ, &
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

      T_EQ ( iV )  =  T ( iV )

      if ( T ( iV ) > 0.0_KDR ) then
        Eta_EQ ( iV )  =  Sign  *  ( Mu_E ( iV )  -  Mu_NP ( iV ) )  &
                          /  T ( iV )
      else
        Eta_EQ ( iV )  =  0.0_KDR
      end if

    end do !-- iV
    !$OMP end parallel do

  end subroutine ComputeEquilibrium_T_Eta_Kernel


  subroutine RegulateKernel &
               ( Xi_J, Xi_H, Xi_N, Chi_J, Chi_H, Chi_N, J, E, R, dT )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      Xi_J, Xi_H, Xi_N, &
      Chi_J, Chi_H, Chi_N
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      E
    real ( KDR ), intent ( in ) :: &
      R, &
      dT

    integer ( KDI ) :: &
      iV, &
      nValues
    real ( KDR ) :: &
      Ratio  

    nValues  =  size ( Xi_J )
    
    !$OMP parallel do private ( iV, Ratio )
    do iV = 1, nValues

      if ( E ( iV ) == 0.0_KDR ) &
        cycle

      Ratio  =  abs ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) * dT  &
                /  E ( iV )

      if ( Ratio > R ) then
         Xi_J ( iV )  =  R / Ratio  *   Xi_J ( iV )
         Xi_H ( iV )  =  R / Ratio  *   Xi_H ( iV )
         Xi_N ( iV )  =  R / Ratio  *   Xi_N ( iV )
        Chi_J ( iV )  =  R / Ratio  *  Chi_J ( iV )
        Chi_H ( iV )  =  R / Ratio  *  Chi_H ( iV )
        Chi_N ( iV )  =  R / Ratio  *  Chi_N ( iV )
call Show ( '>>> Regulating Interactions', CONSOLE % ERROR )
call Show ( PROGRAM_HEADER % Communicator % Rank, '>>> Rank', CONSOLE % ERROR )
call Show ( iV, '>>> iV', CONSOLE % ERROR )
call Show ( R / Ratio, '>>> RegulationFactor', CONSOLE % ERROR )
      end if
    end do !-- iV
    !$OMP end parallel do

  end subroutine RegulateKernel


end module Interactions_NM_1_G__Form
