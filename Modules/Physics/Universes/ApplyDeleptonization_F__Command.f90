module ApplyDeleptonization_F__Command

  !-- M. Liebendorfer, ApJ 633 1024 (2005)

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  public :: &
    ApplyDeleptonization_F

    private :: &
      Apply_DP_DS_Kernel

    real ( KDR ), private :: &
      Y_1         = 0.5_KDR, &     !-- N13
      Y_2         = 0.285_KDR, &   !-- N13
      Y_C         = 0.035_KDR, &   !-- N13
      Rho_1_cgs   = 2.0e7_KDR, &   !-- N13
      Rho_2_cgs   = 2.0e13_KDR, &  !-- N13
      ! Y_1         = 0.5_KDR, &     !-- G15
      ! Y_2         = 0.278_KDR, &   !-- G15
      ! Y_C         = 0.035_KDR      !-- G15
      ! Rho_1       = 3.0e7_KDR, &   !-- G15
      ! Rho_2       = 2.0e13_KDR, &  !-- G15
      RhoTrap_cgs = 2.0e12_KDR, & 
      E_Esc_MeV   = 10.0_KDR, &    
      N_1, N_2, N_Trap, &
      E_Esc

contains


  subroutine ApplyDeleptonization_F &
               ( S, Sources_F, Increment, Fluid, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_F
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iProton, &
      iEntropy
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    N_1     =  Rho_1_cgs    &
                 *  UNIT % MASS_DENSITY_CGS / UNIT % ATOMIC_MASS_UNIT 
    N_2     =  Rho_2_cgs    &
                 *  UNIT % MASS_DENSITY_CGS / UNIT % ATOMIC_MASS_UNIT 
    N_Trap  =  RhoTrap_cgs  &
                 *  UNIT % MASS_DENSITY_CGS / UNIT % ATOMIC_MASS_UNIT
    E_Esc   =  E_Esc_MeV &
                 *  UNIT % MEGA_ELECTRON_VOLT

    select type ( F => Fluid )
    class is ( Fluid_P_HN_Form )

    select type ( S_F => Sources_F )
    class is ( Sources_F_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    call Search ( F % iaConserved, F % CONSERVED_PROTON_DENSITY, iProton )
    call Search ( F % iaConserved, F % CONSERVED_ENTROPY, iEntropy )

    if ( iStage == 1 ) then
      call Clear ( S_F % Value ( :, S_F % RADIATION_DP ) )
      call Clear ( S_F % Value ( :, S_F % RADIATION_DS ) )
    end if

    call Apply_DP_DS_Kernel &
           ( Increment % Value ( :, iProton ), &
             Increment % Value ( :, iEntropy ), &
             S_F % Value ( :, S_F % RADIATION_DP ), &
             S_F % Value ( :, S_F % RADIATION_DS ), &
             Chart % IsProperCell, &
             F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             F % Value ( :, F % TEMPERATURE ), &
             F % Value ( :, F % PROTON_FRACTION ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
             F % Value ( :, F % CHEMICAL_POTENTIAL_E ), &
             F % Value ( :, F % CONSERVED_ENTROPY ), &
             TimeStep, S % B ( iStage ) )

    end select !-- Chart
    end select !-- S_F
    end select !-- F

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyDeleptonization_F


  subroutine Apply_DP_DS_Kernel &
               ( K_DP, K_DS, S_DP, S_DS, IsProperCell, N, T, YP, &
                 Mu_NP, Mu_E, DS, dt, Weight_RK )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_DP, &
      K_DS, &
      S_DP, &
      S_DS
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      T, &
      YP, &
      Mu_NP, &
      Mu_E, &
      DS
    real ( KDR ) :: &
      dt, &
      Weight_RK

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      X, &
      abs_X, &
      YP_bar, &
      dYP, &
      dDP, &
      dDS

    nValues = size ( K_DP )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      if ( .not. IsProperCell ( iV ) ) &
        cycle

      X  =  max ( -1.0_KDR, &
                  min ( 1.0_KDR, &
                        ( 2 * log10 ( N ( iV ) ) &
                          - log10 ( N_2 ) - log10 ( N_1 ) ) &
                        / ( log10 ( N_2 ) - log10 ( N_1 ) ) ) )
      abs_X  =  abs ( X )

      YP_bar  =  0.5_KDR * ( Y_1 + Y_2 )  +  0.5_KDR * X * ( Y_2 - Y_1 ) &
                  +  Y_C * ( 1.0_KDR - abs_X  +  4.0_KDR * abs_X &
                                                 * ( abs_X - 0.5_KDR ) &
                                                 * ( abs_X - 1.0_KDR ) )

      dYP  =  min ( 0.0_KDR, YP_bar - YP ( iV ) )
      dDP  =  N ( iV ) * dYP

      if ( Mu_E ( iV ) - Mu_NP ( iV ) - E_Esc  >  0.0_KDR &
           .and. N ( iV ) < N_Trap ) &
      then
        dDS  =  - dDP * ( Mu_E ( iV ) - Mu_NP ( iV ) - E_Esc )  /  T ( iV )
      else
        dDS  =  0.0_KDR
      end if

      K_DP ( iV )  =  K_DP ( iV )  +  dDP
      S_DP ( iV )  =  S_DP ( iV )  +  Weight_RK * dDP / dt

      K_DS ( iV )  =  K_DS ( iV )  +  dDS
      S_DS ( iV )  =  S_DS ( iV )  +  Weight_RK * dDS / dt

    end do
    !$OMP end parallel do

  end subroutine Apply_DP_DS_Kernel


end module ApplyDeleptonization_F__Command
