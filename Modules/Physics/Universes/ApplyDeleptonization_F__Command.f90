module ApplyDeleptonization_F__Command

  !-- M. Liebendorfer, ApJ 633 1024 (2005)
  !-- But here neutrino stress only implemented in trapped region.

  !-- Actually there is a glitch in velocity caused by the discontinuity
  !   resulting from using the neutrino stress only from the trapped 
  !   region. Therefore commenting out neutrino stress.  
  
  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  public :: &
    ApplyDeleptonization_F

    private :: &
      Apply_DP_DS_Kernel!, &
!      Apply_S_D_Kernel, &
!      ComputeNeutrinoPressure

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
    type ( StorageForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      Fluid
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iD, &
      iProton, &
      iEntropy, &
      iMomentum
    type ( StorageForm ) :: &
      NeutrinoPressure
    type ( TimerForm ), pointer :: &
      Timer
    type ( GradientForm ) :: &
      Gradient

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

    ! call NeutrinoPressure % Initialize &
    !        ( [ F % nValues, 1 ] )
    ! call Gradient % Initialize &
    !        ( 'NeutrinoPressureGradient', [ F % nValues, 1 ] )
    ! call ComputeNeutrinoPressure &
    !        ( NeutrinoPressure % Value ( :, 1 ), &
    !          Chart % IsProperCell, &
    !          F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
    !          F % Value ( :, F % TEMPERATURE ), &
    !          F % Value ( :, F % CHEMICAL_POTENTIAL_N_P ), &
    !          F % Value ( :, F % CHEMICAL_POTENTIAL_E ) )
    ! call S % Current_ASC % Atlas_SC % ApplyBoundaryConditionsFaces &
    !        ( NeutrinoPressure )
    ! do iD = 1, Chart % nDimensions
    !   call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
    !                 iMomentum )
    !   if ( iStage == 1 ) &
    !     call Clear ( S_F % Value ( :, S_F % RADIATION_S_D ( iD ) ) )
    !   call Gradient % Compute ( Chart, NeutrinoPressure, iDimension = iD )
    !   call Apply_S_D_Kernel &
    !          ( Increment % Value ( :, iMomentum ), &
    !            S_F % Value ( :, S_F % RADIATION_S_D ( iD ) ), &
    !            Chart % IsProperCell, &
    !            Gradient % Output % Value ( :, 1 ), &
    !            TimeStep, S % B ( iStage ) )
    ! end do !-- iD

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


  ! subroutine Apply_S_D_Kernel &
  !              ( K_S_D, S_S_D, IsProperCell, dP_nu_dx, dt, Weight_RK )

  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     K_S_D, &
  !     S_S_D
  !   logical ( KDL ), dimension ( : ), intent ( in ) :: &
  !     IsProperCell
  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     dP_nu_dx
  !   real ( KDR ) :: &
  !     dt, &
  !     Weight_RK

  !   integer ( KDI ) :: &
  !     iV, &  !-- iValue
  !     nValues

  !   nValues = size ( K_S_D )

  !   !$OMP parallel do private ( iV )
  !   do iV = 1, nValues

  !     if ( .not. IsProperCell ( iV ) ) &
  !       cycle

  !     K_S_D ( iV )  =  K_S_D ( iV )  -  dP_Nu_dx ( iV ) * dt
  !     S_S_D ( iV )  =  S_S_D ( iV )  -  Weight_RK * dP_Nu_dx ( iV )

  !   end do
  !   !$OMP end parallel do

  ! end subroutine Apply_S_D_Kernel


  ! subroutine ComputeNeutrinoPressure ( P_Nu, IsProperCell, N, T, Mu_NP, Mu_E )

  !   real ( KDR ), dimension ( : ), intent ( inout ) :: &
  !     P_Nu
  !   logical ( KDL ), dimension ( : ), intent ( in ) :: &
  !     IsProperCell
  !   real ( KDR ), dimension ( : ), intent ( in ) :: &
  !     N, &
  !     T, &
  !     Mu_NP, &
  !     Mu_E

  !   integer ( KDI ) :: &
  !     iV, &  !-- iValue
  !     nValues
  !   real ( KDR ) :: &
  !     k, &
  !     a, b, c, d, &
  !     Eta, &
  !     F_3

  !   nValues = size ( P_Nu )

  !   k  =  CONSTANT % BOLTZMANN
  !   a  =  2.0_KDR  *  CONSTANT % PI ** 2
  !   b  =  7.0_KDR / 15.0_KDR   *  CONSTANT % PI ** 4
  !   c  =  7.0_KDR / 120.0_KDR  *  CONSTANT % PI ** 4
  !   d  =  4.0_KDR / 3.0_KDR    *  CONSTANT % PI  &
  !         /  ( 2.0_KDR  *  CONSTANT % PI  *  CONSTANT % PLANCK_REDUCED  &
  !              *  CONSTANT % SPEED_OF_LIGHT ) ** 3

  !   !$OMP parallel do private ( iV )
  !   do iV = 1, nValues

  !     P_Nu ( iV )  =  0.0_KDR

  !     if ( .not. IsProperCell ( iV ) ) &
  !       cycle

  !     if ( N ( iV ) >= N_Trap ) then

  !       Eta  =  ( Mu_E ( iV ) - Mu_NP ( iV ) )  /  ( k * T ( iV ) )
    
  !       F_3  =  ( Eta ** 4  +  a * Eta ** 2  +  b )  /  4.0_KDR  &
  !               -  c * exp ( -Eta )

  !       P_Nu ( iV )  =  d  *  ( k * T ( iV ) ) ** 4  *  F_3

  !     end if

  !   end do
  !   !$OMP end parallel do

  ! end subroutine ComputeNeutrinoPressure


end module ApplyDeleptonization_F__Command
