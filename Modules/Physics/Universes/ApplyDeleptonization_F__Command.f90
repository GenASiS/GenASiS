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
      ApplyProtonNumber

    real ( KDR ), private :: &
      Y_1   = 0.5_KDR, &     !-- N13
      Y_2   = 0.285_KDR, &   !-- N13
      Y_C   = 0.035_KDR, &   !-- N13
      Rho_1 = 2.0e7_KDR, &   !-- N13
      Rho_2 = 2.0e13_KDR, &  !-- N13
      ! Y_1   = 0.5_KDR, &     !-- G15
      ! Y_2   = 0.278_KDR, &   !-- G15
      ! Y_C   = 0.035_KDR      !-- G15
      ! Rho_1 = 2.0e7_KDR, &   !-- G15
      ! Rho_2 = 2.0e13_KDR, &  !-- G15
      N_1, N_2

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
      iBaryon
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    N_1  =  Rho_1  *  UNIT % MASS_DENSITY_CGS / UNIT % ATOMIC_MASS_UNIT 
    N_2  =  Rho_2  *  UNIT % MASS_DENSITY_CGS / UNIT % ATOMIC_MASS_UNIT 

    select type ( F => Fluid )
    class is ( Fluid_P_HN_Form )

    select type ( S_F => Sources_F )
    class is ( Sources_F_Form )

    call Search ( F % iaConserved, F % CONSERVED_PROTON_DENSITY, iProton )
    call Search ( F % iaConserved, F % CONSERVED_BARYON_DENSITY, iBaryon )

    if ( iStage == 1 ) then
      call Clear ( S_F % Value ( :, S_F % RADIATION_DP ) )
    end if

    call ApplyProtonNumber &
           ( Increment % Value ( :, iProton ), &
             S_F % Value ( :, S_F % RADIATION_DP ), &
             F % Value ( :, F % COMOVING_BARYON_DENSITY ), &
             F % Value ( :, F % PROTON_FRACTION ), &
             Increment % Value ( :, iBaryon ), & 
             TimeStep, S % B ( iStage ) )

    end select !-- S_F
    end select !-- F

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyDeleptonization_F


  subroutine ApplyProtonNumber &
               ( K_DP, S_DP, N, YP, K_D, dt, Weight_RK )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K_DP, &
      S_DP
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      N, &
      YP, &
      K_D
    real ( KDR ) :: &
      dt, &
      Weight_RK

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    real ( KDR ) :: &
      X, &
      YP_bar, &
      dYP, &
      dDP

    nValues = size ( K_DP )

    !$OMP parallel do private ( iV )
    do iV = 1, nValues

      X  =  abs ( max ( -1.0_KDR, &
                        min ( 1.0_KDR, &
                              ( 2 * log10 ( N ( iV ) ) &
                                - log10 ( N_2 ) - log10 ( N_1 ) ) &
                              / ( log10 ( N_2 ) - log10 ( N_1 ) ) ) ) )

      YP_bar  =  0.5_KDR * ( Y_1 + Y_2 )  +  0.5_KDR * X * ( Y_2 - Y_1 ) &
                  +  Y_C * ( 1.0_KDR - X  +  4.0_KDR * X &
                                             * ( X - 0.5_KDR ) &
                                             * ( X - 1.0_KDR ) )

      dYP  =  min ( 0.0_KDR, YP_bar - YP ( iV ) )
      dDP  =  N ( iV ) * dYP  +  YP ( iV )  *  K_D ( iV )

      K_DP ( iV )  =  K_DP ( iV )  +  dDP
      S_DP ( iV )  =  S_DP ( iV )  +  Weight_RK * dDP / dt

    end do
    !$OMP end parallel do

  end subroutine ApplyProtonNumber


end module ApplyDeleptonization_F__Command
