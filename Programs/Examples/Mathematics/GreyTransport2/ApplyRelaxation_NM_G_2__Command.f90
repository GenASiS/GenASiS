module ApplyRelaxation_NM_G_2__Command

  use Basics
  use Mathematics
  use NeutrinoMoments_G__Form
  use Sources_RM__Form
  use ApplyRelaxation_RM_2__Command

  implicit none
  private

  public :: &
    ApplyRelaxation_NM_G_2

    private :: &
      ComputeIncrements

contains

  
  subroutine ApplyRelaxation_NM_G_2 &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 NeutrinoMoments_G, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( inout ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      NeutrinoMoments_G
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iV, &
      nV, &
      iNumber

    call Show ( 'ApplyRelaxation_NM_G_2', S % IGNORABILITY + 3 )
    call Show ( NeutrinoMoments_G % Name, 'NeutrinoMoments_G', &
                S % IGNORABILITY + 3 )

    !-- Energy and Momentum

    call ApplyRelaxation_RM_2 &
           ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
             NeutrinoMoments_G, Chart, TimeStep, iStage )

    !-- Number

    nV = size ( IncrementExplicit % Value, dim = 1 )

    select type ( RM => NeutrinoMoments_G )
    class is ( NeutrinoMoments_G_Form )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    call Search ( RM % iaConserved, RM % CONSERVED_NUMBER, iNumber )

    !$OMP parallel do private ( iV )
    do iV = 1, nV

      if ( .not. Chart % IsProperCell ( iV ) ) &
        cycle

      !-- FIXME: total hack, use Sources_RM to accessibly store comoving 
      !          increments 
      call ComputeIncrements &
             ( IncrementExplicit % Value ( iV, iNumber ), &
               SRM % Value ( iV, SRM % NET_EMISSION_N ), &
               I  % Value ( iV, I % EMISSIVITY_N ), &
               I  % Value ( iV, I % OPACITY_N ), &
               RM % Value ( iV, RM % FLUID_VELOCITY_U ( 1 ) ), &
               RM % Value ( iV, RM % COMOVING_NUMBER ), &
               RM % Value ( iV, RM % ENERGY_AVERAGE ), &
               SRM % Value ( iV, SRM % NET_EMISSION_S_D ( 1 ) ), &
               TimeStep )

    end do
    !$OMP end parallel do

    end associate !-- I
    end select !-- Chart
    end select !-- SRM
    end select !-- RM

  end subroutine ApplyRelaxation_NM_G_2


  subroutine ComputeIncrements &
               ( dD, SVNE_N, Xi_N, Chi_N, V_1, N, E_Ave, dH_1, dt )

    real ( KDR ), intent ( inout ) :: &
      dD, &
      SVNE_N
    real ( KDR ), intent ( in ) :: &
      Xi_N, Chi_N, &
      V_1, N, E_Ave, dH_1, &
      dt

    real ( KDR ) :: &
      dN

    dD  =  ( dD  +  ( Xi_N  -  Chi_N * N ) * dt )  &
           /  ( 1.0_KDR  +  Chi_N * dt )

    if ( E_Ave > 0.0_KDR ) then
      dN  =  dD  -  V_1 * dH_1 / E_Ave
    else
      dN  =  dD
    end if

    !-- FIXME: total hack, use this to accessibly store comoving increments 
    SVNE_N  =  dN

  end subroutine ComputeIncrements


end module ApplyRelaxation_NM_G_2__Command
