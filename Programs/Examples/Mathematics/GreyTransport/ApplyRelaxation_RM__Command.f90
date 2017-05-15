module ApplyRelaxation_RM__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    ApplyRelaxation_RM

    private :: &
      ApplyRelaxation_RM_Kernel

contains


  subroutine ApplyRelaxation_RM &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 RadiationMoments, Chart, TimeStep, iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ) :: &
      IncrementExplicit, &
      DampingCoefficient
    class ( CurrentTemplate ), intent ( in ) :: &
      RadiationMoments
    class ( ChartTemplate ), intent ( in ) :: &
      Chart
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iEnergy, &
      iMomentum_1, &
      iMomentum_2, &
      iMomentum_3

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY_DENSITY, &
                  iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 1 ), &
                  iMomentum_1 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 2 ), &
                  iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_DENSITY_D ( 3 ), &
                  iMomentum_3 )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )
      call ApplyRelaxation_RM_Kernel &
             ( IncrementExplicit % Value ( :, iEnergy ), &
               DampingCoefficient % Value ( :, iEnergy ), &
               DampingCoefficient % Value ( :, iMomentum_1 ), &
               DampingCoefficient % Value ( :, iMomentum_2 ), &
               DampingCoefficient % Value ( :, iMomentum_3 ), &
               SRM % Value ( :, SRM % EMISSION_E ), &
               SRM % Value ( :, SRM % EMISSION_S_D ( 1 ) ), &
               SRM % Value ( :, SRM % EMISSION_S_D ( 2 ) ), &
               SRM % Value ( :, SRM % EMISSION_S_D ( 3 ) ), &
               SRM % Value ( :, SRM % ABSORPTION_E ), &
               SRM % Value ( :, SRM % ABSORPTION_S_D ( 1 ) ), &
               SRM % Value ( :, SRM % ABSORPTION_S_D ( 2 ) ), &
               SRM % Value ( :, SRM % ABSORPTION_S_D ( 3 ) ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EQUILIBRIUM_DENSITY ), &
               I % Value ( :, I % EFFECTIVE_OPACITY ), &
               I % Value ( :, I % TRANSPORT_OPACITY ), &
               RM % Value ( :, RM % CONSERVED_ENERGY_DENSITY ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_DENSITY_D ( 1 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_DENSITY_D ( 2 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_DENSITY_D ( 3 ) ), &
               TimeStep, S % B ( iStage ), CONSTANT % SPEED_OF_LIGHT )
    end select !-- SRM

    end associate !-- I

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'ApplyRelaxation_RM__Command', 'module', CONSOLE % ERROR )
      call Show ( 'ApplyRelaxation_RM', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

     end select !-- RM
    
  end subroutine ApplyRelaxation_RM
    

  subroutine ApplyRelaxation_RM_Kernel &
               ( KV_E, DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, &
                 SVE_E, SVE_S_1, SVE_S_2, SVE_S_3, &
                 SVA_E, SVA_S_1, SVA_S_2, SVA_S_3, &
                 IsProperCell, ED, EO, TO, E, S_1, S_2, S_3, dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_E, &
      DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, &
      SVE_E, SVE_S_1, SVE_S_2, SVE_S_3, &
      SVA_E, SVA_S_1, SVA_S_2, SVA_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      ED, &
      EO, &
      TO, &
      E, &
      S_1, S_2, S_3
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK, &
      c

    integer ( KDI ) :: &
      iV, &
      nV

    nV = size ( KV_E )

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle

      KV_E    ( iV )  =  KV_E ( iV )  +  c * EO ( iV ) * ED ( iV ) * dT
      DCV_E   ( iV )  =  c * EO ( iV )
      DCV_S_1 ( iV )  =  c * TO ( iV )
      DCV_S_2 ( iV )  =  c * TO ( iV )
      DCV_S_3 ( iV )  =  c * TO ( iV )

    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle

      SVE_E ( iV )  &
        =  SVE_E ( iV )  +  Weight_RK * c * EO ( iV ) * ED ( iV )

      !-- Approximate since E, S_1, etc. do not include implicit update
      SVA_E ( iV )  &
        =  SVA_E ( iV )  -  Weight_RK * c * EO ( iV ) * E ( iV )
      SVA_S_1 ( iV )  &
        =  SVA_S_1 ( iV )  -  Weight_RK * c * TO ( iV ) * S_1 ( iV )
      SVA_S_2 ( iV )  &
        =  SVA_S_2 ( iV )  -  Weight_RK * c * TO ( iV ) * S_2 ( iV )
      SVA_S_2 ( iV )  &
        =  SVA_S_2 ( iV )  -  Weight_RK * c * TO ( iV ) * S_3 ( iV )

    end do
    !$OMP end parallel do

  end subroutine ApplyRelaxation_RM_Kernel
  

end module ApplyRelaxation_RM__Command
