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

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, &
                  iEnergy )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 2 ), &
                  iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 3 ), &
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
               SRM % Value ( :, SRM % NET_EMISSION_E ), &
               SRM % Value ( :, SRM % NET_EMISSION_S_D ( 1 ) ), &
               SRM % Value ( :, SRM % NET_EMISSION_S_D ( 2 ) ), &
               SRM % Value ( :, SRM % NET_EMISSION_S_D ( 3 ) ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               RM % Value ( :, RM % CONSERVED_ENERGY ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 1 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 2 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 3 ) ), &
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
                 SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3, &
                 IsProperCell, Xi_J, Chi_J, Chi_H, E, S_1, S_2, S_3, &
                 dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_E, &
      DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, &
      SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, &
      Chi_J, &
      Chi_H, &
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

      KV_E    ( iV )  =  KV_E ( iV )  +  c * Xi_J ( iV ) * dT
      DCV_E   ( iV )  =  c * Chi_J ( iV )
      DCV_S_1 ( iV )  =  c * Chi_H ( iV )
      DCV_S_2 ( iV )  =  c * Chi_H ( iV )
      DCV_S_3 ( iV )  =  c * Chi_H ( iV )

    end do
    !$OMP end parallel do

    !$OMP parallel do private ( iV )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle

      !-- Approximate since E, S_1, etc. do not include implicit update

      SVNE_E ( iV )  &
        =  SVNE_E ( iV )  &
           +  Weight_RK * c * ( Xi_J ( iV )  -  Chi_J ( iV ) * E ( iV ) )

      SVNE_S_1 ( iV )  &
        =  SVNE_S_1 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_1 ( iV )
      SVNE_S_2 ( iV )  &
        =  SVNE_S_2 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_2 ( iV )
      SVNE_S_2 ( iV )  &
        =  SVNE_S_2 ( iV )  -  Weight_RK * c * Chi_H ( iV ) * S_3 ( iV )

    end do
    !$OMP end parallel do

  end subroutine ApplyRelaxation_RM_Kernel
  

end module ApplyRelaxation_RM__Command
