module PrepareRelaxation_RM__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    PrepareRelaxation_RM

    private :: &
      PrepareRelaxation_RM_Kernel

contains


  subroutine PrepareRelaxation_RM &
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
    class ( GeometryFlatForm ), pointer :: &
      G

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

    G => Chart % Geometry ( )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )
      call PrepareRelaxation_RM_Kernel &
             ( IncrementExplicit % Value ( :, iEnergy ), &
               IncrementExplicit % Value ( :, iMomentum_1 ), &
               IncrementExplicit % Value ( :, iMomentum_2 ), &
               IncrementExplicit % Value ( :, iMomentum_3 ), &
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
               RM % Value ( :, RM % COMOVING_ENERGY ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               RM % Value ( :, RM % CONSERVED_ENERGY ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 1 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 2 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 3 ) ), &
               RM % Value ( :, RM % FLUID_VELOCITY_U ( 1 ) ), &
               RM % Value ( :, RM % FLUID_VELOCITY_U ( 2 ) ), &
               RM % Value ( :, RM % FLUID_VELOCITY_U ( 3 ) ), &
               G % Value ( :, G % METRIC_DD_22 ), &
               G % Value ( :, G % METRIC_DD_33 ), &
               TimeStep, S % B ( iStage ), CONSTANT % SPEED_OF_LIGHT )
    end select !-- SRM

    end associate !-- I

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'PrepareRelaxation_RM__Command', 'module', CONSOLE % ERROR )
      call Show ( 'PrepareRelaxation_RM', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    end select !-- RM

    nullify ( G )
    
  end subroutine PrepareRelaxation_RM
    

  subroutine PrepareRelaxation_RM_Kernel &
               ( KV_E, KV_S_1, KV_S_2, KV_S_3, &
                 DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, &
                 SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3, &
                 IsProperCell, Xi_J, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                 E, S_1, S_2, S_3, V_1, V_2, V_3, M_DD_22, M_DD_33, &
                 dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_E, KV_S_1, KV_S_2, KV_S_3, &
      DCV_E, DCV_S_1, DCV_S_2, DCV_S_3, &
      SVNE_E, SVNE_S_1, SVNE_S_2, SVNE_S_3
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, Chi_J, &
      Chi_H, &
      J, H_1, H_2, H_3, &
      E, S_1, S_2, S_3, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33
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

      KV_E ( iV )  &
        =  KV_E ( iV )  &
           +  c * dT  &
              * ( Xi_J ( iV )  &
                  -  Chi_J ( iV )  *  ( J ( iV )  -  E ( iV ) )  &
                  -  Chi_H ( iV )  &
                     *  (                     V_1 ( iV ) * H_1 ( iV )  &
                          +  M_DD_22 ( iV ) * V_2 ( iV ) * H_2 ( iV )  &
                          +  M_DD_33 ( iV ) * V_3 ( iV ) * H_3 ( iV ) ) )

      KV_S_1 ( iV )  &
        =  KV_S_1 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  *  ( H_1 ( iV )  -  S_1 ( iV ) )  &
                  +  V_1 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      KV_S_2 ( iV )  &
        =  KV_S_2 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  &
                     *  ( M_DD_22 ( iV ) * H_2 ( iV )  -  S_2 ( iV ) )  &
                  +  M_DD_22 ( iV )  *  V_2 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      KV_S_3 ( iV )  &
        =  KV_S_3 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  &
                     *  ( M_DD_33 ( iV ) * H_3 ( iV )  -  S_3 ( iV ) )  &
                  +  M_DD_33 ( iV )  *  V_3 ( iV )  &
                     *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

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

      !-- Approximate since J, H_1, etc. do not include implicit update

      SVNE_E ( iV )  &
        =  SVNE_E ( iV )  &
           +  Weight_RK * c  &
              * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) &
                  -  Chi_H ( iV )  &
                     *  (                     V_1 ( iV ) * H_1 ( iV )  &
                          +  M_DD_22 ( iV ) * V_2 ( iV ) * H_2 ( iV )  &
                          +  M_DD_33 ( iV ) * V_3 ( iV ) * H_3 ( iV ) ) )
                      
      SVNE_S_1 ( iV )  &
        =  SVNE_S_1 ( iV )  &
           +  Weight_RK * c  &
              * ( -  Chi_H ( iV ) * H_1 ( iV )  &
                  +  V_1 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

      SVNE_S_2 ( iV )  &
        =  SVNE_S_2 ( iV )  &
           +  Weight_RK * c  *  M_DD_22 ( iV )  &
              * ( -  Chi_H ( iV ) * H_2 ( iV )  &
                  +  V_2 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

      SVNE_S_3 ( iV )  &
        =  SVNE_S_3 ( iV )  &
           +  Weight_RK * c  *  M_DD_33 ( iV )  &
              * ( -  Chi_H ( iV ) * H_3 ( iV )  &
                  +  V_3 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

    end do
    !$OMP end parallel do

  end subroutine PrepareRelaxation_RM_Kernel
  

end module PrepareRelaxation_RM__Command
