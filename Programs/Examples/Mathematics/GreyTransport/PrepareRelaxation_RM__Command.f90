module PrepareRelaxation_RM__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    PrepareRelaxation_RM_Energy, &
    PrepareRelaxation_RM_Momentum

    private :: &
      PrepareRelaxation_RM_Energy_Kernel, &
      PrepareRelaxation_RM_Momentum_Kernel

contains


  subroutine PrepareRelaxation_RM_Energy &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 RadiationMoments, Chart, TimeStep, iStage, iEnergy )

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
    integer ( KDI ), intent ( out ) :: &
      iEnergy

    class ( GeometryFlatForm ), pointer :: &
      G

    call Show ( 'PrepareRelaxation_RM_Energy', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, iEnergy )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    G => Chart % Geometry ( )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )
      call PrepareRelaxation_RM_Energy_Kernel &
             ( IncrementExplicit % Value ( :, iEnergy ), &
               DampingCoefficient % Value ( :, iEnergy ), &
               SRM % Value ( :, SRM % INTERACTIONS_J ), &
               Chart % IsProperCell, &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               RM % Value ( :, RM % COMOVING_ENERGY ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               RM % Value ( :, RM % CONSERVED_ENERGY ), &
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
    
  end subroutine PrepareRelaxation_RM_Energy
    

  subroutine PrepareRelaxation_RM_Momentum &
               ( S, Sources_RM, IncrementExplicit, DampingCoefficient, &
                 RadiationMoments, Chart, TimeStep, iStage, &
                 iMomentum_1, iMomentum_2, iMomentum_3 )

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
    integer ( KDI ), intent ( out ) :: &
      iMomentum_1, &
      iMomentum_2, &
      iMomentum_3

    integer ( KDI ) :: &
      iEnergy
    class ( GeometryFlatForm ), pointer :: &
      G

    call Show ( 'PrepareRelaxation_RM_Momentum', S % IGNORABILITY + 3 )
    call Show ( RadiationMoments % Name, 'RadiationMoments', &
                S % IGNORABILITY + 3 )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 1 ), &
                  iMomentum_1 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 2 ), &
                  iMomentum_2 )
    call Search ( RM % iaConserved, RM % CONSERVED_MOMENTUM_D ( 3 ), &
                  iMomentum_3 )
    call Search ( RM % iaConserved, RM % CONSERVED_ENERGY, &
                  iEnergy )

    select type ( Chart )
    class is ( Chart_SL_Template )

    associate ( I => RM % Interactions )

    G => Chart % Geometry ( )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )
      call PrepareRelaxation_RM_Momentum_Kernel &
             ( IncrementExplicit % Value ( :, iMomentum_1 ), &
               IncrementExplicit % Value ( :, iMomentum_2 ), &
               IncrementExplicit % Value ( :, iMomentum_3 ), &
               DampingCoefficient % Value ( :, iMomentum_1 ), &
               DampingCoefficient % Value ( :, iMomentum_2 ), &
               DampingCoefficient % Value ( :, iMomentum_3 ), &
               SRM % Value ( :, SRM % INTERACTIONS_H_D ( 1 ) ), &
               SRM % Value ( :, SRM % INTERACTIONS_H_D ( 2 ) ), &
               SRM % Value ( :, SRM % INTERACTIONS_H_D ( 3 ) ), &
               RM, &
               Chart % IsProperCell, &
               IncrementExplicit % Value ( :, iEnergy ), &
               I % Value ( :, I % EMISSIVITY_J ), &
               I % Value ( :, I % OPACITY_J ), &
               I % Value ( :, I % OPACITY_H ), &
               RM % Value ( :, RM % COMOVING_ENERGY ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
               RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 1 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 2 ) ), &
               RM % Value ( :, RM % CONSERVED_MOMENTUM_D ( 3 ) ), &
               RM % Value ( :, RM % FLUX_FACTOR ), &
               RM % Value ( :, RM % STRESS_FACTOR ), &
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
    
  end subroutine PrepareRelaxation_RM_Momentum
    

  subroutine PrepareRelaxation_RM_Energy_Kernel &
               ( KV_E, DCV_E, SVNE_E, &
                 IsProperCell, Xi_J, Chi_J, Chi_H, J, H_1, H_2, H_3, &
                 E, V_1, V_2, V_3, M_DD_22, M_DD_33, &
                 dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_E, &
      DCV_E, &
      SVNE_E
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      Xi_J, Chi_J, &
      Chi_H, &
      J, H_1, H_2, H_3, &
      E, &
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

      DCV_E ( iV )  =  c * Chi_J ( iV )

    end do
    !$OMP end parallel do

    ! !$OMP parallel do private ( iV )
    ! do iV = 1, nV
    !   if ( .not. IsProperCell ( iV ) ) &
    !     cycle

    !   !-- Approximate since J, H_1, etc. do not include implicit update

    !   SVNE_E ( iV )  &
    !     =  SVNE_E ( iV )  &
    !        +  Weight_RK * c  &
    !           * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) &
    !               -  Chi_H ( iV )  &
    !                  *  (                     V_1 ( iV ) * H_1 ( iV )  &
    !                       +  M_DD_22 ( iV ) * V_2 ( iV ) * H_2 ( iV )  &
    !                       +  M_DD_33 ( iV ) * V_3 ( iV ) * H_3 ( iV ) ) )
                      
    ! end do
    ! !$OMP end parallel do

  end subroutine PrepareRelaxation_RM_Energy_Kernel
  

  subroutine PrepareRelaxation_RM_Momentum_Kernel &
               ( KV_S_1, KV_S_2, KV_S_3, &
                 DCV_S_1, DCV_S_2, DCV_S_3, &
                 SVNE_S_1, SVNE_S_2, SVNE_S_3, &
                 RM, IsProperCell, KV_E, Xi_J, Chi_J, Chi_H, &
                 J, H_1, H_2, H_3, S_1, S_2, S_3, FF, SF, &
                 V_1, V_2, V_3, M_DD_22, M_DD_33, &
                 dT, Weight_RK, c )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KV_S_1, KV_S_2, KV_S_3, &
      DCV_S_1, DCV_S_2, DCV_S_3, &
      SVNE_S_1, SVNE_S_2, SVNE_S_3
    class ( RadiationMomentsForm ), intent ( in ) :: &
      RM
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      KV_E, & 
      Xi_J, Chi_J, &
      Chi_H, &
      J, H_1, H_2, H_3, &
      S_1, S_2, S_3, &
      FF, SF, &
      V_1, V_2, V_3, &
      M_DD_22, M_DD_33
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK, &
      c

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ), dimension ( 3, 3 ) :: &
      K_U_D

    nV = size ( KV_E )

    !$OMP parallel do private ( iV, K_U_D )
    do iV = 1, nV
      if ( .not. IsProperCell ( iV ) ) &
        cycle

!       KV_S_1 ( iV )  &
!         =  KV_S_1 ( iV )  &
!            +  c * dT  &
!               * ( -  Chi_H ( iV )  *  ( H_1 ( iV )  -  S_1 ( iV ) )  &
!                   +  V_1 ( iV )  &
!                      *  ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) ) 

      call RM % ComputeComovingStress_D &
             ( K_U_D ( 1, : ), [ H_1 ( iV ), H_2 ( iV ), H_3 ( iV ) ], &
               J ( iV ), SF ( iV ), M_DD_22 ( iV ), M_DD_33 ( iV ), iD = 1 )
      call RM % ComputeComovingStress_D &
             ( K_U_D ( 2, : ), [ H_1 ( iV ), H_2 ( iV ), H_3 ( iV ) ], &
               J ( iV ), SF ( iV ), M_DD_22 ( iV ), M_DD_33 ( iV ), iD = 2 )
      call RM % ComputeComovingStress_D &
             ( K_U_D ( 3, : ), [ H_1 ( iV ), H_2 ( iV ), H_3 ( iV ) ], &
               J ( iV ), SF ( iV ), M_DD_22 ( iV ), M_DD_33 ( iV ), iD = 3 )

      KV_S_1 ( iV )  &
        =  KV_S_1 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  &
                     *  ( - V_1 ( iV )  &  
                          - ( K_U_D ( 1, 1 )  *  V_1 ( iV )  &
                              +  M_DD_22 ( iV )  &
                                 *  K_U_D ( 2, 1 )  *  V_2 ( iV )  &
                              +  M_DD_33 ( iV )  &
                                 *  K_U_D ( 3, 1 )  *  V_3 ( iV ) )  &
                            /  max ( J ( iV ), tiny ( 0.0_KDR ) )  )  &
                        *  ( J ( iV )  +  KV_E ( iV ) )  &
                  +  V_1 ( iV )  &
                     *  ( Xi_J ( iV )  &
                          -  Chi_J ( iV ) * ( J ( iV )  +  KV_E ( iV ) ) ) ) 

      KV_S_2 ( iV )  &
        =  KV_S_2 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  &
                     *  ( - M_DD_22 ( iV )  *  V_2 ( iV )  &
                          - ( K_U_D ( 1, 2 )  *  V_1 ( iV )  &
                              +  M_DD_22 ( iV )  &
                                 *  K_U_D ( 2, 2 )  *  V_2 ( iV )  &
                              +  M_DD_33 ( iV )  &
                                 *  K_U_D ( 3, 2 )  *  V_3 ( iV ) )  &
                            /  max ( J ( iV ), tiny ( 0.0_KDR ) )  )  &
                        *  ( J ( iV )  +  KV_E ( iV ) )  &
                  +  M_DD_22 ( iV )  *  V_2 ( iV )  &
                     *  ( Xi_J ( iV )  &
                          -  Chi_J ( iV ) * ( J ( iV )  +  KV_E ( iV ) ) ) ) 

      KV_S_3 ( iV )  &
        =  KV_S_3 ( iV )  &
           +  c * dT  &
              * ( -  Chi_H ( iV )  &
                     *  ( - M_DD_33 ( iV )  *  V_3 ( iV )  &
                          - ( K_U_D ( 1, 3 )  *  V_1 ( iV )  &
                              +  M_DD_22 ( iV )  &
                                 *  K_U_D ( 2, 3 )  *  V_2 ( iV )  &
                              +  M_DD_33 ( iV )  &
                                 *  K_U_D ( 3, 3 )  *  V_3 ( iV ) )  &
                            /  max ( J ( iV ), tiny ( 0.0_KDR ) )  )  &
                        *  ( J ( iV )  +  KV_E ( iV ) )  &
                  +  M_DD_33 ( iV )  *  V_3 ( iV )  &
                     *  ( Xi_J ( iV )  &
                          -  Chi_J ( iV ) * ( J ( iV )  +  KV_E ( iV ) ) ) ) 

      DCV_S_1 ( iV )  =  c * Chi_H ( iV )
      DCV_S_2 ( iV )  =  c * Chi_H ( iV )
      DCV_S_3 ( iV )  =  c * Chi_H ( iV )

    end do
    !$OMP end parallel do

    ! !$OMP parallel do private ( iV )
    ! do iV = 1, nV
    !   if ( .not. IsProperCell ( iV ) ) &
    !     cycle

    !   !-- Approximate since J, H_1, etc. do not include implicit update

    !   SVNE_S_1 ( iV )  &
    !     =  SVNE_S_1 ( iV )  &
    !        +  Weight_RK * c  &
    !           * ( -  Chi_H ( iV ) * H_1 ( iV )  &
    !               +  V_1 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

    !   SVNE_S_2 ( iV )  &
    !     =  SVNE_S_2 ( iV )  &
    !        +  Weight_RK * c  *  M_DD_22 ( iV )  &
    !           * ( -  Chi_H ( iV ) * H_2 ( iV )  &
    !               +  V_2 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

    !   SVNE_S_3 ( iV )  &
    !     =  SVNE_S_3 ( iV )  &
    !        +  Weight_RK * c  *  M_DD_33 ( iV )  &
    !           * ( -  Chi_H ( iV ) * H_3 ( iV )  &
    !               +  V_3 ( iV ) * ( Xi_J ( iV )  -  Chi_J ( iV ) * J ( iV ) ) )

    ! end do
    ! !$OMP end parallel do

  end subroutine PrepareRelaxation_RM_Momentum_Kernel
  

end module PrepareRelaxation_RM__Command
