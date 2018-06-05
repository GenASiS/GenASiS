module ApplyCurvilinear_RM__Command

  use Basics
  use Mathematics
  use RadiationMoments_Form
  use Sources_RM__Form

  implicit none
  private

  public :: &
    ApplyCurvilinear_RM

    private :: &
      ApplyCurvilinear_RM_Kernel

contains


  subroutine ApplyCurvilinear_RM &
               ( S, Sources_RM, Increment, RadiationMoments, TimeStep, &
                 iStage )

    class ( Step_RK_C_ASC_Template ), intent ( in ) :: &
      S
    class ( Sources_C_Form ), intent ( inout ) :: &
      Sources_RM
    type ( VariableGroupForm ), intent ( inout ), target :: &
      Increment
    class ( CurrentTemplate ), intent ( in ) :: &
      RadiationMoments
    real ( KDR ), intent ( in ) :: &
      TimeStep
    integer ( KDI ), intent ( in ) :: &
      iStage

    integer ( KDI ) :: &
      iMomentum_1, &
      iMomentum_2
    type ( TimerForm ), pointer :: &
      Timer
    class ( GeometryFlatForm ), pointer :: &
      G

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( RM => RadiationMoments )
    class is ( RadiationMomentsForm )

    call Search &
           ( RM % iaConserved, &
             RM % CONSERVED_MOMENTUM_D ( 1 ), iMomentum_1 )
    call Search &
           ( RM % iaConserved, &
             RM % CONSERVED_MOMENTUM_D ( 2 ), iMomentum_2 )

    select type ( SRM => Sources_RM )
    class is ( Sources_RM_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    if ( trim ( Chart % CoordinateSystem ) == 'RECTANGULAR' ) &
      return

    G => Chart % Geometry ( )
    
    associate &
      ( M_DD_22 => G % Value ( :, G % METRIC_DD_22 ), &
        M_DD_33 => G % Value ( :, G % METRIC_DD_33 ) )
    
    call ApplyCurvilinear_RM_Kernel &
           ( Increment % Value ( :, iMomentum_1 ), &
             Increment % Value ( :, iMomentum_2 ), &
             SRM % Value ( :, SRM % CURVILINEAR_S_D ( 1 ) ), &
             SRM % Value ( :, SRM % CURVILINEAR_S_D ( 2 ) ), &
             Chart % CoordinateSystem, Chart % IsProperCell, &
             RM % Value ( :, RM % COMOVING_ENERGY ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 1 ) ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 2 ) ), &
             RM % Value ( :, RM % COMOVING_MOMENTUM_U ( 3 ) ), &
             RM % Value ( :, RM % STRESS_FACTOR ), &
             M_DD_22, M_DD_33, &
             S % dLogVolumeJacobian_dX ( 1 ) % Value, &
             S % dLogVolumeJacobian_dX ( 2 ) % Value, &
             RM % Value ( :, RM % FLUID_VELOCITY_U ( 2 ) ), &
             RM % Value ( :, RM % FLUID_VELOCITY_U ( 3 ) ), &
             TimeStep, S % B ( iStage ), Chart % nDimensions )

    end associate !-- M_DD_22, etc.
    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'ApplyCurvilinear_RM__Command', 'module', CONSOLE % ERROR )
      call Show ( 'ApplyCurvilinear_RM', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    end select !-- SRM
    end select !-- RM

    nullify ( G )
    
    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyCurvilinear_RM


  subroutine ApplyCurvilinear_RM_Kernel &
               ( KVM_1, KVM_2, SVM_1, SVM_2, CoordinateSystem, IsProperCell, &
                 J, H_1, H_2, H_3, SF, M_DD_22, M_DD_33, &
                 dLVJ_dX1, dLVJ_dX2, V_2, V_3, dT, Weight_RK, nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2, &
      SVM_1, SVM_2
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      J, &
      H_1, H_2, H_3, &
      SF, & 
      M_DD_22, M_DD_33, &
      dLVJ_dX1, dLVJ_dX2, &
      V_2, V_3
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      K_22, K_33, &
      P_22, P_33
    real ( KDR ), dimension ( : ), allocatable :: &
      H

    nV = size ( KVM_1 )

    allocate ( H ( size ( KVM_1 ) ) )
    
    !$OMP parallel do private ( iV )
    do iV = 1, nV
      H ( iV )  =  sqrt ( H_1 ( iV ) ** 2  &
                                +  M_DD_22 ( iV )  *  H_2 ( iV ) ** 2  &
                                +  M_DD_33 ( iV )  *  H_3 ( iV ) ** 2 ) 
    end do
    !$OMP end parallel do
   
    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV, K_33, P_33 )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) cycle

        if ( H ( iV ) > 0.0_KDR ) then
          K_33  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV ) &
                   + 0.5_KDR * ( 3 * SF ( iV ) - 1.0_KDR ) * M_DD_33 ( iV ) &
                     *  H_3 ( iV ) ** 2 / H ( iV )
        else
          K_33  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV )
        end if
        P_33  =  K_33  +  2  *  M_DD_33 ( iV )  *  H_3 ( iV )  *  V_3 ( iV )  

        KVM_1 ( iV )  =  KVM_1 ( iV )  +  P_33 * dLVJ_dX1 ( iV ) * dT
        SVM_1 ( iV )  =  SVM_1 ( iV )  +  Weight_RK * K_33 * dLVJ_dX1 ( iV )
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )

      !$OMP parallel do private ( iV, K_22, K_33, P_22, P_33 )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) cycle
          
        if ( H ( iV ) > 0.0_KDR ) then
          K_22  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV ) &
                   + 0.5_KDR * ( 3 * SF ( iV ) - 1.0_KDR ) * M_DD_22 ( iV ) &
                     *  H_2 ( iV ) ** 2 / H ( iV )
          K_33  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV ) &
                   + 0.5_KDR * ( 3 * SF ( iV ) - 1.0_KDR ) * M_DD_33 ( iV ) &
                     *  H_3 ( iV ) ** 2 / H ( iV )
        else
          K_22  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV )
          K_33  =  0.5_KDR * ( 1.0_KDR - SF ( iV ) ) * J ( iV )
        end if
        P_22  =  K_22  +  2  *  M_DD_22 ( iV )  *  H_2 ( iV )  *  V_2 ( iV )  
        P_33  =  K_33  +  2  *  M_DD_33 ( iV )  *  H_3 ( iV )  *  V_3 ( iV )  

        KVM_1 ( iV ) &
          =  KVM_1 ( iV )  +  0.5_KDR * ( P_22 + P_33 ) * dLVJ_dX1 ( iV ) * dT
        SVM_1 ( iV ) &
          =  SVM_1 ( iV )  &
             +  Weight_RK * 0.5_KDR * ( P_22 + P_33 ) * dLVJ_dX1 ( iV )
        if ( nDimensions > 1 ) then
          KVM_2 ( iV ) &
            =  KVM_2 ( iV )  +  P_33 * dLVJ_dX2 ( iV ) * dT
          SVM_2 ( iV ) &
            =  SVM_2 ( iV )  +  Weight_RK * P_33 * dLVJ_dX2 ( iV )
        end if
      end do
      !$OMP end parallel do
    end select !-- CoordinateSystem

  end subroutine ApplyCurvilinear_RM_Kernel


end module ApplyCurvilinear_RM__Command
