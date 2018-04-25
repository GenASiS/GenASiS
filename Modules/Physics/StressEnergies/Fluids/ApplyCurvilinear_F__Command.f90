module ApplyCurvilinear_F__Command

  use Basics
  use Mathematics
  use Fluid_D__Form
  use Fluid_P__Template
  use Sources_F__Form

  implicit none
  private

  public :: &
    ApplyCurvilinear_F

    private :: &
      ApplyCurvilinear_F_D_Kernel, &
      ApplyCurvilinear_F_P_Kernel

contains


  subroutine ApplyCurvilinear_F &
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
      iMomentum_1, &
      iMomentum_2
    type ( TimerForm ), pointer :: &
      Timer

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( F => Fluid )
    class is ( Fluid_D_Form )
      call Search &
             ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 1 ), iMomentum_1 )
      call Search &
             ( F % iaConserved, F % MOMENTUM_DENSITY_D ( 2 ), iMomentum_2 )
    end select !-- F

    select type ( S_F => Sources_F )
    class is ( Sources_F_Form )

    select type ( Chart => S % Chart )
    class is ( Chart_SL_Template )

    if ( trim ( Chart % CoordinateSystem ) == 'CARTESIAN' ) &
      return

    select type ( F => Fluid )
    class is ( Fluid_D_Form )
      call ApplyCurvilinear_F_D_Kernel &
             ( Increment % Value ( :, iMomentum_1 ), &
               Increment % Value ( :, iMomentum_2 ), &
               S_F % Value ( :, S_F % CURVILINEAR_S_D ( 1 ) ), &
               S_F % Value ( :, S_F % CURVILINEAR_S_D ( 2 ) ), &
               Chart % CoordinateSystem, Chart % IsProperCell, &
               F % Value ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), &
               F % Value ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), &
               F % Value ( :, F % VELOCITY_U ( 2 ) ), &
               F % Value ( :, F % VELOCITY_U ( 3 ) ), &
               S % dLogVolumeJacobian_dX ( 1 ) % Value, &
               S % dLogVolumeJacobian_dX ( 2 ) % Value, &
               TimeStep, S % B ( iStage ), Chart % nDimensions )
    class is ( Fluid_P_Template )
      call ApplyCurvilinear_F_P_Kernel &
             ( Increment % Value ( :, iMomentum_1 ), &
               Increment % Value ( :, iMomentum_2 ), &
               S_F % Value ( :, S_F % CURVILINEAR_S_D ( 1 ) ), &
               S_F % Value ( :, S_F % CURVILINEAR_S_D ( 2 ) ), &
               Chart % CoordinateSystem, Chart % IsProperCell, &
               F % Value ( :, F % PRESSURE ), &
               F % Value ( :, F % MOMENTUM_DENSITY_D ( 2 ) ), &
               F % Value ( :, F % MOMENTUM_DENSITY_D ( 3 ) ), &
               F % Value ( :, F % VELOCITY_U ( 2 ) ), &
               F % Value ( :, F % VELOCITY_U ( 3 ) ), &
               S % dLogVolumeJacobian_dX ( 1 ) % Value, &
               S % dLogVolumeJacobian_dX ( 2 ) % Value, &
               TimeStep, S % B ( iStage ), Chart % nDimensions )
    end select !-- F

    class default
      call Show ( 'Chart type not found', CONSOLE % ERROR )
      call Show ( 'ApplyCurvilinear_F__Command', 'module', CONSOLE % ERROR )
      call Show ( 'ApplyCurvilinear_F', 'subroutine', CONSOLE % ERROR ) 
      call PROGRAM_HEADER % Abort ( )
    end select !-- Chart

    end select !-- S_F

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyCurvilinear_F


  subroutine ApplyCurvilinear_F_D_Kernel &
               ( KVM_1, KVM_2, SVM_1, SVM_2, CoordinateSystem, IsProperCell, &
                 S_2, S_3, V_2, V_3, dLVJ_dX1, dLVJ_dX2, dT, Weight_RK, &
                 nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2, &
      SVM_1, SVM_2
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      S_2, S_3, &
      V_2, V_3, &
      dLVJ_dX1, dLVJ_dX2
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Curvilinear

    nV = size ( KVM_1 )

    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV, Curvilinear )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        Curvilinear   =  ( V_3 ( iV ) * S_3 ( iV ) )  *  dLVJ_dX1 ( iV )
        KVM_1 ( iV )  =  KVM_1 ( iV )  +  Curvilinear * dT
        SVM_1 ( iV )  =  SVM_1 ( iV )  +  Weight_RK * Curvilinear
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )

      !$OMP parallel do private ( iV, Curvilinear )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        Curvilinear  &
          =  ( V_2 ( iV ) * S_2 ( iV )  +  V_3 ( iV ) * S_3 ( iV ) ) &
             * 0.5_KDR * dLVJ_dX1 ( iV )
        KVM_1 ( iV )  =  KVM_1 ( iV )  +  Curvilinear * dT
        SVM_1 ( iV )  =  SVM_1 ( iV )  +  Weight_RK * Curvilinear
      end do
      !$OMP end parallel do

      if ( nDimensions > 1 ) then
        !$OMP parallel do private ( iV, Curvilinear )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          Curvilinear  =  ( V_3 ( iV ) * S_3 ( iV ) ) * dLVJ_dX2 ( iV )
          KVM_2 ( iV )  =  KVM_2 ( iV )  +  Curvilinear * dT
          SVM_2 ( iV )  =  SVM_2 ( iV )  +  Weight_RK * Curvilinear
        end do
        !$OMP end parallel do
      end if

    end select !-- CoordinateSystem

  end subroutine ApplyCurvilinear_F_D_Kernel


  subroutine ApplyCurvilinear_F_P_Kernel &
               ( KVM_1, KVM_2, SVM_1, SVM_2, CoordinateSystem, IsProperCell, &
                 P, S_2, S_3, V_2, V_3, dLVJ_dX1, dLVJ_dX2, dT, Weight_RK, &
                 nDimensions )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      KVM_1, KVM_2, &
      SVM_1, SVM_2
    character ( * ), intent ( in ) :: &
      CoordinateSystem
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      P, &
      S_2, S_3, &
      V_2, V_3, &
      dLVJ_dX1, dLVJ_dX2
    real ( KDR ), intent ( in ) :: &
      dT, &
      Weight_RK
    integer ( KDI ), intent ( in ) :: &
      nDimensions

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      Curvilinear

    nV = size ( KVM_1 )

    select case ( trim ( CoordinateSystem ) )
    case ( 'CYLINDRICAL' )

      !$OMP parallel do private ( iV, Curvilinear )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        Curvilinear &
          =  ( V_3 ( iV ) * S_3 ( iV )  +  P ( iV ) )  *  dLVJ_dX1 ( iV )
        KVM_1 ( iV )  =  KVM_1 ( iV )  +  Curvilinear * dT
        SVM_1 ( iV )  =  SVM_1 ( iV )  +  Weight_RK * Curvilinear
      end do
      !$OMP end parallel do

    case ( 'SPHERICAL' )
      !$OMP parallel do private ( iV, Curvilinear )
      do iV = 1, nV
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        Curvilinear &
          =  (    V_2 ( iV ) * S_2 ( iV )  &
               +  V_3 ( iV ) * S_3 ( iV )  &
               +  2.0 * P ( iV ) ) &
             * 0.5_KDR * dLVJ_dX1 ( iV ) 
        KVM_1 ( iV )  =  KVM_1 ( iV )  +  Curvilinear * dT
        SVM_1 ( iV )  =  SVM_1 ( iV )  +  Weight_RK * Curvilinear
      end do
      !$OMP end parallel do

      if ( nDimensions > 1 ) then
        !$OMP parallel do private ( iV, Curvilinear )
        do iV = 1, nV
          if ( .not. IsProperCell ( iV ) ) &
            cycle
          Curvilinear &
            =  ( V_3 ( iV ) * S_3 ( iV )  +  P ( iV ) )  *  dLVJ_dX2 ( iV )
          KVM_2 ( iV )  =  KVM_2 ( iV )  +  Curvilinear * dT
          SVM_2 ( iV )  =  SVM_2 ( iV )  +  Weight_RK * Curvilinear
        end do
        !$OMP end parallel do
      end if

    end select !-- CoordinateSystem

  end subroutine ApplyCurvilinear_F_P_Kernel


end module ApplyCurvilinear_F__Command
