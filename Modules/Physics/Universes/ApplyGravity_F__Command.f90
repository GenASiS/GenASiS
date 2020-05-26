#include "Preprocessor"

module ApplyGravity_F__Command

  use Basics
  use Mathematics
  use Spaces
  use StressEnergies

  implicit none
  private

  public :: &
    ApplyGravity_F

    private :: &
      ApplyGravityMomentum, &
      ApplyGravityEnergy

contains


  subroutine ApplyGravity_F &
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
      iD, &  !-- iDimension
      iMomentum, &
      iEnergy
    type ( TimerForm ), pointer :: &
      Timer
    class ( Geometry_N_Form ), pointer :: &
      G

    Timer => PROGRAM_HEADER % TimerPointer ( S % iTimerSources )
    if ( associated ( Timer ) ) call Timer % Start ( )

    select type ( PS => S % Current_ASC % Atlas_SC )
    class is ( Atlas_SC_Form )

    select type ( GA => PS % Geometry_ASC )
    class is ( Geometry_ASC_Form )
    G  =>  GA % Geometry_N ( )

    select type ( F => Fluid )
    class is ( Fluid_D_Form )

    select type ( FS => Sources_F )
    class is ( Sources_F_Form )

    select type ( C => PS % Chart )
    class is ( Chart_SLD_Form )
    
    if ( trim ( GA % GeometryType ) == 'NEWTONIAN' ) then
      do iD = 1, C % nDimensions
        call Search ( F % iaConserved, F % MOMENTUM_DENSITY_D ( iD ), &
                      iMomentum )
        if ( iStage == 1 ) &
          call Clear ( FS % Value ( :, FS % GRAVITATIONAL_S_D ( iD ) ), &
                       UseDeviceOption = FS % AllocatedDevice )
        call ApplyGravityMomentum &
               ( Increment % Value ( :, iMomentum ), & 
                 FS % Value ( :, FS % GRAVITATIONAL_S_D ( iD ) ), &
                 C % IsProperCell, &
                 F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
                 G % Value ( :, G % POTENTIAL_GRADIENT_D ( iD ) ), &
                 TimeStep, S % B ( iStage ), &
                 UseDeviceOption = FS % AllocatedDevice )
      end do !-- iD
    end if

    if ( any ( trim ( GA % GeometryType ) &
                 == [ 'NEWTONIAN       ', 'NEWTONIAN_STRESS' ] ) ) &
    then
      select type ( F_P => Fluid )
      class is ( Fluid_P_Template )
        call Search ( F_P % iaConserved, F_P % CONSERVED_ENERGY, iEnergy )
        if ( iStage == 1 ) &
          call Clear ( FS % Value ( :, FS % GRAVITATIONAL_G ), &
                       UseDeviceOption = FS % AllocatedDevice )
        call ApplyGravityEnergy &
               ( Increment % Value ( :, iEnergy ), & 
                 FS % Value ( :, FS % GRAVITATIONAL_G ), &
                 C % IsProperCell, &
                 F % Value ( :, F % BARYON_MASS ), &
                 F % Value ( :, F % CONSERVED_BARYON_DENSITY ), &
                 F % Value ( :, F % VELOCITY_U ( 1 ) ), &
                 F % Value ( :, F % VELOCITY_U ( 2 ) ), &
                 F % Value ( :, F % VELOCITY_U ( 3 ) ), &
                 G % Value ( :, G % POTENTIAL_GRADIENT_D ( 1 ) ), &
                 G % Value ( :, G % POTENTIAL_GRADIENT_D ( 2 ) ), &
                 G % Value ( :, G % POTENTIAL_GRADIENT_D ( 3 ) ), &
                 TimeStep, S % B ( iStage ), &
                 UseDeviceOption = FS % AllocatedDevice )
      end select !-- F_P
    end if
    
    end select !-- C
    end select !-- FS
    end select !-- F
    end select !-- GA
    end select !-- PS
    nullify ( G )

    if ( associated ( Timer ) ) call Timer % Stop ( )

  end subroutine ApplyGravity_F


  subroutine ApplyGravityMomentum &
               ( K, S, IsProperCell, M, N, GradPhi, dt, Weight_RK, &
                 UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K, &
      S
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      GradPhi
    real ( KDR ), intent ( in ) :: &
      dt, &
      Weight_RK
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption
      
    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( K )

    if ( UseDevice ) then
    
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP  schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                  *  dt
        S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                  *  Weight_RK
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    
    else

      !$OMP  parallel do &
      !$OMP  schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                  *  dt
        S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  *  GradPhi ( iV )  &
                                  *  Weight_RK
      end do
      !$OMP  end parallel do
      
    end if

  end subroutine ApplyGravityMomentum


  subroutine ApplyGravityEnergy &
               ( K, S, IsProperCell, M, N, V_1, V_2, V_3, &
                 GradPhi_1, GradPhi_2, GradPhi_3, &
                 dt, Weight_RK, UseDeviceOption )

    real ( KDR ), dimension ( : ), intent ( inout ) :: &
      K, &
      S
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      IsProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      M, &
      N, &
      V_1, V_2, V_3, &
      GradPhi_1, GradPhi_2, GradPhi_3
    real ( KDR ), intent ( in ) :: &
      dt, &
      Weight_RK
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &  !-- iValue
      nValues
    logical ( KDL ) :: &
      UseDevice

    UseDevice = .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice = UseDeviceOption

    nValues = size ( K )
    
    if ( UseDevice ) then
      
      !$OMP  OMP_TARGET_DIRECTIVE parallel do &
      !$OMP  schedule ( OMP_SCHEDULE_TARGET )
      do iV = 1, nValues
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  &
                                  *  (    V_1 ( iV )  *  GradPhi_1 ( iV )  &
                                       +  V_2 ( iV )  *  GradPhi_2 ( iV )  &
                                       +  V_3 ( iV )  *  GradPhi_3 ( iV ) ) &
                                  *  dt
        S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  &
                                  *  (    V_1 ( iV )  *  GradPhi_1 ( iV )  &
                                       +  V_2 ( iV )  *  GradPhi_2 ( iV )  &
                                       +  V_3 ( iV )  *  GradPhi_3 ( iV ) ) &
                                  *  Weight_RK
      end do
      !$OMP  end OMP_TARGET_DIRECTIVE parallel do
    
    else

      !$OMP  parallel do &
      !$OMP  schedule ( OMP_SCHEDULE_HOST )
      do iV = 1, nValues
        if ( .not. IsProperCell ( iV ) ) &
          cycle
        K ( iV )  =  K ( iV )  -  M ( iV )  *  N ( iV )  &
                                  *  (    V_1 ( iV )  *  GradPhi_1 ( iV )  &
                                       +  V_2 ( iV )  *  GradPhi_2 ( iV )  &
                                       +  V_3 ( iV )  *  GradPhi_3 ( iV ) ) &
                                  *  dt
        S ( iV )  =  S ( iV )  -  M ( iV )  *  N ( iV )  &
                                  *  (    V_1 ( iV )  *  GradPhi_1 ( iV )  &
                                       +  V_2 ( iV )  *  GradPhi_2 ( iV )  &
                                       +  V_3 ( iV )  *  GradPhi_3 ( iV ) ) &
                                  *  Weight_RK
      end do
      !$OMP  end parallel do
      
    end if

  end subroutine ApplyGravityEnergy


end module ApplyGravity_F__Command
