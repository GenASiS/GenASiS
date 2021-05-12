#include "Preprocessor"

program Step_RK_CSA__Form_Test

  !-- Slope_DivergenceFiniteVolume_Atlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Steps

  implicit none

  type ( GridImageStreamForm ), allocatable :: &
    GIS
  type ( Atlas_SCG_Form ), allocatable :: &
    A
  type ( Stream_A_Form ), allocatable :: &
    SA
  type ( Geometry_F_A_Form ), allocatable :: &
    GA
  type ( CurrentSet_A_Form ), allocatable :: &
    CSA
  type ( Step_RK_CSA_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Step_RK_CSA__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator, &
           PeriodicOption = [ .true., .true., .true. ] )

  allocate ( SA )
  call SA % Initialize ( A, GIS )

  allocate ( GA )
  call GA % Initialize ( A )

  allocate ( CSA )
  call CSA % Initialize( GA )
  call CSA % SetStream ( SA )

  allocate ( S )
  call S % Initialize ( CSA )

  call   A % Show ( )
  call CSA % Show ( )
  call   S % Show ( )
  call  SA % Show ( )

  call SetWave ( CSA, GA )
  call TestStep ( S )

  deallocate ( S )
  deallocate ( CSA )
  deallocate ( GA )
  deallocate ( SA )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CSA, GA )

    class ( CurrentSet_A_Form ), intent ( inout ) :: &
      CSA
    class ( Geometry_F_A_Form ), intent ( in ), target :: &
      GA

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
    class is ( CurrentSet_C_Form )

    select type ( GC  =>  GA % FieldSet_C ( 1 ) % Element )
    class is ( Geometry_F_C_Form )

    select type ( C  =>  CSC % Chart )
    class is ( Chart_GS_Form )

    nWavelengths  =  0
    nWavelengths ( 1 : C % nDimensions )  =  1
    call PROGRAM_HEADER % GetParameter ( nWavelengths, 'nWavelengths' )

    Offset     =  2.0_KDR
    Amplitude  =  1.0_KDR
    Speed      =  1.0_KDR
    call PROGRAM_HEADER % GetParameter ( Offset, 'Offset' )
    call PROGRAM_HEADER % GetParameter ( Amplitude, 'Amplitude' )
    call PROGRAM_HEADER % GetParameter ( Speed, 'Speed' )

    associate ( BoxSize  =>  C % MaxCoordinate  -  C % MinCoordinate )
    where ( BoxSize  >  0.0_KDR )
      Wavenumber  =  nWavelengths / BoxSize
    elsewhere
      Wavenumber  =  0.0_KDR
    end where
    end associate !-- BoxSize

    associate &
      (     X  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_1 ), &
            Y  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_2 ), &
            Z  =>   GC % Storage_FSC % Storage &
                       % Value ( :, GC % CENTER_U_3 ), &
          Rho  =>  CSC % Storage_FSC % Storage &
                       % Value ( :, CSC % DENSITY_DEFAULT ), &
            V  =>  CSC % VelocityDefault_U, &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V ( 1 )  =  Speed  *  K ( 1 )  /  Abs_K
    V ( 2 )  =  Speed  *  K ( 2 )  /  Abs_K
    V ( 3 )  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end select !-- C
    end select !-- GC
    end select !-- CSC

  end subroutine SetWave


  subroutine TestStep ( S )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    real ( KDR ) :: &
      CourantFactor, &
      TimeStep

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

    CourantFactor  =  0.7_KDR
         TimeStep  =  huge ( 1.0_KDR )
    call ComputeTimeStep ( S, TimeStep, CourantFactor )
    call Show ( TimeStep, 'TimeStep' )

    call S % Compute ( T = 0.0_KDR, dT = TimeStep )

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write ( )
    call GIS % Close ( )

  end subroutine TestStep


  subroutine ComputeTimeStep ( S, TimeStep, CourantFactor )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( inout ) :: &
      TimeStep
    real ( KDR ), intent ( in ) :: &
      CourantFactor

    ! class ( GeometryFlatForm ), pointer :: &
    !   G
    ! class ( CurrentTemplate ), pointer :: &
    !   C

    ! select type ( PS => I % PositionSpace )
    ! class is ( Atlas_SC_Form )

    ! G => CSL % Geometry ( )
    ! C => CA % Current ( )
    associate ( EA  =>  S % RiemannSolver_A % Eigenspeeds_A )

    select type ( EC  =>  EA % FieldSet_C ( 1 ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( GC  =>  CSA % Geometry_A % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )
    select type ( CGS  =>  CSA % Atlas % Chart ( 1 ) % Element )
      class is ( Chart_GS_Form )

    call EC % Compute ( iD = 1 )

    associate &
      ( EV  =>  EC % Storage_FSC % Storage % Value, &
        GV  =>  GC % Storage_FSC % Storage % Value, &
        DeviceMemory  =>  EC % Storage_FSC % DeviceMemory )

    call ComputeTimeStepKernel &
           ( TimeStep, CGS % ProperCell, &
             EV ( :, EC % EIGENSPEED_FAST_PLUS_U ), &
             EV ( :, EC % EIGENSPEED_FAST_PLUS_U ), &
             EV ( :, EC % EIGENSPEED_FAST_PLUS_U ), &
             EV ( :, EC % EIGENSPEED_FAST_MINUS_U ), &
             EV ( :, EC % EIGENSPEED_FAST_MINUS_U ), &
             EV ( :, EC % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, GC % WIDTH_U_1 ), &
             GV ( :, GC % WIDTH_U_2 ), &
             GV ( :, GC % WIDTH_U_3 ), &
             CGS % nDimensions, &
             UseDeviceOption = DeviceMemory )

    end associate !-- EV, etc.
    end select !-- CGS
    end select !-- GC
    end select !-- EC
    end associate !-- EA

    TimeStep  =  CourantFactor  *  TimeStep
    
  end subroutine ComputeTimeStep


  subroutine ComputeTimeStepKernel &
               ( TimeStep, ProperCell, &
                 FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, dX_1, dX_2, dX_3, &
                 nDimensions, UseDeviceOption )

    real ( KDR ), intent ( inout ) :: &
      TimeStep
    logical ( KDL ), dimension ( : ), intent ( in ) :: &
      ProperCell
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3, &
      dX_1, dX_2, dX_3
    integer ( KDI ), intent ( in ) :: &
      nDimensions
    logical ( KDL ), intent ( in ), optional :: &
      UseDeviceOption

    integer ( KDI ) :: &
      iV, &
      nV
    real ( KDR ) :: &
      TimeStepInverse
    logical ( KDL ) :: &
      UseDevice
      
    UseDevice  =  .false.
    if ( present ( UseDeviceOption ) ) &
      UseDevice  =  UseDeviceOption
      
    nV  =  size ( FEP_1 )
    
    TimeStepInverse  =  - huge ( 0.0_KDR )

    select case ( nDimensions )
    case ( 1 )

      !TimeStepInverse &
      !  = maxval ( max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 2 )

      !TimeStepInverse &
      !  = maxval (   max ( FEP_1, -FEM_1 ) / ( dXL_1 + dXR_1 ) &
      !             + max ( FEP_2, -FEM_2 ) / ( Crsn_2 * ( dXL_2 + dXR_2 ) ), &
      !             mask = ProperCell )
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / dX_2 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / dX_2 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    case ( 3 )
      ! TimeStepInverse &
      !   = maxval (   max ( FEP_1, -FEM_1 ) / dX_1 &
      !              + max ( FEP_2, -FEM_2 ) / dX_2 &
      !              + max ( FEP_3, -FEM_3 ) / dX_3, &
      !              mask = ProperCell )
      
      if ( UseDevice ) then
        !$OMP  OMP_TARGET_DIRECTIVE parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_TARGET ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / dX_2 ( iV ) &
                      + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                        / dX_3 ( iV ) )
        end do
        !$OMP  end OMP_TARGET_DIRECTIVE parallel do
      else
        !$OMP  parallel do &
        !$OMP& schedule ( OMP_SCHEDULE_HOST ) private ( iV ) &
        !$OMP& reduction ( max : TimeStepInverse )
        do iV = 1, nV
          if ( ProperCell ( iV ) ) &
            TimeStepInverse &
              = max ( TimeStepInverse, &
                        max ( FEP_1 ( iV ), -FEM_1 ( iV ) ) &
                        / dX_1 ( iV ) &
                      + max ( FEP_2 ( iV ), -FEM_2 ( iV ) ) &
                        / dX_2 ( iV ) &
                      + max ( FEP_3 ( iV ), -FEM_3 ( iV ) ) &
                        / dX_3 ( iV ) )
        end do
        !$OMP  end parallel do
      end if
      
    end select !-- nDimensions

    TimeStepInverse = max ( tiny ( 0.0_KDR ), TimeStepInverse )
    TimeStep = min ( TimeStep, 1.0_KDR / TimeStepInverse )

  end subroutine ComputeTimeStepKernel


end program Step_RK_CSA__Form_Test
