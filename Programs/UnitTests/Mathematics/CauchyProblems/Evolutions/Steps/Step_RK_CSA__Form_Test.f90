#include "Preprocessor"

program Step_RK_CSA__Form_Test

  !-- Slope_DivergenceFiniteVolume_Atlas__Form_Test

  use Basics
  use Manifolds
  use Fields
  use Steps

  implicit none

  integer ( KDI ) :: &
    iD
  character ( 1 ) :: &
    Dimension
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
  class ( FieldSet_A_Element ), dimension ( : ), allocatable :: &
    EA
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

  allocate ( EA ( 3 ) )
  do iD  =  1, 3
    allocate ( Eigenspeeds_F_A_Form :: EA ( iD ) % Element )
    select type ( EA_iD  =>  EA ( iD ) % Element )
      class is ( Eigenspeeds_F_A_Form )
    write ( Dimension, fmt = '(i1.1)' ) iD
    call EA_iD % Initialize &
           ( CSA, &
             NameOption = 'E_' // Dimension // '_' // trim ( CSA % Name ) ) 
    end select !-- EA_iD
  end do !-- iD

  allocate ( S )
  call S % Initialize ( CSA )
  call S % SetStream ( SA, StagesOption = .true. )

  call   A % Show ( )
  call CSA % Show ( )
  call   S % Show ( )
  call  SA % Show ( )

  call SetWave ( CSA, GA )
  call TestStep ( S, EA )

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
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    call CSC % SetVelocityDefault ( Wavenumber, Speed )
    
    end associate !-- Rho, etc.
    end select !-- C
    end select !-- GC
    end select !-- CSC

  end subroutine SetWave


  subroutine TestStep ( S, EA )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    class ( FieldSet_A_Element ), dimension ( : ), intent ( inout ) :: &
      EA

    integer ( KDI ) :: &
      iC, &  !-- iCycle
      nCycles
    real ( KDR ) :: &
      CourantFactor, &
      Time, &
      TimeStep

    Time  =  0.0_KDR

    call GIS % Open ( GIS % ACCESS_CREATE )
    call SA % Write &
           ( TimeOption = Time  *  UNIT % IDENTITY, &
             CycleNumberOption = 0 )
    call GIS % Close ( )

    CourantFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter ( CourantFactor, 'CourantFactor' )

    nCycles  =  1
    call PROGRAM_HEADER % GetParameter ( nCycles, 'nCycles' )

    call Show ( CourantFactor, 'CourantFactor' )
    call Show ( nCycles, 'nCycles' )

    do iC  =  1,  nCycles

      call ComputeTimeStep ( EA, S, CourantFactor, TimeStep )
      call Show ( iC, 'iCycle' )
      call Show ( TimeStep, 'TimeStep' )

      call S % Compute ( T = 0.0_KDR, dT = TimeStep )

      Time  =  Time + TimeStep

      call GIS % Open ( GIS % ACCESS_CREATE )
      call SA % Write &
             ( TimeOption = Time  *  UNIT % IDENTITY, &
               CycleNumberOption = 0 )
      call GIS % Close ( )

    end do !-- iS

  end subroutine TestStep


  subroutine ComputeTimeStep ( EA, S, CourantFactor, TimeStep )

    class ( FieldSet_A_Element ), dimension ( : ), intent ( inout ) :: &
      EA
    class ( Step_RK_CSA_Form ), intent ( in ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      CourantFactor
    real ( KDR ), intent ( out ) :: &
      TimeStep

    type ( CollectiveOperation_R_Form ) :: &
      CO

    call ComputeTimeStepLocal ( EA, S, CourantFactor, TimeStep )

    call CO % Initialize &
           ( PROGRAM_HEADER % Communicator, &
             nOutgoing = [ 1 ], &
             nIncoming = [ 1 ] )
    CO % Outgoing % Value  =  TimeStep

    call CO % Reduce ( REDUCTION % MIN )

    TimeStep  =  CO % Incoming % Value ( 1 )

  end subroutine ComputeTimeStep


  subroutine ComputeTimeStepLocal ( EA, S, CourantFactor, TimeStep )

    class ( FieldSet_A_Element ), dimension ( : ), intent ( inout ) :: &
      EA
    class ( Step_RK_CSA_Form ), intent ( in ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      CourantFactor
    real ( KDR ), intent ( out ) :: &
      TimeStep

    integer ( KDI ) :: &
      iD

    TimeStep  =  huge ( 1.0_KDR )
    
    select type ( EC_1  =>  EA ( 1 ) % Element % FieldSet_C ( 1 ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( EC_2  =>  EA ( 2 ) % Element % FieldSet_C ( 1 ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( EC_3  =>  EA ( 3 ) % Element % FieldSet_C ( 1 ) % Element )
      class is ( Eigenspeeds_F_C_Form )
    select type ( GC  =>  CSA % Geometry_A % FieldSet_C ( 1 ) % Element )
      class is ( Geometry_F_C_Form )
    select type ( CGS  =>  CSA % Atlas % Chart ( 1 ) % Element )
      class is ( Chart_GS_Form )

    call EC_1 % Compute ( iD = 1 )
    call EC_2 % Compute ( iD = 2 )
    call EC_3 % Compute ( iD = 3 )

    associate &
      ( EV_1  =>  EC_1 % Storage_FSC % Storage % Value, &
        EV_2  =>  EC_2 % Storage_FSC % Storage % Value, &
        EV_3  =>  EC_3 % Storage_FSC % Storage % Value, &
        GV    =>  GC   % Storage_FSC % Storage % Value, &
        DeviceMemory  =>  GC % Storage_FSC % DeviceMemory )

    call ComputeTimeStepKernel &
           ( TimeStep, CGS % ProperCell, &
             EV_1 ( :, EC_1 % EIGENSPEED_FAST_PLUS_U ), &
             EV_2 ( :, EC_2 % EIGENSPEED_FAST_PLUS_U ), &
             EV_3 ( :, EC_3 % EIGENSPEED_FAST_PLUS_U ), &
             EV_1 ( :, EC_1 % EIGENSPEED_FAST_MINUS_U ), &
             EV_2 ( :, EC_2 % EIGENSPEED_FAST_MINUS_U ), &
             EV_3 ( :, EC_3 % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, GC % WIDTH_U_1 ), &
             GV ( :, GC % WIDTH_U_2 ), &
             GV ( :, GC % WIDTH_U_3 ), &
             CGS % nDimensions, &
             UseDeviceOption = DeviceMemory )

    end associate !-- EV, etc.
    end select !-- CGS
    end select !-- GC
    end select !-- EC_3
    end select !-- EC_2
    end select !-- EC_1

    TimeStep  =  CourantFactor  *  TimeStep
    
  end subroutine ComputeTimeStepLocal


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
