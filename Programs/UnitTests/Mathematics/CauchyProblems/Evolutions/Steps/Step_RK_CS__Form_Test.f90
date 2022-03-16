#include "Preprocessor"

program Step_RK_CS__Form_Test

  !-- Step_RungeKutta_CurrentSet__Form_Test

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
  type ( StreamForm ), allocatable :: &
    Sm
  type ( Geometry_F_Form ), allocatable :: &
    G
  type ( CurrentSetForm ), allocatable :: &
    CS
  class ( EigenspeedSet_F_Form ), dimension ( : ), allocatable :: &
    ES
  type ( Step_RK_CS_Form ), allocatable :: &
    S

  allocate ( PROGRAM_HEADER )
  call PROGRAM_HEADER % Initialize &
         ( 'Step_RK_CS__Form_Test', DimensionalityOption = '2D' )

  allocate ( GIS )
  call GIS % Initialize &
         ( PROGRAM_HEADER % Name, &
           CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( A )
  call A % Initialize &
         ( CommunicatorOption = PROGRAM_HEADER % Communicator )

  allocate ( Sm )
  call Sm % Initialize ( A, GIS )

  allocate ( G )
  call G % Initialize ( A )
  call G % SetStream ( Sm )

  allocate ( CS )
  call CS % Initialize( G )
  call CS % SetStream ( Sm )
  do iD  =  1, 3
    call CS % SetBoundaryConditionsFace &
           ( [ 'PERIODIC', 'PERIODIC' ], iC = 1, iD = iD )
  end do !-- iD

  allocate ( ES ( 3 ) )
  do iD  =  1, 3
    write ( Dimension, fmt = '(i1.1)' ) iD
    call ES ( iD ) % Initialize &
           ( CS, CS, SuffixOption = 'Egnspd_' // Dimension )
  end do !-- iD

  allocate ( S )
  call S % Initialize ( CS )
  call S % SetStream ( Sm )

  call  A % Show ( )
  call  G % Show ( )
  call CS % Show ( )
  call S % Show ( )
  call Sm % Show ( )

  call SetWave ( CS, G )
  call TestStep ( S, ES )

  deallocate ( S )
  deallocate ( ES )
  deallocate ( CS )
  deallocate ( G )
  deallocate ( Sm )
  deallocate ( A )
  deallocate ( GIS )
  deallocate ( PROGRAM_HEADER )


contains


  subroutine SetWave ( CS, G )

    class ( CurrentSetForm ), intent ( inout ) :: &
      CS
    class ( Geometry_F_Form ), intent ( in ), target :: &
      G

    integer ( KDI ), dimension ( 3 ) :: &
      nWavelengths
    real ( KDR ) :: &
      Offset, &
      Amplitude, &
      Speed
    real ( KDR ), dimension ( 3 ) :: &
      Wavenumber

    select type ( A  =>  CS % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS )

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
      (  GV  =>   G % Storage_GS % Value, &
        CSV  =>  CS % Storage_GS % Value )
    associate &
      (     X  =>   GV ( :,  G % CENTER_U_1 ), &
            Y  =>   GV ( :,  G % CENTER_U_2 ), &
            Z  =>   GV ( :,  G % CENTER_U_3 ), &
          Rho  =>  CSV ( :, CS % DENSITY_CS ), &
          V_1  =>  CSV ( :, CS % VELOCITY_CS_U_1 ), &
          V_2  =>  CSV ( :, CS % VELOCITY_CS_U_2 ), &
          V_3  =>  CSV ( :, CS % VELOCITY_CS_U_3 ), &
            K  =>  Wavenumber, &
        Abs_K  =>  sqrt ( dot_product ( Wavenumber, Wavenumber ) ), &
        TwoPi  =>  2.0_KDR  *  CONSTANT % PI )

    Rho  =  Offset  &
            +  Amplitude  &
               *  sin ( TwoPi * (    K ( 1 ) * X  &
                                  +  K ( 2 ) * Y  &
                                  +  K ( 3 ) * Z  ) )

    V_1  =  Speed  *  K ( 1 )  /  Abs_K
    V_2  =  Speed  *  K ( 2 )  /  Abs_K
    V_3  =  Speed  *  K ( 3 )  /  Abs_K
    
    end associate !-- Rho, etc.
    end associate !-- GV, etc.
    end associate !-- C
    end select !-- A

  end subroutine SetWave


  subroutine TestStep ( S, ES )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( EigenspeedSet_F_Form ), dimension ( : ), intent ( inout ) :: &
      ES

    integer ( KDI ) :: &
      iC, &  !-- iCycle
      nCycles
    real ( KDR ) :: &
      CourantFactor, &
      Time, &
      TimeStep
    type ( TimerForm ), pointer :: &
      T

    Time  =  0.0_KDR

    call GIS % Open ( GIS % ACCESS_CREATE )
    call Sm % Write &
           ( TimeOption = Time  *  UNIT % IDENTITY, &
             CycleNumberOption = 0 )
    call GIS % Close ( )

    CourantFactor  =  0.7_KDR
    call PROGRAM_HEADER % GetParameter ( CourantFactor, 'CourantFactor' )

    nCycles  =  1
    call PROGRAM_HEADER % GetParameter ( nCycles, 'nCycles' )

    call Show ( 'Step computation' )
    call Show ( CourantFactor, 'CourantFactor' )
    call Show ( nCycles, 'nCycles' )

    do iC  =  1,  nCycles

      call ComputeTimeStep ( ES, S, CourantFactor, TimeStep )
      call Show ( iC, 'iCycle' )
      call Show ( Time, 'Time' )
      call Show ( TimeStep, 'TimeStep' )

      T  =>  S % Timer ( Level = 1 )
      call T % Start ( )
      call S % Compute ( T = Time, dT = TimeStep, T_Option = T )
      call T % Stop ( )

      Time  =  Time + TimeStep

      T  =>  Sm % TimerWrite ( Level = 1 )
      call T % Start ( )
      call GIS % Open ( GIS % ACCESS_CREATE )
      call Sm % Write &
             ( TimeOption = Time  *  UNIT % IDENTITY, &
               CycleNumberOption = iC )
      call GIS % Close ( )
      call T % Stop ( )

    end do !-- iS

  end subroutine TestStep


  subroutine ComputeTimeStep ( ES, S, CourantFactor, TimeStep )

    class ( EigenspeedSet_F_Form ), dimension ( : ), intent ( inout ) :: &
      ES
    class ( Step_RK_CS_Form ), intent ( in ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      CourantFactor
    real ( KDR ), intent ( out ) :: &
      TimeStep

    type ( CollectiveOperation_R_Form ) :: &
      CO

    call ComputeTimeStepLocal ( ES, S, CourantFactor, TimeStep )

    call CO % Initialize &
           ( PROGRAM_HEADER % Communicator, &
             nOutgoing = [ 1 ], &
             nIncoming = [ 1 ] )
    CO % Outgoing % Value  =  TimeStep

    call CO % Reduce ( REDUCTION % MIN )

    TimeStep  =  CO % Incoming % Value ( 1 )

  end subroutine ComputeTimeStep


  subroutine ComputeTimeStepLocal ( ES, S, CourantFactor, TimeStep )

    class ( EigenspeedSet_F_Form ), dimension ( : ), intent ( inout ) :: &
      ES
    class ( Step_RK_CS_Form ), intent ( in ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      CourantFactor
    real ( KDR ), intent ( out ) :: &
      TimeStep

    integer ( KDI ) :: &
      iD

    TimeStep  =  huge ( 1.0_KDR )
    
    call ES ( 1 ) % Compute ( iC = 1, iD = 1 )
    call ES ( 2 ) % Compute ( iC = 1, iD = 2 )
    call ES ( 3 ) % Compute ( iC = 1, iD = 3 )

    select type ( A  =>  ES ( 1 ) % Atlas )
      class is ( Atlas_SCG_Form )
    associate &
      ( C  =>  A % Chart_GS, &
        G  =>  ES ( 1 ) % CurrentSet % Geometry )

    associate &
      ( EV_1  =>  ES ( 1 ) % Storage_GS % Value, &
        EV_2  =>  ES ( 2 ) % Storage_GS % Value, &
        EV_3  =>  ES ( 3 ) % Storage_GS % Value, &
        GV    =>   G       % Storage_GS % Value )

    call ComputeTimeStepKernel &
           ( TimeStep, C % ProperCell, &
             EV_1 ( :, ES ( 1 ) % EIGENSPEED_FAST_PLUS_U ), &
             EV_2 ( :, ES ( 2 ) % EIGENSPEED_FAST_PLUS_U ), &
             EV_3 ( :, ES ( 3 ) % EIGENSPEED_FAST_PLUS_U ), &
             EV_1 ( :, ES ( 1 ) % EIGENSPEED_FAST_MINUS_U ), &
             EV_2 ( :, ES ( 2 ) % EIGENSPEED_FAST_MINUS_U ), &
             EV_3 ( :, ES ( 3 ) % EIGENSPEED_FAST_MINUS_U ), &
             GV ( :, G % WIDTH_U_1 ), &
             GV ( :, G % WIDTH_U_2 ), &
             GV ( :, G % WIDTH_U_3 ), &
             C % nDimensions, &
             UseDeviceOption = G % DeviceMemory )

    end associate !-- EV_1, etc.
    end associate !-- C, etc.
    end select !-- A

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


end program Step_RK_CS__Form_Test
