module ConservationLawEvolution_Template

  use Basics
  use DistributedMesh_Form
  use ConservedFields_Template
  use ConservationLawStep_Form

  implicit none
  private

  type, public, abstract :: ConservationLawEvolutionTemplate
    integer ( KDI ) :: &
      iCycle, &
      nRampCycles, &
      nWrite
    real ( KDR ) :: &
      CourantFactor, &
      StartTime, &
      FinishTime, &
      WriteTimeInterval, &
      Time, &
      WriteTime, &
      TimeStep
    type ( MeasuredValueForm ) :: &
      TimeUnit
    character ( LDF ) :: &
      Type = ''
    type ( DistributedMeshForm ) :: &
      DistributedMesh
    class ( ConservedFieldsTemplate ), allocatable :: &
      ConservedFields
    type ( ConservationLawStepForm ) :: &
      ConservationLawStep
  contains
    procedure, public, pass :: &
      Initialize_CLE
    generic, public :: &
      Initialize => Initialize_CLE
    procedure, public, pass :: &
      Evolve
  end type ConservationLawEvolutionTemplate

    private :: &
      ComputeTimeStep

      private :: &
        ComputeTimeStepKernel

contains


  subroutine Initialize_CLE ( CLE, C, BoundaryConditionOption )

    class ( ConservationLawEvolutionTemplate ), intent ( inout ) :: &
      CLE
    type ( CommunicatorForm ), intent ( in ) :: &
      C
    character ( * ), intent ( in ), optional :: &
      BoundaryConditionOption

    call Show ( 'Initializing ' // trim ( CLE % Type ), CONSOLE % INFO_2 )

    associate ( DM => CLE % DistributedMesh )
    call DM % Initialize ( C, BoundaryConditionOption )

    CLE % iCycle = 0
    CLE % nRampCycles = 100
    call PROGRAM_HEADER % GetParameter ( CLE % nRampCycles, 'nRampCycles' )

    select case ( DM % nDimensions )
    case ( 1 )
      CLE % CourantFactor = 0.7
    case ( 2 ) 
      CLE % CourantFactor = 0.4
    case ( 3 )
      CLE % CourantFactor = 0.25
    end select
    call PROGRAM_HEADER % GetParameter &
           ( CLE % CourantFactor, 'CourantFactor' )

    CLE % StartTime  = 0.0_KDR
    CLE % FinishTime = 1.0_KDR
    CLE % TimeUnit   = UNIT % IDENTITY
    call PROGRAM_HEADER % GetParameter &
           ( CLE % StartTime, 'StartTime', InputUnitOption = CLE % TimeUnit )    
    call PROGRAM_HEADER % GetParameter &
           ( CLE % FinishTime, 'FinishTime', InputUnitOption = CLE % TimeUnit )

    CLE % nWrite = 100
    call PROGRAM_HEADER % GetParameter ( CLE % nWrite, 'nWrite' )

    CLE % WriteTimeInterval = ( CLE % FinishTime - CLE % StartTime ) / CLE % nWrite
    call PROGRAM_HEADER % GetParameter &
           ( CLE % WriteTimeInterval, 'WriteTimeInterval' )

    !-- Extensions are responsible for initializing CLE % ConservedFields

    !-- CLE % ConservationLawStep initialized below in CLE % Evolve

    end associate !-- DM

  end subroutine Initialize_CLE


  subroutine Evolve ( CLE )

    class ( ConservationLawEvolutionTemplate ), intent ( inout ) :: &
      CLE
      
    integer ( KDI ) :: &
      iTimerComputation
      
    call PROGRAM_HEADER % AddTimer ( 'Computational', iTimerComputation )
    
    associate &
      ( DM  => CLE % DistributedMesh, &
        CLS => CLE % ConservationLawStep, &
        T => PROGRAM_HEADER % Timer ( iTimerComputation ) )

    call CLS % Initialize ( CLE % ConservedFields )    

    CLE % Time = CLE % StartTime

    call DM % Write &
           ( TimeOption = CLE % Time / CLE % TimeUnit, &
             CycleNumberOption = CLE % iCycle )
    CLE % WriteTime &
      = min ( CLE % Time + CLE % WriteTimeInterval, CLE % FinishTime )

    call Show ( 'Evolving a Fluid', CONSOLE % INFO_2 )
    
    call T % Start ( ) 

    do while ( CLE % Time < CLE % FinishTime )

      call Show ( 'Solving Conservation Equations', CONSOLE % INFO_3 )

      call ComputeTimeStep ( CLE )
      if ( CLE % Time + CLE % TimeStep > CLE % WriteTime ) &
        CLE % TimeStep = CLE % WriteTime - CLE % Time
      call Show ( CLE % TimeStep, CLE % TimeUnit, 'TimeStep', &
                  CONSOLE % INFO_3 )
      
      call CLS % Solve ( CLE % TimeStep )

      CLE % iCycle = CLE % iCycle + 1
      CLE % Time = CLE % Time + CLE % TimeStep
      call Show ( CLE % iCycle, 'iCycle', CONSOLE % INFO_3 )
      call Show ( CLE % Time, CLE % TimeUnit, 'Time', CONSOLE % INFO_3 )

      if ( CLE % Time >= CLE % WriteTime ) then
        
        call T % Stop ( )

        call DM % Write &
               ( TimeOption = CLE % Time / CLE % TimeUnit, &
                 CycleNumberOption = CLE % iCycle )
        CLE % WriteTime &
          = min ( CLE % Time + CLE % WriteTimeInterval, CLE % FinishTime )
        
        call T % Start ( )
        
      end if

    end do
    
    call T % Stop ( )
    
    end associate !-- DM, etc.

  end subroutine Evolve


  subroutine ComputeTimeStep ( CLE )

    class ( ConservationLawEvolutionTemplate ), intent ( inout ) :: &
      CLE

    real ( KDR ), dimension ( :, :, : ), pointer :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ) :: &
      RampFactor
    type ( CollectiveOperation_R_Form ) :: &
      CO

    associate &
      ( DM => CLE % DistributedMesh, &
        CF => CLE % ConservedFields )

    RampFactor &
      = min ( real ( CLE % iCycle + 1, KDR ) / CLE % nRampCycles, 1.0_KDR )

    !-- Only proper cells!

    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_PLUS ( 1 ) ), FEP_1 )
    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_PLUS ( 2 ) ), FEP_2 )
    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_PLUS ( 3 ) ), FEP_3 )
    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_MINUS ( 1 ) ), FEM_1 )
    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_MINUS ( 2 ) ), FEM_2 )
    call DM % SetVariablePointer &
           ( CF % Value ( :, CF % FAST_EIGENSPEED_MINUS ( 3 ) ), FEM_3 )

    call CO % Initialize &
           ( PROGRAM_HEADER % Communicator, &
             nOutgoing = [ 1 ], nIncoming = [ 1 ] )
    associate ( nCPB => DM % nCellsPerBrick )
    call ComputeTimeStepKernel &
           ( FEP_1 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             FEP_2 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             FEP_3 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             FEM_1 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             FEM_2 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             FEM_3 ( 1 : nCPB ( 1 ), 1 : nCPB ( 2 ), 1 : nCPB ( 3 ) ), &
             DM % CellWidth, CO % Outgoing % Value ( 1 ) )
    end associate !-- nCPB
    call CO % Reduce ( REDUCTION % MIN )

    CLE % TimeStep &
      = RampFactor * CLE % CourantFactor * CO % Incoming % Value ( 1 )

    end associate !-- DM, etc.

  end subroutine ComputeTimeStep


  subroutine ComputeTimeStepKernel &
               ( FEP_1, FEP_2, FEP_3, FEM_1, FEM_2, FEM_3, &
                 CellWidth, TimeStepLocal )

    real ( KDR ), dimension ( :, :, : ), intent ( in ) :: &
      FEP_1, FEP_2, FEP_3, &
      FEM_1, FEM_2, FEM_3
    real ( KDR ), dimension ( : ), intent ( in ) :: &
      CellWidth
    real ( KDR ), intent ( out ) :: &
      TimeStepLocal

    TimeStepLocal &
      = minval ( CellWidth ) &
        / maxval ( max ( FEP_1, FEP_2, FEP_3, -FEM_1, -FEM_2, -FEM_3 ) )

  end subroutine ComputeTimeStepKernel


end module ConservationLawEvolution_Template
