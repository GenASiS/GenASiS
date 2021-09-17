module Step_RK_CS__Form

  !-- Step_RungeKutta_CurrentSet_Form

  use Basics
  use Algebra
  use Fields
  use Slopes
  use Step_RK_H__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CS_Form
    type ( FieldSetForm ), allocatable :: &
      Balanced, &
      Intermediate, &
      Solution
    type ( FieldSetElement ), dimension ( : ), allocatable :: &
      SolutionStage
    class ( CurrentSetForm ), pointer :: &
      CurrentSet
    class ( RiemannSolver_HLL_Form ), allocatable :: &
      RiemannSolver
  contains
    procedure, private, pass :: &
      Initialize_CS
    generic, public :: &
      Initialize => Initialize_CS
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_S
    final :: &
      Finalize
    procedure, private, pass :: &
      LoadSolution
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, private, pass :: &
      StoreSolution
  end type Step_RK_CS_Form

    private :: &
      SetSlope_CS

contains


  subroutine Initialize_CS ( S, CS, NameOption, A_Option, B_Option, C_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    real ( KDR ), dimension ( 2 : , : ), intent ( in ), optional :: &
      A_Option
    real ( KDR ), dimension ( : ), intent ( in ), optional :: &
      B_Option
    real ( KDR ), dimension ( 2 : ), intent ( in ), optional :: &
      C_Option

    integer ( KDI ) :: &
      iS  !-- iStage
    character ( 1 ) :: &
      StageNumber
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CS'

    Name  =  'Step_' // trim ( CS % Name )
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    S % CurrentSet  =>  CS

    !-- Balanced

    allocate ( S % Balanced )
    associate ( CS_B  =>  S % Balanced )
    call CS_B % Initialize &
           ( CS, iaSelected = CS % iaBalanced, &
             NameOption = 'Balanced', &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- CS_B

    !-- Intermediate storage

    allocate ( S % Intermediate )
    associate ( Y_I  =>  S % Intermediate )
    call Y_I % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = 'Intermediate', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- Y_I

    !-- Solution storage

    allocate ( S % Solution )
    associate ( Y  =>  S % Solution )
    call Y % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = 'Solution', &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    end associate !-- Y

    !-- RiemannSolver

    if ( .not. allocated ( S % RiemannSolver ) ) then
      allocate ( S % RiemannSolver )
      associate ( RS  =>  S % RiemannSolver )
      call RS % Initialize ( CS )
      end associate !-- RS
    end if !-- allocated RiemannSolver

    !-- Header

    if ( .not. associated ( S % SetSlope ) ) &
      S % SetSlope  =>  SetSlope_CS

    call S % Initialize_H &
           ( CS % Atlas, &
             NameOption = Name, &
             A_Option = A_Option, &
             B_Option = B_Option, &
             C_Option = C_Option )

  end subroutine Initialize_CS


  subroutine SetStream ( S, Sm, StagesOption )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( StreamForm ), intent ( inout ) :: &
      Sm
    logical ( KDL ), intent ( in ), optional :: &
      StagesOption

    integer ( KDI ) :: &
      iS  !-- iStage
    logical ( KDL ) :: &
      Stages
    character ( 1 ) :: &
      StageNumber

    Stages  =  .false.
    if ( present ( StagesOption ) ) &
      Stages  =  StagesOption
    call PROGRAM_HEADER % GetParameter ( Stages, 'StreamStages' )

    call S % SetStream_H ( Sm, StagesOption = Stages )

    if ( Stages ) then
      associate ( RSA  =>  S % RiemannSolver )
      call RSA % SetStream ( Sm, S % nStages )
      end associate !-- RSA
    end if !-- Stages

  end subroutine SetStream


  subroutine Show_S ( S )

    class ( Step_RK_CS_Form ), intent ( in ) :: &
      S

   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord, &
     TypePiece

    call S % Step_RK_H_Form % Show ( )

    call S % Solution % Show ( )
    call S % Intermediate % Show ( )
    call S % RiemannSolver % Show ( )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % RiemannSolver ) ) &
      deallocate ( S % RiemannSolver )
    if ( allocated ( S % SolutionStage ) ) &
      deallocate ( S % SolutionStage )
    if ( allocated ( S % Solution ) ) &
      deallocate ( S % Solution )
    if ( allocated ( S % Intermediate ) ) &
      deallocate ( S % Intermediate )
    if ( allocated ( S % Balanced ) ) &
      deallocate ( S % Balanced )

    nullify ( S % CurrentSet )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    associate &
      ( CS_B  =>  S % Balanced, &
        Y     =>  S % Solution )

    call CS_B % Copy ( Y )

    !-- For diagnostic I/O
    if ( allocated ( S % SolutionStage ) ) then
      associate ( Y_S  =>  S % SolutionStage ( 1 ) % Element )
      call CS_B % Copy ( Y_S )
      end associate !-- Y_S
    end if

    end associate !-- CS_B, etc.

  end subroutine LoadSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    if ( iS  >  1 ) then
      associate &
        ( Y    =>  S % Solution, &
          Y_I  =>  S % Intermediate )

      call Y % Copy ( Y_I )

      end associate !-- Y, etc.
    end if !-- iS > 1

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

    associate &
      ( Y_I  =>  S % Intermediate, &
        K    =>  S % SlopeStage ( iK ) % Element )

    call Y_I % MultiplyAdd ( K, dT * A )

    end associate !-- Y_I, etc.

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, T, iS, T_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart
    type ( TimerForm ), pointer :: &
      T_C, &
      T_EG

    associate ( K  =>  S % SlopeStage ( iS ) % Element )

    if ( iS  >  1 ) then

      associate ( K_1  =>  S % SlopeStage ( 1 ) % Element )
      call K % CloneTimers ( K_1 )
      end associate !-- K_1

      associate &
        (  Y_I  =>  S % Intermediate, &
          CS_B  =>  S % Balanced, &
          CS    =>  S % CurrentSet )

      call Y_I % Copy ( CS_B )
      call CS % ComputeFromBalanced ( )
      call CS % ApplyBoundaryConditions ( )
    
      !-- For diagnostic I/O
      if ( allocated ( S % SolutionStage ) ) then
        associate ( Y_S  =>  S % SolutionStage ( iS ) % Element )
        call Y_I % Copy ( Y_S )
        end associate !-- Y_S
      end if

      end associate !-- Y_I, etc.
 
    end if !-- iStage > 1

    !-- Compute slope

    if ( present ( T_Option ) ) then
      T_C  =>  K % Timer ( LevelOption = T_Option % Level + 1 )
      call T_C % Start ( )
    else
      T_C   =>  null ( )
    end if
    call K % Compute ( T_Option = T_C, iS_Option = iS )
    if ( associated ( T_C ) ) call T_C % Stop ( )

    !-- Slope ghost exchange

    if ( present ( T_Option ) ) then
      T_EG  =>  K % TimerGhost &
                  ( NameRootOption = K % TimerName, &
                    LevelOption = T_Option % Level + 1 )
      call T_EG % Start ( )
    else
      T_EG  =>  null ( )
    end if
    call K % ExchangeGhostData ( )
    if ( associated ( T_EG ) ) call T_EG % Stop ( )

    end associate !-- K

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, dT, iS )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC  !-- iChart

    associate &
      ( Y  =>  S % Solution, &
        K  =>  S % SlopeStage ( iS ) % Element )

    call Y % MultiplyAdd ( K, dT * B )

    end associate !-- Y, etc.

  end subroutine IncrementSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S

    associate &
      (  Y    =>  S % Solution, &
        CS_B  =>  S % Balanced, &
        CS    =>  S % CurrentSet )

    call  Y % Copy ( CS_B )
    call CS % ComputeFromBalanced ( )
    call CS % ApplyBoundaryConditions ( )
  
    end associate !-- Y, etc.
 
  end subroutine StoreSolution


  subroutine SetSlope_CS ( S, K, iS_Option )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    character ( 1 ) :: &
      StageNumber

    allocate ( Slope_DFV_F_Form :: K )
    select type ( K )
      class is ( Slope_DFV_F_Form )
    select type ( S )
      class is ( Step_RK_CS_Form )

    if ( present ( iS_Option ) ) then
      write ( StageNumber, fmt = '(i1.1)' ) iS_Option
      call K % Initialize ( S % RiemannSolver, SuffixOption = StageNumber )
    else
      call K % Initialize ( S % RiemannSolver, &
                            IgnorabilityOption = S % IGNORABILITY )
    end if

    end select !-- S
    end select !-- K

  end subroutine SetSlope_CS


end module Step_RK_CS__Form
