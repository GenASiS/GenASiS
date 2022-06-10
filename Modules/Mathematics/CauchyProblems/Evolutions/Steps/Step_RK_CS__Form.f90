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
     integer ( KDI ) :: &
      iTimer_Crsn    = 0, &  !-- Coarsen
      iTimer_SC      = 0, &  !-- SolutionCopy
      iTimer_CFB     = 0, &  !-- ComputeFromBalanced
      iTimer_BC      = 0     !-- BoundaryCondition
    type ( FieldSetForm ), allocatable :: &
      Balanced, &
      Intermediate, &
      Solution
    class ( CurrentSetForm ), pointer :: &
      CurrentSet
    class ( Coarsening_C_Form ), pointer :: &
      Coarsening => null ( )
    class ( DivergencePart_CS_Form ), allocatable :: &
      DivergenceTotal
    type ( DivergencePartElement ), dimension ( : ), allocatable :: &
      DivergencePart
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
      SetCoarsening
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
      SetSlope_CS, &
      SetSlopeStage_CS, &
      StoreSolution_CS

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
    logical ( KDL ) :: &
      DivergenceParts
    character ( 1 ) :: &
      StageNumber
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CS'

    Name  =  trim ( CS % Name ) // '_Stp' 
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

    !-- DivergenceTotal or DivergencePart

    if (       .not. allocated ( S % DivergenceTotal ) &
         .and. .not. allocated ( S % DivergencePart ) ) then

      DivergenceParts  =  .false.
      call PROGRAM_HEADER % GetParameter ( DivergenceParts, 'DivergenceParts' )

      if ( DivergenceParts ) then
        allocate ( S % DivergencePart ( 1 ) )
        associate ( DP_1D  =>  S % DivergencePart )
        allocate ( DP_1D ( 1 ) % Element )
        associate ( DP  =>  DP_1D ( 1 ) % Element )
        call DP % Initialize ( CS )
        end associate !-- DP
        end associate !-- DP_1D
      else  !-- DivergenceTotal
        allocate ( S % DivergenceTotal )
        associate ( DT  =>  S % DivergenceTotal )
        call DT % Initialize ( CS )
        end associate !-- DT
      end if  !-- DivergenceParts

    end if !-- allocated DivergenceTotal or DivergencePart

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
    if ( .not. associated ( S % SetSlopeStage ) ) &
      S % SetSlopeStage  =>  SetSlopeStage_CS

    call S % Initialize_H &
           ( CS % Atlas, &
             NameOption = Name, &
             A_Option = A_Option, &
             B_Option = B_Option, &
             C_Option = C_Option )

  end subroutine Initialize_CS


  subroutine SetStream ( S, Sm )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( StreamForm ), intent ( inout ) :: &
      Sm

    call S % SetStream_H ( Sm )

  end subroutine SetStream


  subroutine SetCoarsening ( S, C )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    class ( Coarsening_C_Form ), intent ( in ), target :: &
      C

    S % Coarsening  =>  C

  end subroutine SetCoarsening


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
    if ( allocated ( S % DivergencePart ) ) &
      deallocate ( S % DivergencePart )
    if ( allocated ( S % DivergenceTotal ) ) &
      deallocate ( S % DivergenceTotal )
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


  subroutine ComputeStage ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iC  !-- iChart
    type ( TimerForm ), pointer :: &
      T_SS, &  !-- StoreSolution
      T_CS, &  !-- ComputeSlope
      T_EG, &  !-- ExchangeGhost
      T_C,  &  !-- Coarsen
      T_AS     !-- AccumulateSlope

    associate &
      ( K  =>  S % Slope, &
        K_Stage  =>  S % SlopeStage ( iS ) % Element )

    if ( iS  >  1 ) then
      associate ( Y_I  =>  S % Intermediate )
      if ( present ( T_Option ) ) then
        T_SS  =>  S % TimerStoreSolution ( Level = T_Option % Level )
        !-- Pause and then restart ComputeStage timer T_Option to avoid double 
        !   counting time to be attributed to StoreSolution timer T_SS
        call T_Option % Stop ( )  !-- Pause ComputeStage timer
        call T_SS % Start ( )
        call StoreSolution_CS ( S, Y_I, T_Option = T_SS )
        call T_SS % Stop ( )
        call T_Option % Start ( )  !-- Restart ComputeStage timer
      else
        call StoreSolution_CS ( S, Y_I )
      end if
      end associate !-- Y_I, etc.
    end if !-- iStage > 1

    !-- Compute slope

    if ( present ( T_Option ) ) then
      T_CS  =>  K % Timer ( Level = T_Option % Level + 1 )
    else
      T_CS   =>  null ( )
    end if
    if ( associated ( T_CS ) ) call T_CS % Start ( )
    call K % Compute ( T_Option = T_CS )
    if ( associated ( T_CS ) ) call T_CS % Stop ( )

    !-- Coarsening

    if ( associated ( S % Coarsening ) ) then
      if ( present ( T_Option ) ) then
        T_C  =>  PROGRAM_HEADER % Timer &
                   ( Handle = S % iTimer_Crsn, &
                     Name = trim ( K % Name ) // '_Crsn', &
                     Level = T_Option % Level + 1 )
      else
        T_C  =>  null ( )
      end if
      if ( associated ( T_C ) ) call T_C % Start ( )
      call S % Coarsening % Compute ( K )
      if ( associated ( T_C ) ) call T_C % Stop ( )
    end if

    !-- Slope ghost exchange

    if ( present ( T_Option ) ) then
      T_EG  =>  K % TimerGhost ( Level = T_Option % Level + 1 )
    else
      T_EG  =>  null ( )
    end if
    if ( associated ( T_EG ) ) call T_EG % Start ( )
    call K % ExchangeGhostData ( )
    if ( associated ( T_EG ) ) call T_EG % Stop ( )

    !-- Accumulations

    if ( present ( T_Option ) ) then
      T_AS  =>  S % TimerAccumulateSlope ( Level = T_Option % Level + 1 )
    else
      T_AS  =>  null ( )
    end if
    if ( associated ( T_AS ) ) call T_AS % Start ( )

    call K % Copy ( K_Stage )

    call S % AccumulateSlope ( iS ) 

    associate ( CS  =>  S % CurrentSet )
    call CS % AccumulateBoundaryFluence ( dT  *  S % B ( iS ) )
    end associate !-- CS

    if ( associated ( T_AS ) ) call T_AS % Stop ( )

    !-- Cleanup

    end associate !-- K, etc.

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


  subroutine StoreSolution ( S, T_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    associate ( Y  =>  S % Solution )
    call StoreSolution_CS ( S, Y, T_Option )
    end associate !-- Y
  
  end subroutine StoreSolution


  subroutine SetSlope_CS ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    select type ( S )
      class is ( Step_RK_CS_Form )
    if ( allocated ( S % DivergenceTotal ) ) then
      allocate ( Slope_DFV_F_DT_Form :: K )
      select type ( K )
        class is ( Slope_DFV_F_DT_Form )
      call K % Initialize &
             ( S % RiemannSolver, S % DivergenceTotal, &
               IgnorabilityOption = S % IGNORABILITY )
      end select !-- K
    else if ( allocated ( S % DivergencePart ) ) then
      allocate ( Slope_DFV_F_DP_Form :: K )
      select type ( K )
        class is ( Slope_DFV_F_DP_Form )
      call K % Initialize &
             ( S % RiemannSolver, S % DivergencePart, &
               IgnorabilityOption = S % IGNORABILITY )
      end select !-- K
    end if
    end select !-- S

  end subroutine SetSlope_CS


  subroutine SetSlopeStage_CS ( S, K, iS )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( FieldSetForm ), intent ( out ), allocatable :: &
      K
    integer ( KDI ), intent ( in ) :: &
      iS

    character ( 1 ) :: &
      StageNumber

    select type ( S )
      class is ( Step_RK_CS_Form )
    associate &
      ( CS  =>  S % CurrentSet )

    write ( StageNumber, fmt = '(i1.1)' ) iS

    allocate ( K )
    call K % Initialize &
           ( CS % Atlas, &
             FieldOption = CS % Balanced, &
             NameOption = 'Slope_' // StageNumber, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = CS % nBalanced, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

    end associate !-- CS
    end select !-- S

  end subroutine SetSlopeStage_CS


  subroutine StoreSolution_CS ( S, Y, T_Option )

    class ( Step_RK_CS_Form ), intent ( inout ) :: &
      S
    type ( FieldSetForm ), intent ( in ) :: &
      Y
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_SC, &   !-- SolutionCopy
      T_CFB, &  !-- ComputeFromBalanced
      T_BC      !-- BoundaryConditions

    associate &
      ( CS_B  =>  S % Balanced, &
        CS    =>  S % CurrentSet )

    if ( present ( T_Option ) ) then
      T_SC   =>  PROGRAM_HEADER % Timer &
                   ( Handle = S % iTimer_SC, &
                     Name = trim ( S % Name ) // '_CpySltn', &
                     Level = T_Option % Level + 1 )
    else
      T_SC   =>  null ( )
    end if
    if ( associated ( T_SC ) ) call T_SC % Start ( )
    call  Y % Copy ( CS_B )
    if ( associated ( T_SC ) ) call T_SC % Stop ( )

    if ( present ( T_Option ) ) then
      T_CFB  =>  PROGRAM_HEADER % Timer &
                   ( Handle = S % iTimer_CFB, &
                     Name = trim ( S % Name ) // '_FrmBlncd', &
                     Level = T_Option % Level + 1 )
    else
      T_CFB  =>  null ( )
    end if
    if ( associated ( T_CFB ) ) then
      call T_CFB % Start ( )
      call CS % ComputeFromBalanced ( T_Option = T_CFB )
      call T_CFB % Stop ( )
    else
      call CS % ComputeFromBalanced ( )
    end if

    if ( present ( T_Option ) ) then
      T_BC   =>  PROGRAM_HEADER % Timer &
                   ( Handle = S % iTimer_BC, &
                     Name = trim ( S % Name ) // '_BndryCndtns', &
                     Level = T_Option % Level + 1 )
    else
      T_BC   =>  null ( )
    end if
    if ( associated ( T_BC ) ) call T_BC % Start ( )
    call CS % ApplyBoundaryConditions ( )
    if ( associated ( T_BC ) ) call T_BC % Stop ( )
  
    end associate !-- CS_B

  end subroutine StoreSolution_CS


end module Step_RK_CS__Form
