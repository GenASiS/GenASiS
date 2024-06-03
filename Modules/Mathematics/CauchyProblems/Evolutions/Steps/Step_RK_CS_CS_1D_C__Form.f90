module Step_RK_CS_CS_1D_C__Form

  !-- Step_RungeKutta_CurrentSet_CurrentSet_1D_Collected__Form

  use Basics
  use Fields
  use ImplicitDiagnostics_Form
  use Step_RK_H__Form
  use Step_RK_CS__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CS_CS_1D_C_Form
    integer ( KDI ) :: &
      nCurrentSets_1D, &
      MaxImplicitIterations, &
      MaxRelaxationIterations
    real ( KDR ) :: &
      ImplicitTolerance
    class ( ImplicitDiagnosticsForm ), dimension ( :, : ), allocatable :: &
      ImplicitDiagnostics
    class ( Step_RK_CS_Form ), allocatable :: &
      Step_CS
    class ( Step_RK_CS_Form ), dimension ( : ), allocatable :: &
      Step_CS_1D
  contains
    procedure, private, pass :: &
      Initialize_CS_CS_1D_C
    generic, public :: &
      Initialize => Initialize_CS_CS_1D_C
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_S
    final :: &
      Finalize
    procedure, public, pass :: &
      LoadSolution
    procedure, public, pass :: &
      InitializeIntermediate
    procedure, public, pass :: &
      IncrementIntermediate
    procedure, public, pass :: &
      StoreIntermediate
    procedure, public, pass :: &
      ComputeUpdateImplicit
    procedure, public, pass :: &
      ComputeUpdateExplicit
    procedure, public, pass :: &
      IncrementSolution
    procedure, public, pass :: &
      StoreSolution
    procedure, public, pass :: &
      SolveUpdateImplicit
  end type Step_RK_CS_CS_1D_C_Form


contains


  subroutine Initialize_CS_CS_1D_C &
               ( S, CS_1D, CS, NameOption, ImplicitExplicitOption, &
                 nStagesOption )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), dimension ( : ), intent ( in ), target :: &
      CS_1D
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    integer ( KDI ) :: &
      iCS  !-- iCurrentSet
    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CS_CS_1D_C'

    Name  =  'CS_CS_1D_C_Stp' 
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call S % Initialize_H &
           ( CS % Atlas, &
             NameOption = Name, &
             ImplicitExplicitOption = ImplicitExplicitOption, &
             nStagesOption = nStagesOption )

    !-- Storage used in implicit solver

    S % nCurrentSets_1D  =  size ( CS_1D )

    S % MaxImplicitIterations  =  50
    call PROGRAM_HEADER % GetParameter &
           ( S % MaxImplicitIterations, 'MaxImplicitIterations' )

    S % MaxRelaxationIterations  =  10
    call PROGRAM_HEADER % GetParameter &
           ( S % MaxRelaxationIterations, 'MaxRelaxationIterations' )

    S % ImplicitTolerance  =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter &
           ( S % ImplicitTolerance, 'ImplicitTolerance' )

    !-- Steps

    if ( .not. allocated ( S % Step_CS ) ) &
      allocate ( S % Step_CS )
    call S % Step_CS % Initialize &
           ( CS, NameOption, ImplicitExplicitOption, nStagesOption )

    if ( .not. allocated ( S % Step_CS_1D ) ) &
      allocate ( S % Step_CS_1D ( S % nCurrentSets_1D ) )
    do iCS  =  1,  size ( CS_1D )
      call S % Step_CS_1D ( iCS ) % Initialize &
             ( CS_1D ( iCS ), NameOption, ImplicitExplicitOption, &
               nStagesOption )
    end do !-- iCS

  end subroutine Initialize_CS_CS_1D_C


  subroutine SetStream ( S, Sm )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    class ( Stream_BM_Form ), intent ( inout ) :: &
      Sm

    integer ( KDI ) :: &
      iCS, &
      iS

    call S % Step_CS % SetStream ( Sm )

    do iCS  =  1,  S % nCurrentSets_1D

      call S % Step_CS_1D ( iS ) % SetStream ( Sm )

      do iS = 2, S % nStages
        call Sm % AddFieldSet ( S % ImplicitDiagnostics ( iS, iCS ) )
      end do !-- iS
    
    end do !-- iCS

  end subroutine SetStream


  subroutine Show_S ( S )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( in ) :: &
      S

    integer ( KDI ) :: &
      iCS, &
      iS

    call S % Step_RK_H_Form % Show ( )

    call Show ( S % MaxImplicitIterations, 'MaxImplicitIterations', &
                S % IGNORABILITY )
    call Show ( S % MaxRelaxationIterations, 'MaxRelaxationIterations', &
                S % IGNORABILITY )
    call Show ( S % ImplicitTolerance, 'ImplicitTolerance', &
                S % IGNORABILITY )

    call S % Step_CS % Show ( )

    do iCS  =  1,  S % nCurrentSets_1D

      call S % Step_CS_1D ( iCS ) % Show ( )

      do iS = 2, S % nStages
        call S % ImplicitDiagnostics ( iS, iCS ) % Show ( )
      end do !-- iS

    end do !-- iCS

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Step_CS_1D ) ) &
      deallocate ( S % Step_CS_1D )
    if ( allocated ( S % Step_CS ) ) &
      deallocate ( S % Step_CS )
    if ( allocated ( S % ImplicitDiagnostics ) ) &
      deallocate ( S % ImplicitDiagnostics )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % LoadSolution ( )
    
    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % LoadSolution ( )
    end do

  end subroutine LoadSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % InitializeIntermediate ( iS )

    do iCS  =  1, S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % InitializeIntermediate ( iS )
    end do 

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, dT, iS, iK )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS, &
      iK

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % IncrementIntermediate ( dT, iS, iK )

    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % IncrementIntermediate ( dT, iS, iK )
    end do

  end subroutine IncrementIntermediate


  subroutine StoreIntermediate ( S, T_Option )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % StoreIntermediate ( T_Option )

    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % StoreIntermediate ( T_Option )
    end do

  end subroutine StoreIntermediate


  subroutine ComputeUpdateImplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iCS

    if ( iS  ==  1 ) &
      return !-- Member AA below is indexed starting with 2

    associate &
      ( S_CS     =>  S % Step_CS, &
        S_CS_1D  =>  S % Step_CS_1D ( : ), &
        AA       =>  S % Step_CS % AA ( iS ) % Value ( iS ) )

    if ( AA == 0.0_KDR ) &
      return  !-- No use computing an implicit update if it will be added
              !   with zero weight!

    !-- Upon entry, Y_I = Q_(I-1)

    !-- Solve for KK
    call S % SolveUpdateImplicit ( T, dT, iS )

    !-- CS

    associate &
      ( Y_I_CS  =>  S_CS % Intermediate, &
        KK_CS   =>  S_CS % SlopeStageImplicit ( iS ) % Element )

    !-- Assume SlopeImplicit local: no ghost exchange in implicit solve 
    call KK_CS % ExchangeGhostData ( )

!-- FIXME: separate AccumulateSlope implicit and explicit for diagnostic I/O
!    call S_CS % AccumulateSlope ( iS )

    !-- Upon exit, Y_I = Q_(I-1)  +  dT * AA ( iS ) * KK
    call Y_I_CS % MultiplyAdd ( KK_CS, dT * AA )

    call S_CS % StoreIntermediate ( ) !T_Option )

    end associate !-- Y_I_CS, etc.

    !-- CS_1D

    do iCS  =  1,  S % nCurrentSets_1D

      associate &
        ( Y_I_CS_1D  =>  S_CS_1D ( iCS ) % Intermediate, &
           KK_CS_1D  =>  S_CS_1D ( iCS ) % SlopeStageImplicit ( iS ) % Element )

      !-- Assume SlopeImplicit local: no ghost exchange in implicit solve 
      call KK_CS_1D % ExchangeGhostData ( )

  !-- FIXME: separate AccumulateSlope implicit and explicit for diagnostic I/O
  !    call S_CS % AccumulateSlope ( iS )

      !-- Upon exit, Y_I = Q_(I-1)  +  dT * AA ( iS ) * KK
      call Y_I_CS_1D % MultiplyAdd ( KK_CS_1D, dT * AA )

      call S_CS_1D ( iCS ) % StoreIntermediate ( ) !T_Option )

      end associate !-- Y_I_CS, etc.

    end do !-- iCS

    end associate !-- S_CS, etc.

  end subroutine ComputeUpdateImplicit


  subroutine ComputeUpdateExplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % ComputeUpdateExplicit ( T, dT, iS, T_Option )

    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % ComputeUpdateExplicit &
             ( T, dT, iS, T_Option )
    end do 

  end subroutine ComputeUpdateExplicit


  subroutine IncrementSolution ( S, dT, iS )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % IncrementSolution ( dT, iS )

    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % IncrementSolution ( dT, iS )
    end do

  end subroutine IncrementSolution


  subroutine StoreSolution ( S, T_Option )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iCS

    call S % Step_CS % StoreSolution ( T_Option )

    do iCS  =  1,  S % nCurrentSets_1D
      call S % Step_CS_1D ( iCS ) % StoreSolution ( T_Option )
    end do

  end subroutine StoreSolution


  subroutine SolveUpdateImplicit  ( S, T, dT, iS )

    class ( Step_RK_CS_CS_1D_C_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'SolveUpdateImplicit should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_CS_CS_1D_C_Form', 'module', CONSOLE % WARNING )
    call Show ( 'SolveUpdateImplicit', 'subroutine', CONSOLE % WARNING )

  end subroutine SolveUpdateImplicit


end module Step_RK_CS_CS_1D_C__Form
