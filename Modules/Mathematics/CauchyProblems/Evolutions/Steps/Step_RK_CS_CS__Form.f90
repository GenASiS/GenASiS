module Step_RK_CS_CS__Form

  !-- Step_RungeKutta_CurrentSet_CurrentSet_Form

  use Basics
  use Fields
  use Slopes
  use Step_RK_H__Form
  use Step_RK_CS__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CS_CS_Form
    type ( FieldSet_BM_Form ), allocatable :: &
      Implicit_1, &
      Implicit_2   
    class ( Step_RK_CS_Form ), allocatable :: &
      Step_CS_1, &
      Step_CS_2
  contains
    procedure, private, pass :: &
      Initialize_CS_CS
    generic, public :: &
      Initialize => Initialize_CS_CS
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      SetCoarsening
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
  end type Step_RK_CS_CS_Form


contains


  subroutine Initialize_CS_CS &
               ( S, CS_1, CS_2, NameOption, ImplicitExplicitOption, &
                 nStagesOption )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS_1, &
      CS_2
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    character ( LDL ) :: &
      Name

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CS_CS'

    Name  =  trim ( CS_1 % Name ) // '_' // trim ( CS_2 % Name ) // '_Stp' 
    if ( present ( NameOption ) ) &
      Name  =  NameOption

    call S % Initialize_H &
           ( CS_1 % Atlas, &
             NameOption = Name, &
             ImplicitExplicitOption = ImplicitExplicitOption, &
             nStagesOption = nStagesOption )

    !-- Storage used in implicit solver

    allocate ( S % Implicit_1 )
    allocate ( S % Implicit_2 )
    associate &
      ( Y_S_1  =>  S % Implicit_1, &  !-- Y_Star
        Y_S_2  =>  S % Implicit_2 )
    call Y_S_1 % Initialize &
           ( CS_1 % Atlas, &
             FieldOption = CS_1 % Balanced, &
             NameOption = trim ( CS_1 % Name ) // '_Implicit', &
             DeviceMemoryOption = CS_1 % DeviceMemory, &
             DevicesCommunicateOption = CS_1 % DevicesCommunicate, &
             nFieldsOption = CS_1 % nBalanced, &
             IgnorabilityOption = CS_1 % IGNORABILITY + 1 )
    call Y_S_2 % Initialize &
           ( CS_2 % Atlas, &
             FieldOption = CS_2 % Balanced, &
             NameOption = trim ( CS_2 % Name ) // '_Implicit', &
             DeviceMemoryOption = CS_2 % DeviceMemory, &
             DevicesCommunicateOption = CS_2 % DevicesCommunicate, &
             nFieldsOption = CS_2 % nBalanced, &
             IgnorabilityOption = CS_2 % IGNORABILITY + 1 )
    end associate !-- Y_S_1, etc.

    !-- Steps

    if ( .not. allocated ( S % Step_CS_1 ) ) &
      allocate ( S % Step_CS_1 )
    if ( .not. allocated ( S % Step_CS_2 ) ) &
      allocate ( S % Step_CS_2 )

    call S % Step_CS_1 % Initialize &
           ( CS_1, NameOption, ImplicitExplicitOption, nStagesOption )
    call S % Step_CS_2 % Initialize &
           ( CS_2, NameOption, ImplicitExplicitOption, nStagesOption )

  end subroutine Initialize_CS_CS


  subroutine SetStream ( S, Sm )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    class ( Stream_BM_Form ), intent ( inout ) :: &
      Sm

    call S % Step_CS_1 % SetStream ( Sm )
    call S % Step_CS_2 % SetStream ( Sm )

  end subroutine SetStream


  subroutine SetCoarsening ( S, C, iS )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    class ( Coarsening_C_Form ), intent ( in ), target :: &
      C
    integer ( KDI ), intent ( in ) :: &
      iS

    select case ( iS )
    case ( 1 )
      call S % Step_CS_1 % SetCoarsening ( C )
    case ( 2 )
      call S % Step_CS_2 % SetCoarsening ( C )
    end select !-- iS

  end subroutine SetCoarsening


  subroutine Show_S ( S )

    class ( Step_RK_CS_CS_Form ), intent ( in ) :: &
      S

    call S % Step_RK_H_Form % Show ( )

    call S % Step_CS_1 % Show ( )
    call S % Step_CS_2 % Show ( )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Step_CS_2 ) ) &
      deallocate ( S % Step_CS_2 )
    if ( allocated ( S % Step_CS_1 ) ) &
      deallocate ( S % Step_CS_1 )
    if ( allocated ( S % Implicit_1 ) ) &
      deallocate ( S % Implicit_2 )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S

    call S % Step_CS_1 % LoadSolution ( )
    call S % Step_CS_2 % LoadSolution ( )

  end subroutine LoadSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    call S % Step_CS_1 % InitializeIntermediate ( iS )
    call S % Step_CS_2 % InitializeIntermediate ( iS )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, dT, iS, iK )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS, &
      iK

    call S % Step_CS_1 % IncrementIntermediate ( dT, iS, iK )
    call S % Step_CS_2 % IncrementIntermediate ( dT, iS, iK )

  end subroutine IncrementIntermediate


  subroutine StoreIntermediate ( S, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call S % Step_CS_1 % StoreIntermediate ( T_Option )
    call S % Step_CS_2 % StoreIntermediate ( T_Option )

  end subroutine StoreIntermediate


  subroutine ComputeUpdateImplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iNS, &  !-- iNonlinearSolve
      maxNS    

    if ( iS  ==  1 ) &
      return

    associate &
      ( S_1  =>  S % Step_CS_1, &
        S_2  =>  S % Step_CS_2 )
    associate &
      ( AA      =>  S_1 % AA ( iS ) % Value ( iS ), &          
        Y_I_1   =>  S_1 % Intermediate, &
        Y_I_2   =>  S_2 % Intermediate, &
        Y_S_1   =>  S   % Implicit_1, &
        Y_S_2   =>  S   % Implicit_2, &
        CS_B_1  =>  S_1 % Balanced, &
        CS_B_2  =>  S_2 % Balanced, &
        CS_1    =>  S_1 % CurrentSet, &
        CS_2    =>  S_2 % CurrentSet, &
        KK_1    =>  S_1 % SlopeImplicit, &
        KK_2    =>  S_2 % SlopeImplicit, &
        KK_1_Stage  =>  S_1 % SlopeStageImplicit ( iS ) % Element, &
        KK_2_Stage  =>  S_2 % SlopeStageImplicit ( iS ) % Element )

    if ( AA == 0.0_KDR ) &
      return

    !-- Upon entry, Y_I = Q_(I-1)
!-- FIXME: for maxNS > 1, an additional Y_S = Y_Star needed for iteration
    call Y_I_1 % Copy ( Y_S_1 )
    call Y_I_2 % Copy ( Y_S_2 )

    maxNS  =  20
      iNS  =  0
    do 

      iNS  =  iNS + 1
!call Show ( iNS, '>>> iNS' )

      call KK_1 % Compute ( dT )!, T_Option = T_CS )

      call KK_2 % Compute ( dT )!, T_Option = T_CS )
!call Show ( KK_2 % Storage_GS % Value ( :, 5 ), '>>> KK_2' )

      !-- Exit criterion
      if ( iNS  ==  maxNS ) exit

      call Y_S_1 % MultiplyAdd ( Y_I_1, KK_1, dT * AA )      
      call Y_S_1 % Copy ( CS_B_1 )
      call CS_1 % ComputeFromBalanced ( )

      call Y_S_2 % MultiplyAdd ( Y_I_2, KK_2, dT * AA )      
      call Y_S_2 % Copy ( CS_B_2 )
      call CS_2 % ComputeFromBalanced ( )

    end do !-- iNS

    !-- Assume SlopeImplicit local: no ghost exchange in loop above 
    call KK_1 % ExchangeGhostData ( )
    call KK_2 % ExchangeGhostData ( )

    call KK_1 % Copy ( KK_1_Stage )
    call KK_2 % Copy ( KK_2_Stage )

!-- FIXME: separate AccumulateSlope implicit and explicit
!    call S_1 % AccumulateSlope ( iS )
!    call S_2 % AccumulateSlope ( iS )
    
    !-- Upon exit, Y_I = Q_(I-1)  +  dT * AA ( iS ) * KK
    call Y_I_1 % MultiplyAdd ( KK_1, dT * AA )
    call Y_I_2 % MultiplyAdd ( KK_2, dT * AA )

    call S_1 % StoreIntermediate ( ) !T_Option )
    call S_2 % StoreIntermediate ( ) !T_Option )

    end associate !-- KK_1, etc.
    end associate !-- S_1, etc.

  end subroutine ComputeUpdateImplicit


  subroutine ComputeUpdateExplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    call S % Step_CS_1 % ComputeUpdateExplicit ( T, dT, iS, T_Option )
    call S % Step_CS_2 % ComputeUpdateExplicit ( T, dT, iS, T_Option )

  end subroutine ComputeUpdateExplicit


  subroutine IncrementSolution ( S, dT, iS )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % Step_CS_1 % IncrementSolution ( dT, iS )
    call S % Step_CS_2 % IncrementSolution ( dT, iS )

  end subroutine IncrementSolution


  subroutine StoreSolution ( S, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call S % Step_CS_1 % StoreSolution ( T_Option )
    call S % Step_CS_2 % StoreSolution ( T_Option )

  end subroutine StoreSolution


end module Step_RK_CS_CS__Form
