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
      ComputeStageExplicit
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
             nStagesOption = nStagesOption )

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


  subroutine IncrementIntermediate ( S, A, dT, iK, AA_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK
    real ( KDR ), intent ( in ), optional :: &
      AA_Option

    call S % Step_CS_1 % IncrementIntermediate ( A, dT, iK, AA_Option )
    call S % Step_CS_2 % IncrementIntermediate ( A, dT, iK, AA_Option )

  end subroutine IncrementIntermediate


  subroutine StoreIntermediate ( S, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call S % Step_CS_1 % StoreIntermediate ( T_Option )
    call S % Step_CS_2 % StoreIntermediate ( T_Option )

  end subroutine StoreIntermediate


  subroutine ComputeStageExplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    call S % Step_CS_1 % ComputeStageExplicit ( T, dT, iS, T_Option )
    call S % Step_CS_2 % ComputeStageExplicit ( T, dT, iS, T_Option )

  end subroutine ComputeStageExplicit


  subroutine IncrementSolution ( S, B, BE, dT, iS )

    class ( Step_RK_CS_CS_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
       BE, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call S % Step_CS_1 % IncrementSolution ( B, BE, dT, iS )
    call S % Step_CS_2 % IncrementSolution ( B, BE, dT, iS )

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
