module Step_RK_H__Form
  
  !-- Step_RungeKutta_Header_Form

  use Basics
  use Manifolds
  use Fields
  use Slopes
  
  implicit none
  private

  type, public :: Step_RK_H_Form
    integer ( KDI ) :: &
      IGNORABILITY = 0, &
      nStages = 0
    integer ( KDI ) :: &
      iTimer      = 0, &
      iTimer_LI   = 0, &  !-- LoadInitial
      iTimer_II   = 0, &  !-- InitializeIntermediate
      iTimer_II_A = 0, &  !-- IncrementIntermediate
      iTimer_SI   = 0, &  !-- StoreIntermediate
      iTimer_CS   = 0, &  !-- ComputeStage
      iTimer_IS_B = 0, &  !-- IncrementSolution
      iTimer_SS   = 0, &  !-- StoreSolution
      iTimer_AS   = 0     !-- AccumulateSlope
    real ( KDR ), dimension ( : ), allocatable :: &
      C,  CC, &  !-- RungeKutta nodes (explicit, implicit)
      B,  BB, &  !-- RungeKutta weights (explicit, implicit)
      BE, BBE    !-- RungeKutta weights of embedded lower order solution
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      A, AA  !-- RungeKutta matrix (explicit, implicit)
    logical ( KDL ) :: &
      ImplicitExplicit, &
      EmbeddedMethod
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    class ( Slope_H_Form ), allocatable :: &
      Slope, &
      SlopeSum
    type ( FieldSet_BM_Element ), dimension ( : ), allocatable :: &
      SlopeStageImplicit, &  !-- storage of KK for various stages
      SlopeStageExplicit     !-- storage of K for various stages
    procedure ( SS ), pointer :: &
      SetSlopeImplicit => null ( ), &
      SetSlopeExplicit => null ( )
    procedure ( SSS ), pointer :: &
      SetSlopeStage => null ( )  !-- used to initialize either K or KK storage
  contains
    procedure, public, pass :: &
      Initialize_H  !-- Do not overload: needs overriding of SetSlope
    procedure, public, pass :: &
      SetStream_H
    procedure, public, pass :: &
      SetStream => SetStream_H
    procedure, public, pass :: &
      Show => Show_S
    procedure, public, pass :: &
      Timer
    procedure, public, pass :: &
      TimerStoreSolution
    procedure, public, pass :: &
      TimerAccumulateSlope
    procedure, public, pass :: &
      Compute
    procedure, public, pass :: &
      AccumulateSlope
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
      ComputeStageImplicit
    procedure, public, pass :: &
      ComputeStageExplicit
    procedure, public, pass :: &
      IncrementSolution
    procedure, public, pass :: &
      StoreSolution
  end type Step_RK_H_Form

  interface

    subroutine SS ( S, K )
      use Basics
      use Slopes
      import Step_RK_H_Form
      implicit none
      class ( Step_RK_H_Form ), intent ( in ) :: &
        S
      class ( Slope_H_Form ), intent ( out ), allocatable :: &
        K
    end subroutine SS

    subroutine SSS ( S, K, iS )
      use Basics
      use Fields
      import Step_RK_H_Form
      implicit none
      class ( Step_RK_H_Form ), intent ( in ) :: &
        S
      class ( FieldSet_BM_Form ), intent ( out ), allocatable :: &
        K
      integer ( KDI ), intent ( in ) :: &
        iS
    end subroutine SSS

  end interface

    private :: &
      SetCoefficientsImplicitExplicit, &
      SetCoefficientsExplicit, &
      SetSlope_H, &
      SetSlopeStage_H

contains


  subroutine Initialize_H &
               ( S, Atlas, NameOption, ImplicitExplicitOption, nStagesOption )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    class ( Atlas_H_Form ), intent ( in ), target :: &
      Atlas
    character ( * ), intent ( in ), optional :: &
      NameOption
    logical ( KDL ), intent ( in ), optional :: &
      ImplicitExplicitOption
    integer ( KDI ), intent ( in ), optional :: &
      nStagesOption

    integer ( KDI ) :: &
      iS  !-- iStage

    S % IGNORABILITY  =  CONSOLE % INFO_1

    S % ImplicitExplicit  =  .false.
    if ( present ( ImplicitExplicitOption ) ) &
      S % ImplicitExplicit  =  ImplicitExplicitOption

    S % EmbeddedMethod  =  .false.

    if ( S % ImplicitExplicit ) then
      S % nStages  =  3
    else
      S % nStages  =  2
    end if
    if ( present ( nStagesOption ) ) &
      S % nStages  =  nStagesOption
    associate ( nS  =>  S % nStages )

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    S % Name = 'Step'
    if ( present ( NameOption ) ) &
      S % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    if ( S % ImplicitExplicit ) then
      call SetCoefficientsImplicitExplicit ( S )
    else
      call SetCoefficientsExplicit ( S )
    end if

    S % Atlas  =>  Atlas

    if ( .not. associated ( S % SetSlopeExplicit ) ) &
      S % SetSlopeExplicit  =>  SetSlope_H
    if ( .not. associated ( S % SetSlopeStage ) ) &
      S % SetSlopeStage  =>  SetSlopeStage_H

    call S % SetSlopeExplicit ( S % Slope )

    if ( S % ImplicitExplicit ) then
      allocate ( S % SlopeStageImplicit ( nS ) )
      do iS  =  1,  nS
        call S % SetSlopeStage ( S % SlopeStageImplicit ( iS ) % Element, iS )
      end do !-- iS
    end if

    allocate ( S % SlopeStageExplicit ( nS ) )
    do iS  =  1,  nS
      call S % SetSlopeStage ( S % SlopeStageExplicit ( iS ) % Element, iS )
    end do !-- iS

    end associate !-- nS

  end subroutine Initialize_H


  subroutine SetStream_H ( S, Sm )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    class ( Stream_BM_Form ), intent ( inout ) :: &
      Sm

    integer ( KDI ) :: &
      iS  !-- iStage

    if ( .not. allocated ( S % SlopeSum ) ) then
      call S % SetSlopeExplicit ( S % SlopeSum )
      associate ( K_Sum  =>  S % SlopeSum )
      call K_Sum % SetStream ( Sm )
      end associate !-- K_Sum
    end if !-- allocated SlopeSum
    
  end subroutine SetStream_H


  subroutine Show_S ( S )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S

   integer ( KDI ) :: &
     iA, &
     iS  !-- iStage
   character ( 1 ) :: &
     Index
   character ( LDL ), dimension ( : ), allocatable :: &
     TypeWord

    call Split ( S % Type, ' ', TypeWord )
    call Show ( trim ( TypeWord ( 2 ) ) // ' Parameters', S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )
    call Show ( S % ImplicitExplicit, 'ImplicitExplicit', S % IGNORABILITY )
    call Show ( S % EmbeddedMethod, 'EmbeddedMethod', S % IGNORABILITY )

    call Show ( S % nStages, 'nStages', S % IGNORABILITY )

    !-- Explicit Butcher tableau

    do iA  =  2, S % nStages
      write ( Index, fmt = '(i1.1)' ) iA
      call Show ( S % A ( iA ) % Value, 'A ( ' // Index // ' )' )
    end do !-- iA

    call Show ( S % B, 'B' )
    if ( S % EmbeddedMethod ) &
      call Show ( S % BE, 'BE' )

    call Show ( S % C, 'C', lRealOption = 2 )

    if ( S % ImplicitExplicit ) then

      !-- Implicit Butcher tableau

      do iA  =  2, S % nStages
        write ( Index, fmt = '(i1.1)' ) iA
        call Show ( S % AA ( iA ) % Value, 'AA ( ' // Index // ' )' )
      end do !-- iA

      call Show ( S % BB, 'BB' )
      if ( S % EmbeddedMethod ) &
        call Show ( S % BBE, 'BBE' )

      call Show ( S % CC, 'CC', lRealOption = 2 )

    end if !-- ImplicitExplicit

    !-- Slopes

    do iS  =  1, S % nStages
      call S % SlopeStageExplicit ( iS ) % Element % Show ( )
    end do !-- iS
    if ( allocated ( S % SlopeSum ) ) &
      call S % SlopeSum % Show ( )

  end subroutine Show_S


  function Timer ( S, Level ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimer, &
               Name = S % Name, &
               Level = Level )

  end function Timer


  function TimerStoreSolution ( S, Level ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimer_SS, &
               Name = trim ( S % Name ) // '_StrSltn', &
               Level = Level )

  end function TimerStoreSolution


  function TimerAccumulateSlope ( S, Level ) result ( T )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = S % iTimer_AS, &
               Name = trim ( S % Name) // '_AccmltSlp', &
               Level = Level )

  end function TimerAccumulateSlope


  subroutine Compute ( S, T, dT, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T, &
      dT
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iK     !-- iIncrement
    type ( TimerForm ), pointer :: &
      T_LI, &
      T_II, &
      T_II_A, &
      T_SI, &
      T_CS, &
      T_IS_B, &
      T_SS

    call Show ( 'Computing ' // trim ( S % Type ), S % IGNORABILITY + 2 )
    call Show ( S % Name, 'Name', S % IGNORABILITY + 2 )

    if ( present ( T_Option ) ) then
      T_LI    =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_LI, &
                      Name = trim ( S % Name ) // '_LdIntl', &
                      Level = T_Option % Level + 1 )
      T_II    =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_II, &
                      Name = trim ( S % Name ) // '_IntlzIntmdt', &
                      Level = T_Option % Level + 1 )
      T_II_A  =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_II_A, &
                      Name = trim ( S % Name ) // '_IncrmntIntmdt', &
                      Level = T_Option % Level + 1 )
      T_SI    =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_SI, &
                      Name = trim ( S % Name ) // '_StrIntmdt', &
                      Level = T_Option % Level + 1 )
    else
      T_LI    =>  null ( )
      T_II    =>  null ( )
      T_II_A  =>  null ( )
      T_SI    =>  null ( )
    end if

    !-- Set  Y  =  Y_N  (old value)
    if ( associated ( T_LI ) ) call T_LI % Start ( )
    call S % LoadSolution ( )
    if ( associated ( T_LI ) ) call T_LI % Stop ( )

    !-- Compute stages

    do iS = 1, S % nStages

      call Show ( 'Computing a stage', S % IGNORABILITY + 3 )
      call Show ( iS, 'iStage', S % IGNORABILITY + 3 )

      !-- Set  Q_(I-1)  =  Y
      if ( associated ( T_II ) ) call T_II % Start ( )
      call S % InitializeIntermediate ( iS )
      if ( associated ( T_II ) ) call T_II % Stop ( )

      !-- Increment Q_(I-1) with previous updates
      do iK = 1, iS - 1
        !-- Set Q_(I-1)  =  Q_(I-1)  +  dT * A  * K  ( iK )  
        !                            +  dT * AA * KK ( iK ) 
        associate ( A  =>  S % A ( iS ) % Value ( iK ) )
        if ( associated ( T_II_A ) ) call T_II_A % Start ( )
        if ( S % ImplicitExplicit ) then
          associate ( AA  =>  S % AA ( iS ) % Value ( iK ) )          
          call S % IncrementIntermediate ( A, dT, iK, AA_Option = AA )
          end associate !-- AA
        else
          call S % IncrementIntermediate ( A, dT, iK )
        end if
        if ( associated ( T_II_A ) ) call T_II_A % Stop ( )
        end associate !-- A
      end do !-- iK

      !-- Obtain Y_(I) and KK ( iS ) = dY/dT_Implicit ( Y_(I) )
      !   from nonlinear solve Y_(I) = Q_(I-1) + dt A_II KK ( iS )

      !-- Store Y_(I) back to Y
      if ( associated ( T_SI ) ) call T_SI % Start ( )
        call S % StoreIntermediate ( T_Option = T_SI )
      if ( associated ( T_SI ) ) call T_SI % Stop ( )

      !-- Compute K ( iS )  =  dY/dT_Explicit ( Y_(I) )
      if ( present ( T_Option ) ) then
        T_CS  =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_CS, &
                      Name = trim ( S % Name ) // '_CmptStg', &
                      Level = T_Option % Level + 1 )
        call T_CS % Start ( )
        call S % ComputeStageExplicit ( T, dT, iS, T_Option = T_CS )
        call T_CS % Stop ( )
      else
        call S % ComputeStageExplicit ( T, dT, iS )
      end if

    end do !-- iS

    !-- Assemble stages

    if ( present ( T_Option ) ) then
      T_IS_B  =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_IS_B, &
                      Name = trim ( S % Name ) // '_IncrmntSltn', &
                      Level = T_Option % Level + 1 )
    else
      T_IS_B  =>  null ( )
    end if
    if ( associated ( T_IS_B ) ) call T_IS_B % Start ( )
    do iS = 1, S % nStages
      associate &
        ( B   =>  S % B ( iS ), &
          BE  =>  S % BE ( iS ) )
      !-- Set Y  =  Y  +  dT * B * K ( iS )
      if ( B /=  0.0_KDR ) & 
        call S % IncrementSolution ( B, BE, dT, iS )
      end associate !-- B, BE
    end do !-- iS
    if ( associated ( T_IS_B ) ) call T_IS_B % Stop ( )

    !-- Set Y_(N+1)  =  Y

    if ( present ( T_Option ) ) then
      T_SS  =>  S % TimerStoreSolution ( Level = T_Option % Level + 1 )
      call T_SS % Start ( )
      call S % StoreSolution ( T_Option = T_SS )
      call T_SS % Stop ( )
    else
      call S % StoreSolution ( )
    end if

  end subroutine Compute


  subroutine AccumulateSlope ( S, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    if ( .not. allocated ( S % SlopeSum ) ) &
      return

    associate &
      ( K_Sum  =>  S % SlopeSum, &
        K      =>  S % Slope )

    if ( iS == 1 ) &
      call K_Sum % ClearRecursive ( )
    call K_Sum % MultiplyAddRecursive ( K, S % B ( iS ) )

    end associate !-- K_Sum, etc.

  end subroutine AccumulateSlope


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % SlopeStageExplicit ) ) &
      deallocate ( S % SlopeStageExplicit )
    if ( allocated ( S % SlopeStageImplicit ) ) &
      deallocate ( S % SlopeStageImplicit )
    if ( allocated ( S % SlopeSum ) ) &
      deallocate ( S % SlopeSum )
    if ( allocated ( S % Slope ) ) &
      deallocate ( S % Slope )
    if ( allocated ( S % AA ) ) &
      deallocate ( S % AA )
    if ( allocated ( S % A ) ) &
      deallocate ( S % A )
    if ( allocated ( S % BBE ) ) &
      deallocate ( S % BBE )
    if ( allocated ( S % BE ) ) &
      deallocate ( S % BE )
    if ( allocated ( S % BB ) ) &
      deallocate ( S % BB )
    if ( allocated ( S % B ) ) &
      deallocate ( S % B )
    if ( allocated ( S % CC ) ) &
      deallocate ( S % CC )
    if ( allocated ( S % C ) ) &
      deallocate ( S % C )

    if ( S % Name == '' ) &
      return

    call Show ( 'Finalizing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    call Show ( 'LoadSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'LoadSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine LoadSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'InitializeIntermediate should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'InitializeIntermediate', 'subroutine', CONSOLE % WARNING )

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK, AA_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK
    real ( KDR ), intent ( in ), optional :: &
      AA_Option

    call Show ( 'IncrementIntermediate should be overridden', &
                CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'IncrementIntermediate', 'subroutine', CONSOLE % WARNING )

  end subroutine IncrementIntermediate


  subroutine StoreIntermediate ( S, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call Show ( 'StoreIntermediate should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'StoreIntermediate', 'subroutine', CONSOLE % WARNING )

  end subroutine StoreIntermediate


  subroutine ComputeStageImplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    call Show ( 'ComputeStageImplicit should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'ComputeStageImplicit', 'subroutine', CONSOLE % WARNING )

  end subroutine ComputeStageImplicit


  subroutine ComputeStageExplicit ( S, T, dT, iS, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    call Show ( 'ComputeStageExplicit should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'ComputeStageExplicit', 'subroutine', CONSOLE % WARNING )

  end subroutine ComputeStageExplicit


  subroutine IncrementSolution ( S, B, BE, dT, iS )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
       BE, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    call Show ( 'IncrementSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'IncrementSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine IncrementSolution


  subroutine StoreSolution ( S, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    call Show ( 'StoreSolution should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'StoreSolution', 'subroutine', CONSOLE % WARNING )

  end subroutine StoreSolution


  subroutine SetCoefficientsExplicit ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iStage
    real ( KDR ), dimension ( :, : ), allocatable :: &
      A  !-- RungeKutta matrix, explicit
    real ( KDR ), dimension ( : ), allocatable :: &
      B, &  !-- RungeKutta weights, explicit
      BE    !-- RungeKutta weights, embedded lower order, explicit
    real ( KDR ), dimension ( : ), allocatable :: &
      C  !-- RungeKutta nodes, explicit

    associate ( nS  =>  S % nStages )

    select case ( nS )
    case ( 2 )
      allocate ( A ( 2 : 2, 1 : 1 ) )
      A           =  0.0_KDR
      A ( 2, 1 )  =  1.0_KDR
      allocate ( B ( 1 : 2 ) )
      B ( 1 )  =  0.5_KDR
      B ( 2 )  =  0.5_KDR
      allocate ( C ( 2 : 2 ) )
      C ( 2 )  =  1.0_KDR
      allocate ( BE ( 1 : 2 ) )
      BE ( 1 )  =  1.0_KDR
      BE ( 2 )  =  0.0_KDR
      S % EmbeddedMethod  =  .true.
    case ( 3 )
      allocate ( A ( 2 : 3, 1 : 2 ) )
      A           =   0.0_KDR
      A ( 2, 1 )  =   0.5_KDR
      A ( 3, 1 )  =  -1.0_KDR
      A ( 3, 2 )  =   2.0_KDR
      allocate ( B ( 1 : 3 ) )
      B ( 1 )  =  1.0_KDR / 6.0_KDR
      B ( 2 )  =  2.0_KDR / 3.0_KDR
      B ( 3 )  =  1.0_KDR / 6.0_KDR
      allocate ( C ( 2 : 3 ) )
      C ( 2 )  =  0.5_KDR
      C ( 3 )  =  1.0_KDR
      allocate ( BE ( 1 : 3 ) )
      BE  =  0.0_KDR
      S % EmbeddedMethod  =  .false.
    case ( 4 )
      allocate ( A ( 2 : 4, 1 : 3 ) )
      A           =  0.0_KDR
      A ( 2, 1 )  =  0.5_KDR
      A ( 3, 1 )  =  0.0_KDR
      A ( 4, 1 )  =  0.0_KDR
      A ( 3, 2 )  =  0.5_KDR
      A ( 4, 2 )  =  0.0_KDR
      A ( 4, 3 )  =  1.0_KDR
      allocate ( B ( 1 : 4 ) )
      B ( 1 )  =  1.0_KDR / 6.0_KDR
      B ( 2 )  =  1.0_KDR / 3.0_KDR
      B ( 3 )  =  1.0_KDR / 3.0_KDR
      B ( 4 )  =  1.0_KDR / 6.0_KDR
      allocate ( C ( 2 : 4 ) )
      C ( 2 )  =  0.5_KDR
      C ( 3 )  =  0.5_KDR
      C ( 4 )  =  1.0_KDR
      allocate ( BE ( 1 : 4 ) )
      BE  =  0.0_KDR
      S % EmbeddedMethod  =  .false.
    case default
      call Show ( 'RungeKutta order not implemented', CONSOLE % ERROR )
      call Show ( S % nStages, 'nStages', CONSOLE % ERROR )
      call Show ( 'Step_RK_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCoefficientsExplicit', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- nStages

    allocate ( S % A  ( 2 : nS ) )
    do iS  =  2,  nS
      call S % A  ( iS ) % Initialize ( iS - 1 )
      S % A  ( iS ) % Value  =  A ( iS, 1 : iS - 1 )
    end do !-- iS

    allocate ( S % B  ( nS ) )
    S % B   =  B

    allocate ( S % C  ( 2 : nS ) )
    S % C   =  C

    allocate ( S % BE  ( nS ) )
    S % BE   =  BE
 
    end associate !-- nS

  end subroutine SetCoefficientsExplicit


  subroutine SetCoefficientsImplicitExplicit ( S )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iStage
    real ( KDR ) :: &
      A_32
    real ( KDR ), dimension ( :, : ), allocatable :: &
      A, AA  !-- RungeKutta matrix, explicit
    real ( KDR ), dimension ( : ), allocatable :: &
      B,  BB, &  !-- RungeKutta weights, explicit
      BE, BBE    !-- RungeKutta weights, embedded lower order, explicit
    real ( KDR ), dimension ( : ), allocatable :: &
      C, CC  !-- RungeKutta nodes, explicit

    associate ( nS  =>  S % nStages )

    select case ( nS )
    case ( 3 )
      !-- Giraldo et al. 2013, SIAM J. Sci. Comput. 35, B1162
      allocate ( A ( 2 : 3, 1 : 2 ) )
      A_32        =  ( 3.0_KDR  +  2.0_KDR * sqrt ( 2.0_KDR ) )  /  6.0_KDR
      A           =  0.0_KDR
      A ( 2, 1 )  =  2.0_KDR  -  sqrt ( 2.0_KDR )
      A ( 3, 1 )  =  1.0_KDR  -  A_32
      A ( 3, 2 )  =  A_32
      allocate ( AA ( 2 : 3, 1 : 3 ) )
      AA           =  0.0_KDR
      AA ( 2, 1 )  =  1.0_KDR  -  1.0_KDR / sqrt ( 2.0_KDR )
      AA ( 2, 2 )  =  1.0_KDR  -  1.0_KDR / sqrt ( 2.0_KDR )
      AA ( 3, 1 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      AA ( 3, 2 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      AA ( 3, 3 )  =  1.0_KDR  -  1.0_KDR / sqrt ( 2.0_KDR )
      allocate ( B ( 1 : 3 ) )
      B ( 1 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      B ( 2 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      B ( 3 )  =  1.0_KDR  -  1.0_KDR / sqrt ( 2.0_KDR )
      allocate ( BB ( 1 : 3 ) )
      BB ( 1 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      BB ( 2 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      BB ( 3 )  =  1.0_KDR  -  1.0_KDR / sqrt ( 2.0_KDR )
      allocate ( C ( 2 : 3 ) )
      C ( 2 )  =  2.0_KDR  -  sqrt ( 2.0_KDR )
      C ( 3 )  =  1.0_KDR
      allocate ( CC ( 2 : 3 ) )
      CC ( 2 )  =  2.0_KDR  -  sqrt ( 2.0_KDR )
      CC ( 3 )  =  1.0_KDR
      allocate ( BE ( 1 : 3 ) )
      BE ( 1 )  =  ( 4.0_KDR - sqrt ( 2.0_KDR ) ) / 8.0_KDR
      BE ( 2 )  =  ( 4.0_KDR - sqrt ( 2.0_KDR ) ) / 8.0_KDR
      BE ( 3 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      allocate ( BBE ( 1 : 3 ) )
      BBE ( 1 )  =  ( 4.0_KDR - sqrt ( 2.0_KDR ) ) / 8.0_KDR
      BBE ( 2 )  =  ( 4.0_KDR - sqrt ( 2.0_KDR ) ) / 8.0_KDR
      BBE ( 3 )  =  1.0_KDR / ( 2.0_KDR * sqrt ( 2.0_KDR ) )
      S % EmbeddedMethod  =  .true.
    case default
      call Show ( 'RungeKutta order not implemented', CONSOLE % ERROR )
      call Show ( S % nStages, 'nStages', CONSOLE % ERROR )
      call Show ( 'Step_RK_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCoefficientsImplicitExplicit', 'subroutine', &
                  CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- nStages

    allocate ( S % A  ( 2 : nS ) )
    allocate ( S % AA ( 2 : nS ) )
    do iS  =  2,  nS
      call S % A  ( iS ) % Initialize ( iS - 1 )
      call S % AA ( iS ) % Initialize ( iS )
      S % A  ( iS ) % Value  =  A  ( iS, 1 : iS - 1 )
      S % AA ( iS ) % Value  =  AA ( iS, 1 : iS )
    end do !-- iS

    allocate ( S % B  ( nS ) )
    allocate ( S % BB ( nS ) )
    S % B   =  B
    S % BB  =  BB

    allocate ( S % C  ( 2 : nS ) )
    allocate ( S % CC ( 2 : nS ) )
    S % C   =  C
    S % CC  =  CC

    allocate ( S % BE  ( nS ) )
    allocate ( S % BBE ( nS ) )
    S % BE   =  BE
    S % BBE  =  BBE

    end associate !-- nS

  end subroutine SetCoefficientsImplicitExplicit


  subroutine SetSlope_H ( S, K )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( Slope_H_Form ), intent ( out ), allocatable :: &
      K

    allocate ( Slope_H_Form :: K )
    associate ( A  =>  S % Atlas )
    call K % Initialize &
           ( A, &
             NameOption = 'Slope', &
             IgnorabilityOption = A % IGNORABILITY )
    end associate !-- A

  end subroutine SetSlope_H


  subroutine SetSlopeStage_H ( S, K, iS )

    class ( Step_RK_H_Form ), intent ( in ) :: &
      S
    class ( FieldSet_BM_Form ), intent ( out ), allocatable :: &
      K
    integer ( KDI ), intent ( in ) :: &
      iS

    character ( 1 ) :: &
      StageNumber

    write ( StageNumber, fmt = '(i1.1)' ) iS

    allocate ( K )
    associate ( A  =>  S % Atlas )
    call K % Initialize &
           ( A, &
             NameOption = 'Slope_' // StageNumber, &
             IgnorabilityOption = A % IGNORABILITY )
    end associate !-- A

  end subroutine SetSlopeStage_H


end module Step_RK_H__Form
