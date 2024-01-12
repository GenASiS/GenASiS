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
      C, &  !-- RungeKutta nodes
      B, &  !-- RungeKutta weights
      BE    !-- RungeKutta weights of embedded lower order solution for error
    type ( Real_1D_Form ), dimension ( : ), allocatable :: &
      A  !-- RungeKutta matrix
    character ( LDF ) :: &
      Type = '', &
      Name = ''
    class ( Atlas_H_Form ), pointer :: &
      Atlas => null ( )
    class ( Slope_H_Form ), allocatable :: &
      Slope, &
      SlopeSum
    type ( FieldSet_BM_Element ), dimension ( : ), allocatable :: &
      SlopeStage
    procedure ( SS ), pointer :: &
      SetSlope => null ( )
    procedure ( SSS ), pointer :: &
      SetSlopeStage => null ( )
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
      ComputeStage
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
      SetCoefficientsExplicit, &
      SetSlope_H, &
      SetSlopeStage_H

contains


  subroutine Initialize_H ( S, Atlas, NameOption, OrderOption )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    class ( Atlas_H_Form ), intent ( in ), target :: &
      Atlas
    character ( * ), intent ( in ), optional :: &
      NameOption
    integer ( KDI ), intent ( in ), optional :: &
      OrderOption

    integer ( KDI ) :: &
      iS  !-- iStage

    S % IGNORABILITY  =  CONSOLE % INFO_1

    S % nStages  =  2
    if ( present ( OrderOption ) ) &
      S % nStages  =  OrderOption
    associate ( nS  =>  S % nStages )

    if ( S % Type == '' ) &
      S % Type = 'a Step_RK' 

    S % Name = 'Step'
    if ( present ( NameOption ) ) &
      S % Name  =  NameOption

    call Show ( 'Initializing ' // trim ( S % Type ), S % IGNORABILITY )
    call Show ( S % Name, 'Name', S % IGNORABILITY )

    call SetCoefficientsExplicit ( S )

    S % Atlas  =>  Atlas

    if ( .not. associated ( S % SetSlope ) ) &
      S % SetSlope  =>  SetSlope_H
    if ( .not. associated ( S % SetSlopeStage ) ) &
      S % SetSlopeStage  =>  SetSlopeStage_H

     call S % SetSlope ( S % Slope )

    allocate ( S % SlopeStage ( nS ) )
    do iS  =  1,  nS
      call S % SetSlopeStage ( S % SlopeStage ( iS ) % Element, iS )
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
      call S % SetSlope ( S % SlopeSum )
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

    call Show ( S % nStages, 'nStages', S % IGNORABILITY )

    do iA  =  2, S % nStages
      write ( Index, fmt = '(i1.1)' ) iA
      call Show ( S % A ( iA ) % Value, 'A ( ' // Index // ' )' )
    end do !-- iA

    call Show ( S % B, 'B' )
    call Show ( S % C, 'C', lRealOption = 2 )

    do iS  =  1, S % nStages
      call S % SlopeStage ( iS ) % Element % Show ( )
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

      !-- Set  Y_I  =  Y
      if ( associated ( T_II ) ) call T_II % Start ( )
      call S % InitializeIntermediate ( iS )
      if ( associated ( T_II ) ) call T_II % Stop ( )

      do iK = 1, iS - 1
        associate ( A  =>  S % A ( iS ) % Value ( iK ) )
        !-- Set Y_I  =  Y_I  +  dT * A * K ( iK )
        if ( associated ( T_II_A ) ) call T_II_A % Start ( )
        if ( A  /=  0.0_KDR ) & 
          call S % IncrementIntermediate ( A, dT, iK )
        if ( associated ( T_II_A ) ) call T_II_A % Stop ( )
        end associate !-- A
        if ( iK  ==  iS - 1 ) then
          if ( associated ( T_SI ) ) call T_SI % Start ( )
          call S % StoreIntermediate ( T_Option = T_SI )
          if ( associated ( T_SI ) ) call T_SI % Stop ( )
        end if
      end do !-- iK

      !-- Compute K ( iS )  =  dY/dT ( Y_I )
      if ( present ( T_Option ) ) then
        T_CS  =>  PROGRAM_HEADER % Timer &
                    ( Handle = S % iTimer_CS, &
                      Name = trim ( S % Name ) // '_CmptStg', &
                      Level = T_Option % Level + 1 )
        call T_CS % Start ( )
        call S % ComputeStage ( T, dT, iS, T_Option = T_CS )
        call T_CS % Stop ( )
      else
        call S % ComputeStage ( T, dT, iS )
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

    if ( allocated ( S % SlopeStage ) ) &
      deallocate ( S % SlopeStage )
    if ( allocated ( S % SlopeSum ) ) &
      deallocate ( S % SlopeSum )
    if ( allocated ( S % Slope ) ) &
      deallocate ( S % Slope )
    if ( allocated ( S % A ) ) &
      deallocate ( S % A )
    if ( allocated ( S % BE ) ) &
      deallocate ( S % BE )
    if ( allocated ( S % B ) ) &
      deallocate ( S % B )
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


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

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


  subroutine ComputeStage ( S, T, dT, iS, T_Option )

    class ( Step_RK_H_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       T, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage
    type ( TimerForm ), intent ( inout ), optional :: &
      T_Option

    call Show ( 'ComputeStage should be overridden', CONSOLE % WARNING )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % WARNING )
    call Show ( 'ComputeStage', 'subroutine', CONSOLE % WARNING )

  end subroutine ComputeStage


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
    case ( 3 )
      allocate ( A ( 3 : 3, 2 : 2 ) )
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
    case ( 4 )
      allocate ( A ( 4 : 4, 3 : 3 ) )
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
    case default
      call Show ( 'RungeKutta order not implemented', CONSOLE % ERROR )
      call Show ( S % nStages, 'Order', CONSOLE % ERROR )
      call Show ( 'Step_RK_H__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetCoefficientsExplicit', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- nStages

    allocate ( S % A ( 2 : nS ) )
    do iS  =  2,  nS
      call S % A ( iS ) % Initialize ( iS - 1 )
      S % A ( iS ) % Value  =  A ( iS, 1 : iS - 1 )
    end do !-- iS

    allocate ( S % B ( nS ) )
    S % B  =  B

    allocate ( S % C ( 2 : nS ) )
    S % C  =  C

    allocate ( S % BE ( nS ) )
    S % BE  =  BE

    end associate !-- nS

  end subroutine SetCoefficientsExplicit


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
