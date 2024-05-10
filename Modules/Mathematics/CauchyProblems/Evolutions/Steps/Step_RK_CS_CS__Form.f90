module Step_RK_CS_CS__Form

  !-- Step_RungeKutta_CurrentSet_CurrentSet_Form

  use Basics
  use Manifolds
  use Fields
  use Slopes
  use Step_RK_H__Form
  use Step_RK_CS__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CS_CS_Form
    integer ( KDI ) :: &
      IMPLICIT_UNSET     = 1, &
      IMPLICIT_POOR      = 2, &
      IMPLICIT_FAIR      = 3, &
      IMPLICIT_GOOD      = 4, &
      IMPLICIT_EXCELLENT = 5
    integer ( KDI ) :: &
      MaxImplicitIterations
    integer ( KDI ), dimension ( : ), allocatable :: &
      nImplicitIterations
    integer ( KDI ), dimension ( :, : ), allocatable :: &
      ImplicitQuality_1, &
      ImplicitQuality_2
    real ( KDR ) :: &
      ImplicitTolerance
    real ( KDR ), dimension ( :, :, : ), allocatable :: &
      ImplicitError_1, &
      ImplicitError_2
    character ( LDL ), dimension ( 5 ) :: &
      QUALITY  =  [ 'UNSET    ', 'POOR     ', 'FAIR     ', 'GOOD     ', &
                    'EXCELLENT' ]
!    type ( CollectiveOperation_I_Form ), allocatable :: &
!      CO_Quality
    type ( FieldSet_BM_Form ), allocatable :: &
      Iteration_1, &
      Iteration_2, &
      IterationPrevious_1, &
      IterationPrevious_2, &
      Residual_1, &
      Residual_2
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

    private :: &
      TestImplicitQuality

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

    integer ( KDI ) :: &
      iS  !-- iStage
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

    S % MaxImplicitIterations  =  20
    call PROGRAM_HEADER % GetParameter &
           ( S % MaxImplicitIterations, 'MaxImplicitIterations' )

    S % ImplicitTolerance  =  1.0e-6_KDR
    call PROGRAM_HEADER % GetParameter &
           ( S % ImplicitTolerance, 'ImplicitTolerance' )

    associate &
      ( nS    =>  S % nStages, &
        mII   =>  S % MaxImplicitIterations, &
        nB_1  =>  CS_1 % nBalanced, &
        nB_2  =>  CS_2 % nBalanced )

    allocate ( S % nImplicitIterations ( 2 : nS ) )
    allocate ( S % ImplicitQuality_1 ( nB_1, 2 : nS ) )
    allocate ( S % ImplicitQuality_2 ( nB_2, 2 : nS ) )
    allocate ( S % ImplicitError_1 ( mII, nB_1, 2 : nS ) )
    allocate ( S % ImplicitError_2 ( mII, nB_2, 2 : nS ) )

    S % ImplicitQuality_1  =  S % IMPLICIT_UNSET
    S % ImplicitQuality_2  =  S % IMPLICIT_UNSET

    ! !-- FIXME: Assumes single chart
    ! select type ( A  =>  S % Atlas )
    ! class is ( Atlas_SCG_Form )
    !   allocate ( S % CO_Quality )
    !   call S % CO_Quality % Initialize &
    !          ( A % Chart_GS % Communicator, &
    !            nOutgoing = [ ( nB_1 + nB_2 ) * ( nS - 1 ) ], &
    !            nIncoming = [ ( nB_1 + nB_2 ) * ( nS - 1 ) ] )
    ! class default
    !   call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
    !   call Show ( 'Step_RK_CS_CS_Form', 'module', CONSOLE % ERROR )
    !   call Show ( 'Initialize_CS_CS', 'subroutine', CONSOLE % ERROR )
    !   call PROGRAM_HEADER % Abort ( )
    ! end select !-- A

    end associate !-- nS

    !-- Iterated field values
    allocate ( S % Iteration_1 )
    allocate ( S % Iteration_2 )
    associate &
      ( Y_1  =>  S % Iteration_1, &
        Y_2  =>  S % Iteration_2 )
    call Y_1 % Initialize &
           ( CS_1 % Atlas, &
             FieldOption = CS_1 % Balanced, &
             NameOption = trim ( CS_1 % Name ) // '_Implicit', &
             DeviceMemoryOption = CS_1 % DeviceMemory, &
             DevicesCommunicateOption = CS_1 % DevicesCommunicate, &
             nFieldsOption = CS_1 % nBalanced, &
             IgnorabilityOption = CS_1 % IGNORABILITY + 1 )
    call Y_2 % Initialize &
           ( CS_2 % Atlas, &
             FieldOption = CS_2 % Balanced, &
             NameOption = trim ( CS_2 % Name ) // '_Implicit', &
             DeviceMemoryOption = CS_2 % DeviceMemory, &
             DevicesCommunicateOption = CS_2 % DevicesCommunicate, &
             nFieldsOption = CS_2 % nBalanced, &
             IgnorabilityOption = CS_2 % IGNORABILITY + 1 )
    end associate !-- Y_1, etc.

    !-- Implicit slope computed in the previous iteration 
    allocate ( S % IterationPrevious_1 )
    allocate ( S % IterationPrevious_2 )
    associate &
      ( KK_P_1  =>  S % IterationPrevious_1, &
        KK_P_2  =>  S % IterationPrevious_2 )
    call KK_P_1 % Initialize &
           ( CS_1 % Atlas, &
             FieldOption = CS_1 % Balanced, &
             NameOption = trim ( CS_1 % Name ) // '_IterationPrevious', &
             DeviceMemoryOption = CS_1 % DeviceMemory, &
             DevicesCommunicateOption = CS_1 % DevicesCommunicate, &
             nFieldsOption = CS_1 % nBalanced, &
             IgnorabilityOption = CS_1 % IGNORABILITY + 1 )
    call KK_P_2 % Initialize &
           ( CS_2 % Atlas, &
             FieldOption = CS_2 % Balanced, &
             NameOption = trim ( CS_2 % Name ) // '_IterationPrevious', &
             DeviceMemoryOption = CS_2 % DeviceMemory, &
             DevicesCommunicateOption = CS_2 % DevicesCommunicate, &
             nFieldsOption = CS_2 % nBalanced, &
             IgnorabilityOption = CS_2 % IGNORABILITY + 1 )
    end associate !-- KK_P_1, etc.

    !-- Relative difference in implicit slopes in this and previous iterations
    allocate ( S % Residual_1 )
    allocate ( S % Residual_2 )
    associate &
      ( R_1  =>  S % Residual_1, &
        R_2  =>  S % Residual_2 )
    call R_1 % Initialize &
           ( CS_1 % Atlas, &
             FieldOption = CS_1 % Balanced, &
             NameOption = trim ( CS_1 % Name ) // '_Residual', &
             DeviceMemoryOption = CS_1 % DeviceMemory, &
             DevicesCommunicateOption = CS_1 % DevicesCommunicate, &
             nFieldsOption = CS_1 % nBalanced, &
             IgnorabilityOption = CS_1 % IGNORABILITY + 1 )
    call R_2 % Initialize &
           ( CS_2 % Atlas, &
             FieldOption = CS_2 % Balanced, &
             NameOption = trim ( CS_2 % Name ) // '_Residual', &
             DeviceMemoryOption = CS_2 % DeviceMemory, &
             DevicesCommunicateOption = CS_2 % DevicesCommunicate, &
             nFieldsOption = CS_2 % nBalanced, &
             IgnorabilityOption = CS_2 % IGNORABILITY + 1 )
    end associate !-- RD_1, etc.

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

    call Show ( S % MaxImplicitIterations, 'MaxImplicitIterations', &
                S % IGNORABILITY )
    call Show ( S % ImplicitTolerance, 'ImplicitTolerance', &
                S % IGNORABILITY )

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
    if ( allocated ( S % Residual_2 ) ) &
      deallocate ( S % Residual_2 )
    if ( allocated ( S % Residual_1 ) ) &
      deallocate ( S % Residual_1 )
    if ( allocated ( S % IterationPrevious_2 ) ) &
      deallocate ( S % IterationPrevious_2 )
    if ( allocated ( S % IterationPrevious_1 ) ) &
      deallocate ( S % IterationPrevious_1 )
    if ( allocated ( S % Iteration_2 ) ) &
      deallocate ( S % Iteration_2 )
    if ( allocated ( S % Iteration_1 ) ) &
      deallocate ( S % Iteration_1 )
!    if ( allocated ( S % CO_Quality ) ) &
!      deallocate ( S % CO_Quality )
    if ( allocated ( S % ImplicitError_2 ) ) &
      deallocate ( S % ImplicitError_2 )
    if ( allocated ( S % ImplicitError_1 ) ) &
      deallocate ( S % ImplicitError_1 )
    if ( allocated ( S % ImplicitQuality_2 ) ) &
      deallocate ( S % ImplicitQuality_2 )
    if ( allocated ( S % ImplicitQuality_1 ) ) &
      deallocate ( S % ImplicitQuality_1 )
    if ( allocated ( S % nImplicitIterations ) ) &
      deallocate ( S % nImplicitIterations )

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
  iV
logical ( KDL ), dimension ( : ), allocatable :: &
  Freeze

    if ( iS  ==  1 ) &
      return

    associate &
      ( S_1  =>  S % Step_CS_1, &
        S_2  =>  S % Step_CS_2 )
    associate &
      ( AA      =>  S_1 % AA ( iS ) % Value ( iS ), &          
        Y_I_1   =>  S_1 % Intermediate, &
        Y_I_2   =>  S_2 % Intermediate, &
        CS_B_1  =>  S_1 % Balanced, &
        CS_B_2  =>  S_2 % Balanced, &
        CS_1    =>  S_1 % CurrentSet, &
        CS_2    =>  S_2 % CurrentSet, &
        KK_I_1  =>  S_1 % SlopeImplicitIterate, &
        KK_I_2  =>  S_2 % SlopeImplicitIterate, &
        KK_1    =>  S_1 % SlopeImplicit, &
        KK_2    =>  S_2 % SlopeImplicit, &
        KK_1_S  =>  S_1 % SlopeStageImplicit ( iS ) % Element, &
        KK_2_S  =>  S_2 % SlopeStageImplicit ( iS ) % Element, &
        Y_1     =>  S % Iteration_1, &
        Y_2     =>  S % Iteration_2, &
        Y_P_1   =>  S % IterationPrevious_1, &
        Y_P_2   =>  S % IterationPrevious_2, &
         iII    =>  S % nImplicitIterations ( iS ) )

    if ( AA == 0.0_KDR ) &
      return

    !-- Upon entry, Y_I = Q_(I-1)
    call Y_I_1 % Copy ( Y_1 )
    call Y_I_2 % Copy ( Y_2 )

    !-- Iterate on subset of balanced fields

! call Show ( '>>> Stage' )
! call Show ( iS, '>>> iS' )

! associate &
!   ( Q_R_E    => Y_I_1 % Storage_GS % Value ( :, 1 ), &
!     Q_F_E    => Y_I_2 % Storage_GS % Value ( :, 5 ), &
!     Q_R_N    => Y_I_1 % Storage_GS % Value ( :, 5 ), &
!     Q_F_N    => Y_I_2 % Storage_GS % Value ( :, 6 ), &
!     J_Eq     => CS_1 % Storage_GS % Value ( :, 2 ) )
! !do iV  =  1, 10
!   iV  =  5
!   call Show ( iV, '>>> iV' )
!   call Show ( J_Eq ( iV ),  '>>> J_Eq' )
!   call Show ( Q_R_E ( iV ), '>>> Q_R_E' )
!   call Show ( Q_F_E ( iV ), '>>> Q_F_E' )
!   call Show ( Q_R_N ( iV ), '>>> Q_R_N' )
!   call Show ( Q_F_N ( iV ), '>>> Q_F_N' )
! !end do
! end associate

    iII     =  0
    allocate ( Freeze ( CS_1 % Storage_GS % nValues ) )
    Freeze  =  .false.
    do 

      iII  =  iII + 1
!call Show ( '>>> Iteration' )
!call Show ( iII, '>>> iII' )

!      call KK_1 % Compute ( dT )!, T_Option = T_CS )
!      call KK_2 % Compute ( dT )!, T_Option = T_CS )
      call KK_I_1 % Compute ( dT )!, T_Option = T_CS )
      call KK_I_2 % Compute ( dT )!, T_Option = T_CS )

! associate &
!   ( J_Eq     => CS_1 % Storage_GS % Value ( :, 2 ) )
! call Show ( J_Eq ( iV ),  '>>> J_Eq' )
! end associate

! associate &
!   ( Q_R_E    => Y_I_1 % Storage_GS % Value ( :, 1 ), &
!     Q_R_N    => Y_I_1 % Storage_GS % Value ( :, 5 ), &
!     KK_R_E  => KK_I_1 % Storage_GS % Value ( :, 1 ), &
!     KK_R_N  => KK_I_1 % Storage_GS % Value ( :, 5 ), &
!     KK_F_E  => KK_I_2 % Storage_GS % Value ( :, 5 ), &
!     KK_F_N  => KK_I_2 % Storage_GS % Value ( :, 6 ) )
! call Show ( dT * KK_R_E ( iV ), '>>> dT * KK_R_E raw' )
! call Show ( dT * KK_F_E ( iV ), '>>> dT * KK_F_E raw' )
! call Show ( dT * KK_R_N ( iV ), '>>> dT * KK_R_N raw' )
! call Show ( dT * KK_F_N ( iV ), '>>> dT * KK_F_N raw' )
! where ( abs ( dT * KK_R_E )  >  0.1 * Q_R_E )
!   KK_R_E  =  sign ( 0.1 * Q_R_E / dT, KK_R_E )
!   KK_F_E  =  - KK_R_E
! end where
! where ( abs ( dT * KK_R_N )  >  0.1 * Q_R_N )
!   KK_R_N  =  sign ( 0.1 * Q_R_N / dT, KK_R_N )
!   KK_F_N  =  - KK_R_N
! end where
! call Show ( dT * KK_R_E ( iV ), '>>> dT * KK_R_E edit' )
! call Show ( dT * KK_F_E ( iV ), '>>> dT * KK_F_E edit' )
! call Show ( dT * KK_R_N ( iV ), '>>> dT * KK_R_N edit' )
! call Show ( dT * KK_F_N ( iV ), '>>> dT * KK_F_N edit' )
! end associate

! associate &
!   ( Q_R_E    => Y_I_1 % Storage_GS % Value ( :, 1 ), &
!     Q_R_N    => Y_I_1 % Storage_GS % Value ( :, 5 ), &
!     Q_F_E    => Y_I_2 % Storage_GS % Value ( :, 5 ), &
!     Q_F_N    => Y_I_2 % Storage_GS % Value ( :, 6 ), &
!     Y_R_E    => Y_1 % Storage_GS % Value ( :, 1 ), &
!     Y_R_N    => Y_1 % Storage_GS % Value ( :, 5 ), &
!     Y_F_E    => Y_2 % Storage_GS % Value ( :, 5 ), &
!     Y_F_N    => Y_2 % Storage_GS % Value ( :, 6 ), &
!     KK_R_E  => KK_I_1 % Storage_GS % Value ( :, 1 ), &
!     KK_R_N  => KK_I_1 % Storage_GS % Value ( :, 5 ), &
!     KK_F_E  => KK_I_2 % Storage_GS % Value ( :, 5 ), &
!     KK_F_N  => KK_I_2 % Storage_GS % Value ( :, 6 ) )
! call Show ( Q_R_E ( iV )  +  dT * KK_R_E ( iV ), '>>> Q_R_E + dT * KK_R_E' )
! call Show ( Q_F_E ( iV )  +  dT * KK_F_E ( iV ), '>>> Q_F_E + dT * KK_F_E' )
! call Show ( Q_R_N ( iV )  +  dT * KK_R_N ( iV ), '>>> Q_R_N + dT * KK_R_N' )
! call Show ( Q_F_N ( iV )  +  dT * KK_F_N ( iV ), '>>> Q_F_N + dT * KK_F_N' )
! call Show ( Y_R_E ( iV )  +  dT * KK_R_E ( iV ), '>>> Y_R_E + dT * KK_R_E' )
! call Show ( Y_F_E ( iV )  +  dT * KK_F_E ( iV ), '>>> Y_F_E + dT * KK_F_E' )
! call Show ( Y_R_N ( iV )  +  dT * KK_R_N ( iV ), '>>> Y_R_N + dT * KK_R_N' )
! call Show ( Y_F_N ( iV )  +  dT * KK_F_N ( iV ), '>>> Y_F_N + dT * KK_F_N' )
! where ( .not. Freeze &
!         .and. ( Q_R_E  +  dT * AA * KK_R_E  >  0.0_KDR &
!                 .or. Q_R_N  +  dT * AA * KK_R_N  >  0.0_KDR ) )
!   Y_R_E  =  Q_R_E  +  dT * AA * KK_R_E
!   Y_F_E  =  Q_F_E  +  dT * AA * KK_F_E
!   Y_R_N  =  Q_R_N  +  dT * AA * KK_R_N
!   Y_F_N  =  Q_F_N  +  dT * AA * KK_F_N
! else where ( .not. Freeze )
!   Y_R_E  =  Y_R_E  +  dT * AA * KK_R_E
!   Y_F_E  =  Y_F_E  +  dT * AA * KK_F_E
!   Y_R_N  =  Y_R_N  +  dT * AA * KK_R_N
!   Y_F_N  =  Y_F_N  +  dT * AA * KK_F_N
!   Freeze  =  .true.
! end where
! end associate

    call Y_1 % MultiplyAdd ( Y_I_1, KK_I_1, dT * AA )      
    call Y_2 % MultiplyAdd ( Y_I_2, KK_I_2, dT * AA )      

! associate &
!   ( Y_R_E    => Y_1   % Storage_GS % Value ( :, 1 ), &
!     Y_F_E    => Y_2   % Storage_GS % Value ( :, 5 ), &
!     Y_R_N    => Y_1   % Storage_GS % Value ( :, 5 ), &
!     Y_F_N    => Y_2   % Storage_GS % Value ( :, 6 ) )
! !do iV  =  1, 10
!   call Show ( Y_R_E ( iV ), '>>> Y_R_E' )
!   call Show ( Y_F_E ( iV ), '>>> Y_F_E' )
!   call Show ( Y_R_N ( iV ), '>>> Y_R_N' )
!   call Show ( Y_F_N ( iV ), '>>> Y_F_N' )
! !end do
! end associate

! ! !-- Where negative, reset
! ! if ( any ( Y_R_E  <  0.0_KDR )  .or.  any ( Y_R_N  <  0.0_KDR ) ) then
! ! call Show ( '>>> Negative density encountered' )
! ! where ( Y_R_E  <  0.0_KDR )
! !   Y_R_E  =  Q_R_E  +  0.02 * iII * Q_R_E
! !   Y_F_E  =  Q_F_E  -  0.02 * iII * Q_R_E
! ! end where
! ! where ( Y_R_N  <  0.0_KDR )
! !   Y_R_N  =  Q_R_N  +  0.02 * iII * Q_R_N
! !   Y_F_N  =  Q_F_N  -  0.02 * iII * Q_R_N
! ! end where
! ! !do iV  =  1, 10
! !   iV  =  3
! !   call Show ( iV, '>>> iV' )
! !   call Show ( Y_R_E ( iV ), '>>> Y_R_E edited' )
! !   call Show ( Y_F_E ( iV ), '>>> Y_F_E edited' )
! !   call Show ( Y_R_N ( iV ), '>>> Y_R_N edited' )
! !   call Show ( Y_F_N ( iV ), '>>> Y_F_N edited' )
! ! !end do
! ! end if

! end associate

      !-- Fill out fields for next iteration of subset of balanced fields, 
      !   or update of all balanced fields after exit

      call Y_1 % Copy ( CS_B_1 )
      call CS_1 % ComputeFromBalanced ( )

      call Y_2 % Copy ( CS_B_2 )
      call CS_2 % ComputeFromBalanced ( )

! call Show ( iII, '>>> iII' )
! call Show ( Y_1 % Storage ( 1 ) % Value ( 3, 1 ), '>>> Y_1' )
! call Show ( Y_2 % Storage ( 1 ) % Value ( 3, 5 ), '>>> Y_2' )

      if ( iII  >  1 ) then

        !-- Exit conditions

        associate &
          ( R_1   =>  S % Residual_1, &
            R_2   =>  S % Residual_2, &
            IE_1  =>  S % ImplicitError_1 ( :, :, iS ), &
            IE_2  =>  S % ImplicitError_2 ( :, :, iS ), &
            IQ_1  =>  S % ImplicitQuality_1 ( :, iS ), &
            IQ_2  =>  S % ImplicitQuality_2 ( :, iS ) )

        call TestImplicitQuality ( S, R_1, IE_1, IQ_1, Y_1, Y_P_1, Y_I_1, iII )
        call TestImplicitQuality ( S, R_2, IE_2, IQ_2, Y_2, Y_P_2, Y_I_2, iII )

        if (     any ( IQ_1  ==  S % IMPLICIT_POOR ) &
            .or. any ( IQ_2  ==  S % IMPLICIT_POOR ) ) &
            ! any ( IQ_2  ==  S % IMPLICIT_POOR ) ) &
        then
call Show ( '>>> Exit diverging', CONSOLE % WARNING )
call Show ( iS, '>>> iS' )
call Show ( iE_1 ( : iII, : ), '>>> iE_1', CONSOLE % WARNING )
call Show ( iE_2 ( : iII, : ), '>>> iE_2', CONSOLE % WARNING )
call PROGRAM_HEADER % Abort ( )
! call Show ( '>>> Reducing Limiter in implicit solver', CONSOLE % WARNING )

!           !-- Reset to Q ( i-1 ) and start iteration over with smaller dT_IS,
!           !   effectively weakening the right-hand side

!           call Y_I_1 % Copy ( CS_B_1 )
!           call CS_1 % ComputeFromBalanced ( )

!           call Y_I_2 % Copy ( CS_B_2 )
!           call CS_2 % ComputeFromBalanced ( )

!           iII  =  0
!           Limiter  =  Limiter  /  2.0_KDR

! !call Show ( Limiter, '>>> Limiter (implicit solve)' )

!           cycle

! associate &
!   ( Q_1_E    => Y_I_1 % Storage_GS % Value ( :, 1 ), &
!     Y_P_1_E  => Y_P_1 % Storage_GS % Value ( :, 1 ), &
!     Y_1_E    => Y_1   % Storage_GS % Value ( :, 1 ), &
!     Q_2_E    => Y_I_2 % Storage_GS % Value ( :, 5 ), &
!     Y_P_2_E  => Y_P_2 % Storage_GS % Value ( :, 5 ), &
!     Y_2_E    => Y_2   % Storage_GS % Value ( :, 5 ), &
!     Q_1_N    => Y_I_1 % Storage_GS % Value ( :, 5 ), &
!     Y_P_1_N  => Y_P_1 % Storage_GS % Value ( :, 5 ), &
!     Y_1_N    => Y_1   % Storage_GS % Value ( :, 5 ), &
!     Q_2_N    => Y_I_2 % Storage_GS % Value ( :, 6 ), &
!     Y_P_2_N  => Y_P_2 % Storage_GS % Value ( :, 6 ), &
!     Y_2_N    => Y_2   % Storage_GS % Value ( :, 6 ) )
! do iV  =  1, 20
!   call Show ( iV, '>>> iV' )
!   call Show ( [ Q_1_E ( iV ), Y_P_1_E ( iV ), Y_1_E ( iV ) ], '>>> E_R' )
!   call Show ( [ Q_2_E ( iV ), Y_P_2_E ( iV ), Y_2_E ( iV ) ], '>>> E_F' )
!   call Show ( [ Q_1_N ( iV ), Y_P_1_N ( iV ), Y_1_N ( iV ) ], '>>> N_R' )
!   call Show ( [ Q_2_N ( iV ), Y_P_2_N ( iV ), Y_2_N ( iV ) ], '>>> N_F' )
! end do
! end associate
           exit  !-- Diverging
        end if

        if (      all ( IQ_1  /=  S % IMPLICIT_UNSET ) &
            .and. all ( IQ_2  /=  S % IMPLICIT_UNSET ) ) &
            !all ( IQ_2  /=  S % IMPLICIT_UNSET ) ) &
        then
!call Show ( '>>> Exit converging' )
!           if (     any ( IQ_1  /=  S % IMPLICIT_EXCELLENT ) &
!               .or. any ( IQ_2  /=  S % IMPLICIT_EXCELLENT ) ) &
!           then
! call Show ( '>>> Exit slowly converging', CONSOLE % WARNING )
! call Show ( iS, '>>> iS' )
! call Show ( iE_1 ( : iII, : ), '>>> iE_1', CONSOLE % WARNING )
! call Show ( iE_2 ( : iII, : ), '>>> iE_2', CONSOLE % WARNING )
!           end if
          exit  !-- Converged, or converging
        end if

        end associate !-- RD_1, etc.

      end if !-- iII > 1

      ! !-- Set up next iteration

      ! call Y_1 % Copy ( Y_P_1 )
      ! call Y_2 % Copy ( Y_P_2 )

      ! call Y_1 % Copy ( CS_B_1 )
      ! call CS_1 % ComputeFromBalanced ( )

      ! call Y_2 % Copy ( CS_B_2 )
      ! call CS_2 % ComputeFromBalanced ( )

      !-- Set up for exit test next iteration

      call Y_1 % Copy ( Y_P_1 )
      call Y_2 % Copy ( Y_P_2 )

    end do !-- iII

    !-- After convergence of subset, compute update of all fields

    call KK_1 % Compute ( dT )!, T_Option = T_CS )
    call KK_2 % Compute ( dT )!, T_Option = T_CS )

! associate &
!   ( KK_R_V  => KK_1 % Storage_GS % Value, &
!     KK_F_V  => KK_2 % Storage_GS % Value )
! KK_R_V  =  Limiter * KK_R_V
! KK_R_V  =  Limiter * KK_R_V
! end associate

    !-- Assume SlopeImplicit local: no ghost exchange in loop above 
    call KK_1 % ExchangeGhostData ( )
    call KK_2 % ExchangeGhostData ( )

    !-- Copy to stage storage
    call KK_1 % Copy ( KK_1_S )
    call KK_2 % Copy ( KK_2_S )

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

    integer ( KDI ) :: &
      iS, &
      iF

    call S % Step_CS_1 % StoreSolution ( T_Option )
    call S % Step_CS_2 % StoreSolution ( T_Option )

!     !-- Reduce ImplicitQuality over base manifold domains
!     associate ( CO => S % CO_Quality )
!     CO % Outgoing % Value  =  S % ImplicitQuality
!     call CO % Reduce ( REDUCTION % MAX )
!     S % ImplicitQuality  =  CO % Incoming % Value
!     end associate !-- CO

! do iS  =  2, S % nStages
!   call Show ( iS, '>>> iS' )
!   associate ( nII  =>  S % nImplicitIterations ( iS ) )
!   call Show ( nII, '>>> nImplicitIterations' )
!   call Show ( S % ImplicitError_1 ( nII, :, iS ), '>>> ImplicitError_1' )
!   call Show ( S % ImplicitError_2 ( nII, :, iS ), '>>> ImplicitError_2' )
!   end associate !-- nII
!   call Show ( [ ( S % QUALITY ( S % ImplicitQuality_1 ( iF, iS ) ), &
!                   iF = 1, size ( S % ImplicitQuality_1, dim = 1 ) ) ], &
!                 '>>> ImplicitQuality_1' )
!   call Show ( [ ( S % QUALITY ( S % ImplicitQuality_2 ( iF, iS ) ), &
!                   iF = 1, size ( S % ImplicitQuality_2, dim = 1 ) ) ], &
!                 '>>> ImplicitQuality_2' )
! end do

  end subroutine StoreSolution


  subroutine TestImplicitQuality ( S, R, IE, IQ, Y, Y_P, Y_I, iII )

    class ( Step_RK_CS_CS_Form ), intent ( in ) :: &
      S
    class ( FieldSet_BM_Form ), intent ( inout ) :: &
      R  !-- Residual
    real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      IE  !-- ImplicitError
    integer ( KDI ), dimension ( : ), intent ( inout ) :: &
      IQ  !-- ImplicitQuality
    class ( FieldSet_BM_Form ), intent ( in ) :: &
      Y, &    !-- current iteration
      Y_P, &  !-- previous iteration
      Y_I     !-- upon entry to implicit solver
    integer ( KDI ), intent ( in ) :: &
      iII  !-- iImplicitIterations

    integer ( KDI ) :: &
      iF

    !-- First tested iteration
    if ( iII  ==  2 ) then
      IE  =  0.0_KDR
      IQ  =  S % IMPLICIT_UNSET
    end if

    call R % RelativeDifference ( Y, Y_P, Y_I )

    do iF  =  1, Y % nFields

      associate &
        ( IEV   =>  IE ( iII,     iF ), &
          IEPV  =>  IE ( iII - 1, iF ), &
          IQV   =>  IQ ( iF ) )

      !-- FIXME: Assumes single chart
      select type ( A  =>  S % Atlas )
      class is ( Atlas_SCG_Form )
        IEV  =  maxval ( R % Storage_GS % Value ( :, iF ) )
! if ( iF == 1 ) then
!   call Show ( iF, '>>> Test exit iF' )
!   call Show ( maxloc ( R % Storage_GS % Value ( :, iF ) ), '>>> maxloc' )
!   call Show ( IEV, '>>> maxval' )
! end if
      class default
        call Show ( 'Atlas type not recognized', CONSOLE % ERROR )
        call Show ( 'Step_RK_CS_CS_Form', 'module', CONSOLE % ERROR )
        call Show ( 'TestImplicitQuality', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select !-- A

! !-- For testing, exclude momentum
! if ( iF == 2 .or. iF == 3 .or. iF == 4 ) then
!   IQV  =  S % IMPLICIT_EXCELLENT
!   cycle
! end if

      if ( IQV  /=  S % IMPLICIT_UNSET )  &
        cycle

      if ( IEV  <  S % ImplicitTolerance ) then
        !-- Converged
        if ( iII  <  S % MaxImplicitIterations / 2 ) then
          IQV  =  S % IMPLICIT_EXCELLENT
        else
          IQV  =  S % IMPLICIT_GOOD
        end if
!      else if ( iII  >  2 .and. IEV  >  IEPV ) then  
      else if ( iII  >  2 .and. IEV  >  IEPV .and. IEV  >  0.01_KDR ) then  
!      else if ( IEV  >  IEPV .and. iII  ==  S % MaxImplicitIterations ) then  
!      else if ( iII  ==  S % MaxImplicitIterations ) then  
        !-- Diverging
        IQV  =  S % IMPLICIT_POOR
      else if ( iII  ==  S % MaxImplicitIterations ) then  
        !-- Converging slowly
        IQV  =  S % IMPLICIT_FAIR
      end if

      end associate !-- IEV, etc.

    end do !-- iF

  end subroutine TestImplicitQuality


end module Step_RK_CS_CS__Form
