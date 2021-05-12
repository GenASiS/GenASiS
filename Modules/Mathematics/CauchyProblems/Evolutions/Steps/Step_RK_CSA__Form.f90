module Step_RK_CSA__Form

  !-- Step_RungeKutta_FieldSetAtlas_Form

  use Basics
  use Algebra
  use Fields
  use Slopes
  use Step_RK_H__Form

  implicit none
  private

  type, public, extends ( Step_RK_H_Form ) :: Step_RK_CSA_Form
    type ( FieldSet_A_Form ), allocatable :: &
      Solution_A, &
      Intermediate_A
    class ( CurrentSet_A_Form ), pointer :: &
      CurrentSet_A
    class ( RiemannSolver_HLL_A_Form ), allocatable :: &
      RiemannSolver_A
    type ( FieldSet_A_Element ), dimension ( : ), allocatable :: &
      Slope_A  !-- Some extension of FieldSet_A_Form that computes slopes
  contains
    procedure, private, pass :: &
      Initialize_CSA
    generic, public :: &
      Initialize => Initialize_CSA
    procedure, private, pass :: &
      Show_S
    final :: &
      Finalize
    procedure, private, pass :: &
      LoadSolution
    procedure, private, pass :: &
      StoreSolution
    procedure, private, pass :: &
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, private, pass :: &
      ComputeStage
    procedure, private, pass :: &
      IncrementSolution
    procedure, public, nopass :: &
      LoadSolution_C
    procedure, public, nopass :: &
      StoreSolution_C
    procedure, public, nopass :: &
      InitializeIntermediate_C
    procedure, public, nopass :: &
      IncrementIntermediate_C
    procedure, public, nopass :: &
      ComputeStage_C
    procedure, public, nopass :: &
      IncrementSolution_C
  end type Step_RK_CSA_Form


contains


  subroutine Initialize_CSA ( S, CSA, NameOption )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    integer ( KDI ) :: &
      iS  !-- iStage
    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B
    character ( 1 ) :: &
      StageNumber

    if ( S % Type  ==  '' ) &
      S % Type  =  'a Step_RK_CSA'

    S % CurrentSet_A  =>  CSA

    call Clear ( A )
    A ( 2, 1 ) = 1.0_KDR

    B ( 1 ) = 0.5_KDR
    B ( 2 ) = 0.5_KDR

    C ( 2 ) = 1.0_KDR
    
    select type ( CSC  =>  CSA % FieldSet_C ( 1 ) % Element )
    class is ( CurrentSet_C_Form )

    associate &
      (         nEquations  =>  CSC % nBalanced, &
                  Equation  =>  CSC % Balanced, & 
              DeviceMemory  =>  CSC % Storage_FSC % DeviceMemory, &
              PinnedMemory  =>  CSC % Storage_FSC % PinnedMemory, &
        DevicesCommunicate  =>  CSC % GhostExchange_FSC % DevicesCommunicate )

    call S % Step_RK_H_Form % Initialize &
           ( A, B, C, nEquations, NameOption )

    !-- Solution storage

    allocate ( S % Solution_A )
    associate ( SA  =>  S % Solution_A )
    call SA % Initialize &
           ( CSA % Atlas, &
             FieldOption = Equation, &
             NameOption = 'Solution', &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nEquations )
    end associate !-- SA 

    !-- Intermediate storage

    allocate ( S % Intermediate_A )
    associate ( IA  =>  S % Intermediate_A )
    call IA % Initialize &
           ( CSA % Atlas, &
             FieldOption = Equation, &
             NameOption = 'Intermediate', &
             DeviceMemoryOption = DeviceMemory, &
             PinnedMemoryOption = PinnedMemory, &
             DevicesCommunicateOption = DevicesCommunicate, &
             nFieldsOption = nEquations )
    end associate !-- SA 

    !-- RiemannSolver

    if ( .not. allocated ( S % Slope_A ) ) then
      allocate ( S % RiemannSolver_A )
      associate ( RSA  =>  S % RiemannSolver_A )
      call RSA % Initialize ( CSA )
      end associate !-- RSA
    end if !-- allocated RiemannSolver

    !-- Slopes

    if ( .not. allocated ( S % Slope_A ) ) then
      associate ( nS  =>  S % nStages )
      allocate ( S % Slope_A ( nS ) )
      do iS  =  1,  nS
        write ( StageNumber, fmt = '(i1.1)' ) iS
        allocate ( Slope_DFV_A_Form :: S % Slope_A ( iS ) % Element )
        select type ( SA  =>  S % Slope_A ( iS ) % Element )
          class is ( Slope_DFV_A_Form )
        associate ( RSA  =>  S % RiemannSolver_A )
        call SA % Initialize &
               ( RSA, &
                 NameOption = 'SDFV_' // StageNumber // '_' &
                                // trim ( CSA % Name ) )
        end associate !-- RSA
        end select !-- SA
      end do !-- iS
      end associate !-- nS
    end if !-- allocated Slope

    !-- Cleanup

    end associate !-- nEquations, etc.
    end select !-- CSC

  end subroutine Initialize_CSA


  subroutine Show_S ( S )

    class ( Step_RK_CSA_Form ), intent ( in ) :: &
      S

    integer ( KDI ) :: &
      iS  !-- iStage

    call S % Step_RK_H_Form % Show ( )

    call S % Solution_A % Show ( )
    call S % Intermediate_A % Show ( )
    call S % RiemannSolver_A % Show ( )
    do iS  =  1, S % nStages
      call S % Slope_A ( iS ) % Element % Show ( )
    end do !-- iS

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Slope_A ) ) &
      deallocate ( S % Slope_A )
    if ( allocated ( S % RiemannSolver_A ) ) &
      deallocate ( S % RiemannSolver_A )

    nullify ( S % CurrentSet_A )

    if ( allocated ( S % Intermediate_A ) ) &
      deallocate ( S % Intermediate_A )
    if ( allocated ( S % Solution_A ) ) &
      deallocate ( S % Solution_A )

  end subroutine Finalize


  subroutine LoadSolution ( S )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      associate &
        ( Solution_C  =>  S % Solution_A % FieldSet_C ( iC ) % Element )
      select type &
          ( CurrentSet_C  =>  S % CurrentSet_A % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )

      call S % LoadSolution_C ( Solution_C, CurrentSet_C )

      end select !-- CSC
      end associate !-- Solution_C

    end do !-- iC
    end associate !-- nC

  end subroutine LoadSolution


  subroutine StoreSolution ( S )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      select type &
          ( CurrentSet_C  =>  S % CurrentSet_A % FieldSet_C ( iC ) % Element )
        class is ( CurrentSet_C_Form )
      associate &
        ( Solution_C  =>  S % Solution_A % FieldSet_C ( iC ) % Element )

      call S % StoreSolution_C ( CurrentSet_C, Solution_C )

      end associate !-- Solution_C
      end select !-- CSC

    end do !-- iC
    end associate !-- nC

  end subroutine StoreSolution


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    integer ( KDI ) :: &
      iC  !-- iChart

    if ( iS  >  1 ) then

      associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
      do iC  =  1,  nC

        associate &
          ( Intermediate_C  =>  S % Intermediate_A % FieldSet_C ( iC ) &
                                  % Element, &
                Solution_C  =>  S %     Solution_A % FieldSet_C ( iC ) &
                                  % Element )

        call S % InitializeIntermediate_C ( Intermediate_C, Solution_C )

        end associate !-- Intermediate_C, etc.

      end do !-- iC
      end associate !-- nC

    end if !-- iStage > 1

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, dT, iK )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       A, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iK

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      associate &
        ( Intermediate_C  =>  S % Intermediate_A % FieldSet_C ( iC ) &
                                % Element, &
                 Slope_C  =>  S % Slope_A ( iK ) % Element &
                                % FieldSet_C ( iC ) % Element )

      call S % IncrementIntermediate_C ( Intermediate_C, Slope_C, A, dT )

      end associate !-- Intermediate_C, etc.

    end do !-- iC
    end associate !-- nC

  end subroutine IncrementIntermediate


  subroutine ComputeStage ( S, T, iS )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      T
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    integer ( KDI ) :: &
      iC  !-- iChart

    if ( iS  >  1 ) then
      select type ( S )
      type is ( Step_RK_CSA_Form )
        call S % StoreSolution ( )
      class default
        call Show ( 'Wrong Step type', CONSOLE % ERROR )
        call Show ( 'Step_RK_CSA_Form', 'module', CONSOLE % ERROR )
        call Show ( 'ComputeStage', 'subroutine', CONSOLE % ERROR )
        call PROGRAM_HEADER % Abort ( )
      end select
    end if !-- iStage > 1

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      select type &
        ( Slope_C  =>  S % Slope_A ( iS ) % Element &
                         % FieldSet_C ( iC ) % Element )
      class is ( Slope_DFV_C_Form )

      call S % ComputeStage_C ( Slope_C )

      end select !-- Slope_C

    end do !-- iC
    end associate !-- nC

  end subroutine ComputeStage


  subroutine IncrementSolution ( S, B, dT, iS )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
       B, &
      dT
    integer ( KDI ), intent ( in ) :: &
      iS

    integer ( KDI ) :: &
      iC  !-- iChart

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      associate &
        ( Solution_C  =>  S % Solution_A % FieldSet_C ( iC ) % Element, &
             Slope_C  =>  S % Slope_A ( iS ) % Element &
                            % FieldSet_C ( iC ) % Element )

      call S % IncrementSolution_C ( Solution_C, Slope_C, B, dT )

      end associate !-- Solution_C, etc.

    end do !-- iC
    end associate !-- nC

  end subroutine IncrementSolution


  subroutine LoadSolution_C ( Solution_C, CurrentSet_C )

    type ( FieldSet_C_Form ), intent ( inout ) :: &
      Solution_C
    class ( CurrentSet_C_Form ), intent ( in ) :: &
      CurrentSet_C

    integer ( KDI ) :: &
      iB  !-- iBalanced
      
    associate ( iaB  =>  CurrentSet_C % iaBalanced )
    do iB  =  1,  CurrentSet_C % nBalanced
      
      associate &
        ( CV  => CurrentSet_C % Storage_FSC % Storage &
                   % Value ( :, iaB ( iB ) ), &
          SV  => Solution_C % Storage_FSC % Storage &
                   % Value ( :, iB ) )
      
      call Copy ( CV, SV, &
                  UseDeviceOption = CurrentSet_C % Storage_FSC % DeviceMemory )
      
      end associate !-- CV, etc.
      
    end do !-- iB
    end associate !-- iaB

  end subroutine LoadSolution_C


  subroutine StoreSolution_C ( CurrentSet_C, Solution_C )

    class ( CurrentSet_C_Form ), intent ( inout ) :: &
      CurrentSet_C
    type ( FieldSet_C_Form ), intent ( in ) :: &
      Solution_C

    integer ( KDI ) :: &
      iB  !-- iBalanced
      
    associate ( iaB  =>  CurrentSet_C % iaBalanced )
    do iB  =  1,  CurrentSet_C % nBalanced
      
      associate &
        ( SV  => Solution_C % Storage_FSC % Storage &
                   % Value ( :, iB ), &
          CV  => CurrentSet_C % Storage_FSC % Storage &
                   % Value ( :, iaB ( iB ) ) )
      
      call Copy ( SV, CV, &
                  UseDeviceOption = Solution_C % Storage_FSC % DeviceMemory )
      
      end associate !-- CV, etc.
      
    end do !-- iB
    end associate !-- iaB

    call CurrentSet_C % ComputeFromConserved ( )

  end subroutine StoreSolution_C


  subroutine InitializeIntermediate_C ( Intermediate_C, Solution_C )

    type ( FieldSet_C_Form ), intent ( inout ) :: &
      Intermediate_C
    type ( FieldSet_C_Form ), intent ( in ) :: &
      Solution_C

      associate &
        ( SV  => Solution_C % Storage_FSC % Storage % Value, &
          YV  => Intermediate_C % Storage_FSC % Storage % Value )

    call Copy ( SV, YV, &
                UseDeviceOption = Intermediate_C % Storage_FSC % DeviceMemory )

    end associate !-- SV, etc.

  end subroutine InitializeIntermediate_C


  subroutine IncrementIntermediate_C ( Intermediate_C, Slope_C, A, dT )

    type ( FieldSet_C_Form ), intent ( inout ) :: &
      Intermediate_C
    class ( FieldSet_C_Form ), intent ( in ) :: &
      Slope_C
    real ( KDR ), intent ( in ) :: &
       A, &
      dT

    associate &
      ( YV  =>  Intermediate_C % Storage_FSC % Storage % Value, &
        KV  =>  Slope_C % Storage_FSC % Storage % Value )
    
    call MultiplyAdd &
           ( YV, KV, dT * A, &
             UseDeviceOption = Intermediate_C % Storage_FSC % DeviceMemory )
    
    end associate !-- YV, etc.

  end subroutine IncrementIntermediate_C


  subroutine ComputeStage_C ( Slope_C )

    class ( Slope_DFV_C_Form ), intent ( inout ) :: &
      Slope_C

    call Slope_C % Clear ( )
    call Slope_C % Compute ( )
    call Slope_C % ExchangeGhostData ( )

  end subroutine ComputeStage_C


  subroutine IncrementSolution_C ( Solution_C, Slope_C, B, dT )

    type ( FieldSet_C_Form ), intent ( inout ) :: &
      Solution_C
    class ( FieldSet_C_Form ), intent ( in ) :: &
      Slope_C
    real ( KDR ), intent ( in ) :: &
       B, &
      dT

    associate &
      ( SV  => Solution_C % Storage_FSC % Storage % Value, &
        KV  =>    Slope_C % Storage_FSC % Storage % Value )
    
    call MultiplyAdd &
           ( SV, KV, dT * B, &
             UseDeviceOption = Solution_C % Storage_FSC % DeviceMemory )
    
    end associate !-- SV, etc.

  end subroutine IncrementSolution_C


end module Step_RK_CSA__Form
