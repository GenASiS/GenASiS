module Step_RK_CSA__Form

  !-- Step_RungeKutta_FieldSetAtlas_Form

  use Basics
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
    class ( FieldSet_A_Form ), allocatable :: &
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
      InitializeIntermediate
    procedure, private, pass :: &
      IncrementIntermediate
    procedure, public, nopass :: &
      LoadSolution_C
    procedure, public, nopass :: &
      InitializeIntermediate_C
    procedure, public, nopass :: &
      IncrementIntermediate_C
  end type Step_RK_CSA_Form


contains


  subroutine Initialize_CSA ( S, CSA, NameOption )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    class ( CurrentSet_A_Form ), intent ( in ), target :: &
      CSA
    character ( * ), intent ( in ), optional :: &
      NameOption

    real ( KDR ), dimension ( 2 : 2, 1 : 1 ) :: &
      A
    real ( KDR ), dimension ( 2 : 2 ) :: &
      C
    real ( KDR ), dimension ( 1 : 2 ) :: &
      B

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

    !-- Slope

    if ( .not. allocated ( S % Slope_A ) ) then
      allocate ( Slope_DFV_A_Form :: S % Slope_A )
      select type ( SA  =>  S % Slope_A )
      class is ( Slope_DFV_A_Form )
      call SA % Initialize ( CSA )
      end select !-- SA
    end if !-- allocated Slope

    !-- Cleanup

    end associate !-- nEquations, etc.
    end select !-- CSC

  end subroutine Initialize_CSA


  subroutine Show_S ( S )

    class ( Step_RK_CSA_Form ), intent ( in ) :: &
      S

    call S % Step_RK_H_Form % Show ( )
    call S % Solution_A % Show ( )
    call S % Intermediate_A % Show ( )
    call S % Slope_A % Show ( )

  end subroutine Show_S


  impure elemental subroutine Finalize ( S )

    type ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S

    if ( allocated ( S % Slope_A ) ) &
      deallocate ( S % Slope_A )

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


  subroutine InitializeIntermediate ( S, iS )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      iS  !-- iStage

    integer ( KDI ) :: &
      iC  !-- iChart

    if ( iS  ==  1 ) &
      return

    associate ( nC  =>  S % CurrentSet_A % Atlas % nCharts )
    do iC  =  1,  nC

      associate &
        ( Intermediate_C  =>  S % Solution_A % FieldSet_C ( iC ) % Element, &
              Solution_C  =>  S % Solution_A % FieldSet_C ( iC ) % Element )

      call S % InitializeIntermediate_C ( Intermediate_C, Solution_C )

      end associate !-- Intermediate_C

    end do !-- iC
    end associate !-- nC

  end subroutine InitializeIntermediate


  subroutine IncrementIntermediate ( S, A, iK )

    class ( Step_RK_CSA_Form ), intent ( inout ) :: &
      S
    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    call Show ( 'IncrementIntermediate must be overridden', CONSOLE % ERROR )
    call Show ( 'Step_RK_H_Form', 'module', CONSOLE % ERROR )
    call Show ( 'IncrementIntermediate', 'subroutine', CONSOLE % ERROR )
    call PROGRAM_HEADER % Abort ( )

  end subroutine IncrementIntermediate


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


  subroutine IncrementIntermediate_C ( A, iK )

    real ( KDR ), intent ( in ) :: &
      A
    integer ( KDI ), intent ( in ) :: &
      iK

    associate &
      ( YV  => S % Y % Value, &
        KV  => S % K ( iK ) % Field % Value )
    
    call MultiplyAdd ( YV, KV, A, UseDeviceOption = S % Y % AllocatedDevice )
    
    end associate !-- YV, etc.

  end subroutine IncrementIntermediate_C


end module Step_RK_CSA__Form
