module RiemannSolver_HLL__Form

  !-- RiemannSolver_HartenLaxVanLeer__Form

  use Basics
  use Fields
  use Reconstruction_Form

  implicit none
  private

    integer ( KDI ), private, parameter :: &
      N_SOLVER_SPEEDS_HLL  =  2

  type, public, extends ( FieldSetForm ) :: RiemannSolver_HLL_Form
    integer ( KDI ) :: &
      N_SOLVER_SPEEDS_HLL = N_SOLVER_SPEEDS_HLL
    integer ( KDI ) :: &
      ALPHA_PLUS_U    = 0, &
      ALPHA_MINUS_U   = 0, &
      N_SOLVER_SPEEDS = 0
    integer ( KDI ) :: &
      iTimer       = 0, &
      iTimerKernel = 0
    class ( FieldSetForm ), allocatable :: &
      BalancedSet
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
    class ( FluxSetForm ), allocatable :: &
      FluxSet
    class ( EigenspeedSet_F_Form ), allocatable :: &
      EigenspeedSet
    class ( ReconstructionForm ), allocatable :: &
      Reconstruction_BS, &
      Reconstruction_FS, &
      Reconstruction_ES
    type ( FieldSetElement ), dimension ( :, : ), allocatable :: &
      StageDimension
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    generic, public :: &
      Initialize => InitializeAllocate_RS
    procedure, public, pass :: &
      SetStream
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Timer
    procedure, private, pass :: &
      TimerKernel
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLL_Form

    private :: &
      ComputeKernel

    interface
      
      module subroutine ComputeKernel &
               ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
                 F_I, AP_I, AM_I, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          F_IL, F_IR, &
          U_IL, U_IR
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          EP_IL, EP_IR, &
          EM_IL, EM_IR
        real ( KDR ), dimension ( :, : ), intent ( out ) :: &
          F_I
        real ( KDR ), dimension ( : ), intent ( out ) :: &
          AP_I, AM_I
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_RS &
               ( RS, CS, FieldOption, PrefixOption, nFieldsOption )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      PrefixOption
    integer ( KDI ), intent ( in ), optional :: &
      nFieldsOption

    integer ( KDI ) :: &
      nFields
    character ( LDL ) :: &
      Name
    character ( LDL ), dimension ( : ), allocatable :: &
      Field

    if ( RS % Type  ==  '' ) &
      RS % Type  =  'a RiemannSolver_HLL' 
    
    Name  =  'RS_' // trim ( CS % Name )
    if ( present ( PrefixOption ) ) &
      Name  =  trim ( PrefixOption ) // '_' // trim ( CS % Name )

    RS % CurrentSet  =>  CS

    allocate &
      ( RS % BalancedSet, &
        RS % FluxSet, &
        RS % EigenspeedSet )
    associate &
      ( BS  =>  RS % BalancedSet, &
        FS  =>  RS % FluxSet, &
        ES  =>  RS % EigenspeedSet )
    call BS % Initialize &
           ( CS, CS % iaBalanced, &
             NameOption = 'B_' // trim ( CS % Name ), &
             IgnorabilityOption = CS % IGNORABILITY + 1 )
    call FS % Initialize ( CS )
    call ES % Initialize ( CS )

    allocate &
      ( RS % Reconstruction_BS, &
        RS % Reconstruction_FS, &
        RS % Reconstruction_ES )
    associate &
      ( RBS  =>  RS % Reconstruction_BS, &
        RFS  =>  RS % Reconstruction_FS, &
        RES  =>  RS % Reconstruction_ES )
    call RBS % Initialize ( CS % Geometry, BS )
    call RFS % Initialize ( CS % Geometry, FS )
    call RES % Initialize ( CS % Geometry, ES )

    associate ( nB  =>  CS % nBalanced )

    !-- Field indices

    if ( present ( nFieldsOption ) ) then
      nFields  =  nFieldsOption
    else
      RS % N_SOLVER_SPEEDS  =  RS % N_SOLVER_SPEEDS_HLL
      nFields  =  nB  +  RS % N_SOLVER_SPEEDS
    end if

    RS % ALPHA_PLUS_U   =  nB  +  1
    RS % ALPHA_MINUS_U  =  nB  +  2

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( : nB )  =  CS % Balanced

    Field ( nB + 1 : nB + RS % N_SOLVER_SPEEDS ) &
      =  [ 'AlphaPlus_U ', &
           'AlphaMinus_U' ]
          
    !-- FieldSet

    call RS % FieldSetForm % Initialize &
           ( CS % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             nFieldsOption = nFields, &
             IgnorabilityOption = CS % IGNORABILITY + 1 )

    end associate !-- nB
    end associate !-- RBS, etc.
    end associate !-- BS, etc.

  end subroutine InitializeAllocate_RS


  subroutine SetStream ( RS, S, nS )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      nS  !-- nStages

    integer ( KDI ) :: &
      iS, &  !-- iStage
      iD     !-- iDimension
    character ( 1 ) :: &
      StageNumber, &
      DimensionNumber

    associate ( nD  =>  3 )

    allocate ( RS % StageDimension ( nS, nD ) )
    do iS  =  1, nS
      do iD  =  1, nD
        write ( StageNumber, fmt = '(i1.1)' ) iS
        write ( DimensionNumber, fmt = '(i1.1)' ) iD
        allocate ( RS % StageDimension ( iS, iD ) % Element )
        associate ( SDC  =>  RS % StageDimension ( iS, iD ) % Element )
        call SDC % Initialize &
               ( RS % Atlas, &
                 FieldOption = RS % Field, &
                 NameOption = trim ( RS % Name ) // '_' // StageNumber // '_' &
                              // DimensionNumber, &
                 DeviceMemoryOption = RS % DeviceMemory, &
                 DevicesCommunicateOption = RS % DevicesCommunicate, &
                 nFieldsOption = RS % nFields )
        call S % AddFieldSet ( SDC )
        end associate !-- SDC
      end do !-- iD
    end do !-- iS

    end associate !-- nD

    associate &
      ( RBS  =>  RS % Reconstruction_BS, &
        RFS  =>  RS % Reconstruction_FS, &
        RES  =>  RS % Reconstruction_ES )
    call RBS % SetStream ( S, nS )
    call RFS % SetStream ( S, nS )
    call RES % SetStream ( S, nS )
    end associate !-- RBS, etc.

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      FS

    call FS % FieldSetForm % Show ( )
    call FS % FluxSet % Show ( )
    call FS % EigenspeedSet % Show ( )
    call FS % Reconstruction_BS % Show ( )
    call FS % Reconstruction_FS % Show ( )
    call FS % Reconstruction_ES % Show ( )

  end subroutine Show_FS


  function Timer ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  RS % iTimer )

    if ( iT == 0 ) then
      TimerName  =  RS % Name
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer


  function TimerKernel ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDF ) :: &
      TimerName

    associate ( iT  =>  RS % iTimerKernel )

    if ( iT == 0 ) then
      TimerName  =  trim ( RS % Name ) // '_Kernel' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function TimerKernel


  subroutine Compute ( RS, iD, T_Option, iS_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iD  !-- iDimensions
    type ( TimerForm ), intent ( in ), pointer, optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option

    type ( TimerForm ), pointer :: &
      T_FS, &
      T_ES, &
      T_RBS, &
      T_RFS, &
      T_RES, &
      T_K

    call Show ( 'Computing ' // trim ( RS % Type ), RS % IGNORABILITY + 4 )
    call Show ( RS % Name, 'Name', RS % IGNORABILITY + 4 )

    associate &
      (  CS  =>  RS % CurrentSet, &
         FS  =>  RS % FluxSet, &
         ES  =>  RS % EigenspeedSet, &
        RBS  =>  RS % Reconstruction_BS, &
        RFS  =>  RS % Reconstruction_FS, &
        RES  =>  RS % Reconstruction_ES )

    if ( present ( T_Option ) ) then
      T_FS   =>   FS % Timer ( LevelOption = T_Option % Level + 1 )
      T_ES   =>   ES % Timer ( LevelOption = T_Option % Level + 1 )
      T_RBS  =>  RBS % Timer ( LevelOption = T_Option % Level + 1 )
      T_RFS  =>  RFS % Timer ( LevelOption = T_Option % Level + 1 )
      T_RES  =>  RES % Timer ( LevelOption = T_Option % Level + 1 )
      T_K    =>   RS % TimerKernel ( LevelOption = T_Option % Level + 1 )
    else
      T_FS   =>  null ( )
      T_ES   =>  null ( )
      T_RBS  =>  null ( )
      T_RFS  =>  null ( )
      T_RES  =>  null ( )
      T_K    =>  null ( )
    end if !-- T_Option

    if ( associated ( T_FS ) ) call T_FS % Start ( )
    call FS % Compute ( iD )
    if ( associated ( T_FS ) ) call T_FS % Stop ( )

    if ( associated ( T_ES ) ) call T_ES % Start ( )
    call ES % Compute ( iD )
    if ( associated ( T_ES ) ) call T_ES % Stop ( )

    if ( associated ( T_RBS ) ) call T_RBS % Start ( )
    call RBS % Compute ( iD )
    if ( associated ( T_RBS ) ) call T_RBS % Stop ( )

    if ( associated ( T_RFS ) ) call T_RFS % Start ( )
    call RFS % Compute ( iD )
    if ( associated ( T_RFS ) ) call T_RFS % Stop ( )

    if ( associated ( T_RES ) ) call T_RES % Start ( )
    call RES % Compute ( iD )
    if ( associated ( T_RES ) ) call T_RES % Stop ( )

    if ( associated ( T_K ) ) call T_K % Start ( )
    associate &
      ( RSV     =>  RS % Storage ( 1 ) % Value, &
        RBV_IL  =>  RBS % Output_IL % Storage ( 1 ) % Value, &
        RBV_IR  =>  RBS % Output_IR % Storage ( 1 ) % Value, &
        RFV_IL  =>  RFS % Output_IL % Storage ( 1 ) % Value, &
        RFV_IR  =>  RFS % Output_IR % Storage ( 1 ) % Value, &
        REV_IL  =>  RES % Output_IL % Storage ( 1 ) % Value, &
        REV_IR  =>  RES % Output_IR % Storage ( 1 ) % Value )
    associate &
      (  F_I   =>  RSV ( :, 1 : CS % nBalanced ), &
        AP_I   =>  RSV ( :, RS % ALPHA_PLUS_U ), &
        AM_I   =>  RSV ( :, RS % ALPHA_MINUS_U ), &
         F_IL  =>  RFV_IL ( :, : ), &
         F_IR  =>  RFV_IR ( :, : ), &
         U_IL  =>  RBV_IL ( :, : ), &
         U_IR  =>  RBV_IR ( :, : ), &
        EP_IL  =>  REV_IL ( :, ES % EIGENSPEED_FAST_PLUS_U ), &
        EP_IR  =>  REV_IR ( :, ES % EIGENSPEED_FAST_PLUS_U ), &
        EM_IL  =>  REV_IL ( :, ES % EIGENSPEED_FAST_MINUS_U ), &
        EM_IR  =>  REV_IR ( :, ES % EIGENSPEED_FAST_MINUS_U ) )

    call ComputeKernel &
           ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
             F_I, AP_I, AM_I, UseDeviceOption = RS % DeviceMemory )

    end associate !-- F_I, etc.
    end associate !-- RSV, etc.
    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- CS, etc.

    if ( allocated ( RS % StageDimension ) .and. present ( iS_Option ) ) then
      associate ( SDC  =>  RS % StageDimension ( iS_Option, iD ) % Element )
      call RS % Copy ( SDC )
      end associate !-- SDC
    end if

  end subroutine Compute


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS

    if ( allocated ( RS % StageDimension ) ) &
      deallocate ( RS % StageDimension )
    if ( allocated ( RS % Reconstruction_ES ) ) &
      deallocate ( RS % Reconstruction_ES )
    if ( allocated ( RS % Reconstruction_FS ) ) &
      deallocate ( RS % Reconstruction_FS )
    if ( allocated ( RS % Reconstruction_BS ) ) &
      deallocate ( RS % Reconstruction_BS )
    if ( allocated ( RS % EigenspeedSet ) ) &
      deallocate ( RS % EigenspeedSet )
    if ( allocated ( RS % FluxSet ) ) &
      deallocate ( RS % FluxSet )
    if ( allocated ( RS % BalancedSet ) ) &
      deallocate ( RS % BalancedSet )

    nullify ( RS % CurrentSet )

  end subroutine Finalize


end module RiemannSolver_HLL__Form
