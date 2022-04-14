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
      ALPHA_MINUS_U   = 0
    integer ( KDI ) :: &
      iTimer_P     = 0, &  !-- Prepare
      iTimer_C     = 0, &  !-- Compute
      iTimer_R     = 0, &  !-- Reconstruction
      iTimer_CFP   = 0, &  !-- ComputeFromPrimitive
      iTimer_E     = 0, &  !-- Eigenspeeds
      iTimer_A     = 0, &  !-- Alpha
      iTimer_F     = 0, &  !-- Flux
      iTimer_K_HLL = 0     !-- Kernel_HLL
    character ( LDL ) :: &
      ReconstructedSet = ''
    class ( FieldSetForm ), allocatable :: &
      PrimitiveSet, &
      CurrentSet_IL, CurrentSet_IR, &
      FluxSet, &
      FluxSet_IL, FluxSet_IR
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
    class ( EigenspeedSet_F_Form ), allocatable :: &
      EigenspeedSet_IL, EigenspeedSet_IR
    class ( ReconstructionForm ), allocatable :: &
      Reconstruction_PS
  contains
    procedure, private, pass :: &
      InitializeAllocate_RS
    generic, public :: &
      Initialize => InitializeAllocate_RS
    procedure, public, pass :: &
      Show => Show_FS
    procedure, public, pass :: &
      Timer_P  !-- Prepare
    procedure, public, pass :: &
      Timer_C  !-- Compute
    procedure, public, pass :: &
      Prepare
    procedure, public, pass :: &
      ComputeFlux
    procedure, public, pass :: &
      ComputeDiffusion
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLL_Form

    private :: &
      PrepareKernel, &
      ComputeFluxKernel, &
      ComputeKernel

    interface
      
      module subroutine PrepareKernel &
               ( EP_IL, EP_IR, EM_IL, EM_IR, iAP, iAM, RSV, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( : ), intent ( in ) :: &
          EP_IL, EP_IR, &
          EM_IL, EM_IR
        integer ( KDI ), intent ( in ) :: &
          iAP, &  !-- iAlphaPlus
          iAM     !-- iAlphaMinus
        real ( KDR ), dimension ( :, : ), intent ( out ) :: &
          RSV     !-- RiemannSolver Value
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine PrepareKernel

      module subroutine ComputeFluxKernel &
               ( RSV, F_IL, F_IR, iaFluxes, iAP, iAM, UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          RSV     !-- RiemannSolver Value
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          F_IL, F_IR
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaFluxes
        integer ( KDI ), intent ( in ) :: &
          iAP, &  !-- iAlphaPlus
          iAM     !-- iAlphaMinus
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeFluxKernel

      module subroutine ComputeDiffusionKernel &
               ( RSV, U_IL, U_IR, iaBalanced, iaFluxes, iAP, iAM, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          RSV     !-- RiemannSolver Value
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          U_IL, U_IR
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaBalanced, &
          iaFluxes
        integer ( KDI ), intent ( in ) :: &
          iAP, &  !-- iAlphaPlus
          iAM     !-- iAlphaMinus
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeDiffusionKernel

      module subroutine ComputeKernel &
               ( RSV, F_IL, F_IR, U_IL, U_IR, iaBalanced, iaFluxes, iAP, iAM, &
                 UseDeviceOption )
        use Basics
        implicit none
        real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
          RSV     !-- RiemannSolver Value
        real ( KDR ), dimension ( :, : ), intent ( in ) :: &
          F_IL, F_IR, &
          U_IL, U_IR
        integer ( KDI ), dimension ( : ), intent ( in ) :: &
          iaBalanced, &
          iaFluxes
        integer ( KDI ), intent ( in ) :: &
          iAP, &  !-- iAlphaPlus
          iAM     !-- iAlphaMinus
        logical ( KDL ), intent ( in ), optional :: &
          UseDeviceOption
      end subroutine ComputeKernel

    end interface


contains


  subroutine InitializeAllocate_RS &
               ( RS, CS, FieldOption, ReconstructedSetOption, PrefixOption, &
                 nFieldsOption )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    class ( CurrentSetForm ), intent ( in ), target :: &
      CS
    character ( * ), dimension ( : ), intent ( in ), optional :: &
      FieldOption
    character ( * ), intent ( in ), optional :: &
      ReconstructedSetOption, &
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
    
    Name  =  trim ( CS % Name ) // '_RmnnSlvr'
    if ( present ( PrefixOption ) ) &
      Name  =  trim ( PrefixOption ) // '_' // trim ( CS % Name )

    RS % CurrentSet  =>  CS

    RS % ReconstructedSet  =  'PRIMITIVE'
    if ( present ( ReconstructedSetOption ) ) &
      RS % ReconstructedSet  =  ReconstructedSetOption
    call PROGRAM_HEADER % GetParameter &
      ( RS % ReconstructedSet, 'ReconstructedSet' )

    select case ( trim ( RS % ReconstructedSet ) )
    case ( 'PRIMITIVE' )

      allocate ( RS % PrimitiveSet )
      associate ( PS  =>  RS % PrimitiveSet )
      call PS % Initialize &
             ( CS, CS % iaPrimitive, &
               NameOption = trim ( CS % Name ) // '_Prmtv', &
               IgnorabilityOption = CS % IGNORABILITY + 1 )

      allocate &
        ( RS % CurrentSet_IL, &
          RS % CurrentSet_IR )
      associate &
        ( CS_IL  =>  RS % CurrentSet_IL, &
          CS_IR  =>  RS % CurrentSet_IR )
      call CS_IL % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Field, &
               NameOption = trim ( CS % Name ) // '_IL', &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               UnitOption = CS % Unit, &
               nFieldsOption = CS % nFields, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )
      call CS_IR % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Field, &
               NameOption = trim ( CS % Name ) // '_IR', &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               UnitOption = CS % Unit, &
               nFieldsOption = CS % nFields, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )

      allocate ( RS % Reconstruction_PS )
      associate ( RPS  =>  RS % Reconstruction_PS )
      call RPS % Initialize &
             ( CS % Geometry, PS, CS_IL, CS_IR, CS % iaPrimitive )

      allocate &
        ( RS % FluxSet, &
          RS % FluxSet_IL, &
          RS % FluxSet_IR, &
          RS % EigenspeedSet_IL, &
          RS % EigenspeedSet_IR )
      associate &
        ( FS     =>  RS % FluxSet, &
          FS_IL  =>  RS % FluxSet_IL, &
          FS_IR  =>  RS % FluxSet_IR, &
          ES_IL  =>  RS % EigenspeedSet_IL, &
          ES_IR  =>  RS % EigenspeedSet_IR )
      call FS % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Balanced, &
               NameOption = 'F_' // trim ( CS % Name ), &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               nFieldsOption = CS % nBalanced, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )               
      call FS_IL % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Balanced, &
               NameOption = 'F_IL_' // trim ( CS % Name ), &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               nFieldsOption = CS % nBalanced, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )               
      call FS_IR % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Balanced, &
               NameOption = 'F_IR_' // trim ( CS % Name ), &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               nFieldsOption = CS % nBalanced, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )               
      call ES_IL % Initialize ( CS, CS_IL )
      call ES_IR % Initialize ( CS, CS_IR )

      end associate !-- FS, etc.
      end associate !-- RPS
      end associate !-- CS_IL, etc.
      end associate !-- PS

    case default
      call Show ( 'ReconstructedSet not recognized', CONSOLE % ERROR )
      call Show ( RS % ReconstructedSet, 'ReconstructedSet', CONSOLE % ERROR )
      call Show ( 'RiemannSolver_HLL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'InitializeAllocate_RS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- ReconstructedSet

    associate ( nB  =>  CS % nBalanced )

    !-- Field indices

    nFields  =  nB  +  RS % N_SOLVER_SPEEDS_HLL
    if ( present ( nFieldsOption ) )  &
      nFields  =  nFieldsOption

    RS % ALPHA_PLUS_U   =  nB  +  1
    RS % ALPHA_MINUS_U  =  nB  +  2

    !-- Field names

    if ( present ( FieldOption ) ) then
      allocate ( Field, source = FieldOption )
    else
      allocate ( Field ( nFields ) )
    end if !-- FieldOption

    Field ( : nB )  =  CS % Balanced

    Field ( nB + 1 : nB + RS % N_SOLVER_SPEEDS_HLL ) &
      =  [ 'AlphaPlus_U ', &
           'AlphaMinus_U' ]
          
    !-- FieldSet

    call RS % FieldSetForm % Initialize &
           ( CS % Atlas, &
             FieldOption = Field, &
             NameOption = Name, &
             DeviceMemoryOption = CS % DeviceMemory, &
             DevicesCommunicateOption = CS % DevicesCommunicate, &
             AssociateFieldsOption = .false., &
             nFieldsOption = nFields, &
             IgnorabilityOption = CS % IGNORABILITY )

    end associate !-- nB

  end subroutine InitializeAllocate_RS


  subroutine Show_FS ( FS )

    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      FS

    call FS % FieldSetForm % Show ( )
    call Show ( FS % ReconstructedSet, 'ReconstructedSet', FS % IGNORABILITY )

    select case ( trim ( FS % ReconstructedSet ) )
    case ( 'PRIMITIVE' )
      call FS % Reconstruction_PS % Show ( )
      call FS % FluxSet_IL % Show ( )
      call FS % FluxSet_IR % Show ( )
      call FS % EigenspeedSet_IL % Show ( )
      call FS % EigenspeedSet_IR % Show ( )
    case default
      call Show ( 'ReconstructedSet not recognized', CONSOLE % ERROR )
      call Show ( FS % ReconstructedSet, 'ReconstructedSet', CONSOLE % ERROR )
      call Show ( 'RiemannSolver_HLL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'Show_FS', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- ReconstructedSet

  end subroutine Show_FS


  function Timer_P ( RS, Level ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = RS % iTimer_P, &
               Name = trim ( RS % Name ) // '_Prpr', &
               Level = Level )

  end function Timer_P


  function Timer_C ( RS, Level ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      Level
    type ( TimerForm ), pointer :: &
      T

    T  =>  PROGRAM_HEADER % Timer &
             ( Handle = RS % iTimer_C, &
               Name = trim ( RS % Name ) // '_Cmpt', &
               Level = Level )

  end function Timer_C


  subroutine Prepare ( RS, iC, iD, T_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    type ( TimerForm ), pointer :: &
      T_RPS, &
      T_CFP, &
      T_ES, &
      T_A
    
    associate &
      (  CS     =>  RS % CurrentSet, &
         CS_IL  =>  RS % CurrentSet_IL, &
         CS_IR  =>  RS % CurrentSet_IR, &
         ES_IL  =>  RS % EigenspeedSet_IL, &
         ES_IR  =>  RS % EigenspeedSet_IR, &
        RPS     =>  RS % Reconstruction_PS )

    if ( present ( T_Option ) ) then
      T_RPS  =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_R, &
                     Name = trim ( RS % Name ) // '_Rcnstrctn', &
                     Level = T_Option % Level + 1 )
      T_CFP  =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_CFP, &
                     Name = trim ( RS % Name ) // '_FrmPrmtv', &
                     Level = T_Option % Level + 1 )
      T_ES   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_E, &
                     Name = trim ( RS % Name ) // '_Egnspds', &
                     Level = T_Option % Level + 1 )
      T_A    =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_A, &
                     Name = trim ( RS % Name ) // '_Alpha', &
                     Level = T_Option % Level + 1 )
    else
      T_RPS  =>  null ( )
      T_CFP  =>  null ( )
      T_ES   =>  null ( )
      T_A    =>  null ( )
    end if !-- T_Option

    if ( associated ( T_RPS ) ) call T_RPS % Start ( )
    call RPS % Compute ( iC, iD )
    if ( associated ( T_RPS ) ) call T_RPS % Stop ( )

    if ( associated ( T_CFP ) ) call T_CFP % Start ( )
    call CS % ComputeFromPrimitive ( CS_IL )
    call CS % ComputeFromPrimitive ( CS_IR )
    if ( associated ( T_CFP ) ) call T_CFP % Stop ( )

    if ( associated ( T_ES ) ) call T_ES % Start ( )
    call ES_IL % Compute ( iC, iD )
    call ES_IR % Compute ( iC, iD )
    if ( associated ( T_ES ) ) call T_ES % Stop ( )

    if ( associated ( T_A ) ) call T_A % Start ( )
    associate &
      ( RSS     =>  RS % Storage ( iC ), &
        RES_IL  =>  RS % EigenspeedSet_IL % Storage ( iC ), &
        RES_IR  =>  RS % EigenspeedSet_IR % Storage ( iC ) )
    associate &
      ( RSV     =>  RSS % Value, &
         EP_IL  =>  RES_IL % Value ( :, ES_IL % EIGENSPEED_FAST_PLUS_U ), &
         EP_IR  =>  RES_IR % Value ( :, ES_IR % EIGENSPEED_FAST_PLUS_U ), &
         EM_IL  =>  RES_IL % Value ( :, ES_IL % EIGENSPEED_FAST_MINUS_U ), &
         EM_IR  =>  RES_IR % Value ( :, ES_IR % EIGENSPEED_FAST_MINUS_U ) )
    
    call PrepareKernel &
           ( EP_IL, EP_IR, EM_IL, EM_IR, &
             RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, RSV, &
             UseDeviceOption = RS % DeviceMemory )

    end associate !-- RSV, etc.
    end associate !-- RSS, etc.
    if ( associated ( T_A ) ) call T_A % Stop ( )

    end associate !-- CS, etc.

  end subroutine Prepare


  subroutine ComputeFlux ( RS, DP, iC, iD, T_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iF 
    integer ( KDI ), dimension ( RS % CurrentSet % nBalanced ) :: &
      iaFluxes
    type ( TimerForm ), pointer :: &
      T_F, &
      T_K
    
    call Show ( 'Computing ' // trim ( RS % Type ), RS % IGNORABILITY + 3 )
    call Show ( RS % Name, 'Name', RS % IGNORABILITY + 3 )

    associate &
      (  CS     =>  RS % CurrentSet, &
         CS_IL  =>  RS % CurrentSet_IL, &
         CS_IR  =>  RS % CurrentSet_IR, &
         FS     =>  RS % FluxSet, &
         FS_IL  =>  RS % FluxSet_IL, &
         FS_IR  =>  RS % FluxSet_IR, &
        RPS     =>  RS % Reconstruction_PS )

    if ( present ( T_Option ) ) then
      T_F   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_F, &
                     Name = trim ( RS % Name ) // '_Flx', &
                     Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_K_HLL, &
                     Name = trim ( RS % Name ) // '_Krnl_HLL', &
                     Level = T_Option % Level + 1 )
    else
      T_F   =>  null ( )
      T_K   =>  null ( )
    end if !-- T_Option

    if ( associated ( T_F ) ) call T_F % Start ( )
    call DP % ComputeFluxes ( FS,    CS,    iC, iD )
    call DP % ComputeFluxes ( FS_IL, CS_IL, iC, iD )
    call DP % ComputeFluxes ( FS_IR, CS_IR, iC, iD )
    if ( associated ( T_F ) ) call T_F % Stop ( )

    if ( associated ( T_K ) ) call T_K % Start ( )
    associate &
      ( RSS     =>  RS % Storage ( iC ), &
        RFS_IL  =>  RS % FluxSet_IL % Storage ( iC ), &
        RFS_IR  =>  RS % FluxSet_IR % Storage ( iC ) )
    associate &
      (   RSV  =>  RSS % Value, &
         F_IL  =>  RFS_IL % Value, &
         F_IR  =>  RFS_IR % Value )
    
    iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
    call RFS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call RFS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    
    call ComputeFluxKernel &
           ( RSV, F_IL, F_IR, iaFluxes, RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, &
             UseDeviceOption = RS % DeviceMemory )

    call RFS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call RFS_IL % ReassociateHost ( AssociateVariablesOption = .true. )

    end associate !-- F_I, etc.
    end associate !-- RSV, etc.
    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- CS, etc.

  end subroutine ComputeFlux


  subroutine ComputeDiffusion ( RS, iC, iD, T_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iF 
    integer ( KDI ), dimension ( RS % CurrentSet % nBalanced ) :: &
      iaFluxes
    type ( TimerForm ), pointer :: &
      T_F, &
      T_K
    
    call Show ( 'Computing ' // trim ( RS % Type ), RS % IGNORABILITY + 3 )
    call Show ( RS % Name, 'Name', RS % IGNORABILITY + 3 )

    associate &
      (  CS     =>  RS % CurrentSet, &
         CS_IL  =>  RS % CurrentSet_IL, &
         CS_IR  =>  RS % CurrentSet_IR, &
        RPS     =>  RS % Reconstruction_PS )

    if ( present ( T_Option ) ) then
      T_K   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_K_HLL, &
                     Name = trim ( RS % Name ) // '_Krnl_HLL', &
                     Level = T_Option % Level + 1 )
    else
      T_K   =>  null ( )
    end if !-- T_Option

    if ( associated ( T_K ) ) call T_K % Start ( )
    associate &
      ( RSS     =>  RS % Storage ( iC ), &
        CSS_IL  =>  RS % CurrentSet_IL % Storage ( iC ), &
        CSS_IR  =>  RS % CurrentSet_IR % Storage ( iC ) )
    associate &
      (   RSV  =>  RSS % Value, &
         U_IL  =>  CSS_IL % Value, &
         U_IR  =>  CSS_IR % Value )
    
    iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    
    call ComputeDiffusionKernel &
           ( RSV, U_IL, U_IR, CS % iaBalanced, iaFluxes, &
             RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, &
             UseDeviceOption = RS % DeviceMemory )

    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .true. )

    end associate !-- F_I, etc.
    end associate !-- RSV, etc.
    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- CS, etc.

  end subroutine ComputeDiffusion


  subroutine Compute ( RS, DP, iC, iD, T_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ), target :: &
      RS
    class ( DivergencePart_CS_Form ), intent ( inout ) :: &
      DP
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option

    integer ( KDI ) :: &
      iF 
    integer ( KDI ), dimension ( RS % CurrentSet % nBalanced ) :: &
      iaFluxes
    type ( TimerForm ), pointer :: &
      T_F, &
      T_K
    
    call Show ( 'Computing ' // trim ( RS % Type ), RS % IGNORABILITY + 3 )
    call Show ( RS % Name, 'Name', RS % IGNORABILITY + 3 )

    associate &
      (  CS     =>  RS % CurrentSet, &
         CS_IL  =>  RS % CurrentSet_IL, &
         CS_IR  =>  RS % CurrentSet_IR, &
         FS_IL  =>  RS % FluxSet_IL, &
         FS_IR  =>  RS % FluxSet_IR )

    if ( present ( T_Option ) ) then
      T_F   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_F, &
                     Name = trim ( RS % Name ) // '_Flx', &
                     Level = T_Option % Level + 1 )
      T_K   =>  PROGRAM_HEADER % Timer &
                   ( Handle = RS % iTimer_K_HLL, &
                     Name = trim ( RS % Name ) // '_Krnl_HLL', &
                     Level = T_Option % Level + 1 )
    else
      T_F   =>  null ( )
      T_K   =>  null ( )
    end if !-- T_Option

    if ( associated ( T_F ) ) call T_F % Start ( )
    call DP % ComputeFluxes ( FS_IL, CS_IL, iC, iD )
    call DP % ComputeFluxes ( FS_IR, CS_IR, iC, iD )
    if ( associated ( T_F ) ) call T_F % Stop ( )

    if ( associated ( T_K ) ) call T_K % Start ( )
    associate &
      ( RSS     =>  RS % Storage ( iC ), &
        CSS_IL  =>  RS % CurrentSet_IL % Storage ( iC ), &
        CSS_IR  =>  RS % CurrentSet_IR % Storage ( iC ), &
        FSS_IL  =>  RS % FluxSet_IL % Storage ( iC ), &
        FSS_IR  =>  RS % FluxSet_IR % Storage ( iC ) )
    associate &
      (   RSV  =>  RSS % Value, &
         F_IL  =>  FSS_IL % Value, &
         F_IR  =>  FSS_IR % Value, &
         U_IL  =>  CSS_IL % Value, &
         U_IR  =>  CSS_IR % Value )
    
    iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    call FSS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call FSS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    
    call ComputeKernel &
           ( RSV, F_IL, F_IR, U_IL, U_IR, CS % iaBalanced, iaFluxes, &
             RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, &
             UseDeviceOption = RS % DeviceMemory )

    call FSS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call FSS_IL % ReassociateHost ( AssociateVariablesOption = .true. )
    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .true. )

    end associate !-- RSV, etc.
    end associate !-- RSS, etc.
    if ( associated ( T_K ) ) call T_K % Stop ( )

    end associate !-- CS, etc.

  end subroutine Compute


  impure elemental subroutine Finalize ( RS )

    type ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS

    if ( allocated ( RS % Reconstruction_PS ) ) &
      deallocate ( RS % Reconstruction_PS )
    if ( allocated ( RS % EigenspeedSet_IR ) ) &
      deallocate ( RS % EigenspeedSet_IR )
    if ( allocated ( RS % EigenspeedSet_IL ) ) &
      deallocate ( RS % EigenspeedSet_IL )
    if ( allocated ( RS % FluxSet_IR ) ) &
      deallocate ( RS % FluxSet_IR )
    if ( allocated ( RS % FluxSet_IL ) ) &
      deallocate ( RS % FluxSet_IL )
    if ( allocated ( RS % FluxSet ) ) &
      deallocate ( RS % FluxSet )
    if ( allocated ( RS % CurrentSet_IR ) ) &
      deallocate ( RS % CurrentSet_IR )
    if ( allocated ( RS % CurrentSet_IL ) ) &
      deallocate ( RS % CurrentSet_IL )
    if ( allocated ( RS % PrimitiveSet ) ) &
      deallocate ( RS % PrimitiveSet )

    nullify ( RS % CurrentSet )

  end subroutine Finalize


end module RiemannSolver_HLL__Form
