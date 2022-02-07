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
      iFiducialDensityFlux = 0
    integer ( KDI ) :: &
      iTimer     = 0, &
      iTimer_CFP = 0, &
      iTimer_A   = 0, &
      iTimer_K   = 0
    integer ( KDI ), dimension ( : ), allocatable :: &
      iaAbundances, &
      iaAbundanceFluxes
    character ( LDL ) :: &
      ReconstructedSet = ''
    class ( FieldSetForm ), allocatable :: &
      PrimitiveSet, &
      BalancedSet, &
      CurrentSet_IL, CurrentSet_IR, &
      FluxSet, &
      FluxSet_IL, FluxSet_IR
    class ( CurrentSetForm ), pointer :: &
      CurrentSet => null ( )
    class ( EigenspeedSet_F_Form ), allocatable :: &
      EigenspeedSet, &
      EigenspeedSet_IL, EigenspeedSet_IR
    class ( ReconstructionForm ), allocatable :: &
      Reconstruction_PS, &
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
      Timer_CFP
    procedure, private, pass :: &
      Timer_A
    procedure, private, pass :: &
      Timer_K
    procedure, public, pass :: &
      Prepare
    procedure, public, pass :: &
      Compute
    final :: &
      Finalize
  end type RiemannSolver_HLL_Form

    ! private :: &
    !   ComputeWithReconstructedFluxes, &
    !   ComputeWithReconstructedPrimitive

      private :: &
        ComputeAlphaKernel, &
        ComputeFluxesKernel!, &
!        ComputeAbundancesKernel

    interface
      
      module subroutine ComputeAlphaKernel &
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
      end subroutine ComputeAlphaKernel

      module subroutine ComputeFluxesKernel &
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
      end subroutine ComputeFluxesKernel

      ! module subroutine ComputeAbundancesKernel &
      !          ( RSV, CS_IL, CS_IR, iaAF, iaA, iAP, iAM, iFDF, UseDeviceOption )
      !   use Basics
      !   implicit none
      !   real ( KDR ), dimension ( :, : ), intent ( inout ) :: &
      !     RSV     !-- RiemannSolver Value
      !   real ( KDR ), dimension ( :, : ), intent ( in ) :: &
      !     CS_IL, CS_IR
      !   integer ( KDI ), dimension ( : ), intent ( in ) :: &
      !     iaAF, &  !-- iaAbundanceFluxes
      !     iaA      !-- iaAbundances
      !   integer ( KDI ), intent ( in ) :: &
      !     iAP, &  !-- iAlphaPlus
      !     iAM, &  !-- iAlphaMinus
      !     iFDF    !-- iFiducialDensityFlux
      !   logical ( KDL ), intent ( in ), optional :: &
      !     UseDeviceOption
      ! end subroutine ComputeAbundancesKernel

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
    
    Name  =  'RS_' // trim ( CS % Name )
    if ( present ( PrefixOption ) ) &
      Name  =  trim ( PrefixOption ) // '_' // trim ( CS % Name )

    RS % CurrentSet  =>  CS

    RS % ReconstructedSet  =  'PRIMITIVE'
    if ( present ( ReconstructedSetOption ) ) &
      RS % ReconstructedSet  =  ReconstructedSetOption
    call PROGRAM_HEADER % GetParameter &
      ( RS % ReconstructedSet, 'ReconstructedSet' )

    select case ( trim ( RS % ReconstructedSet ) )
    case ( 'FLUXES' )

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
      call FS % Initialize &
             ( CS % Atlas, &
               FieldOption = CS % Balanced, &
               NameOption = 'F_' // trim ( CS % Name ), &
               DeviceMemoryOption = CS % DeviceMemory, &
               DevicesCommunicateOption = CS % DevicesCommunicate, &
               nFieldsOption = CS % nBalanced, &
               IgnorabilityOption = CS % IGNORABILITY + 1 )               
      call ES % Initialize ( CS, CS )

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

      end associate !-- RBS, etc.
      end associate !-- BS, etc.

    case ( 'PRIMITIVE' )

      allocate ( RS % PrimitiveSet )
      associate ( PS  =>  RS % PrimitiveSet )
      call PS % Initialize &
             ( CS, CS % iaPrimitive, &
               NameOption = 'P_' // trim ( CS % Name ), &
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
        ( RS % FluxSet_IL, &
          RS % FluxSet_IR, &
          RS % EigenspeedSet_IL, &
          RS % EigenspeedSet_IR )
      associate &
        ( FS_IL  =>  RS % FluxSet_IL, &
          FS_IR  =>  RS % FluxSet_IR, &
          ES_IL  =>  RS % EigenspeedSet_IL, &
          ES_IR  =>  RS % EigenspeedSet_IR )
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

      end associate !-- FS_IL, etc.
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
             AssociateFieldsOption = .false., &
             nFieldsOption = nFields, &
             IgnorabilityOption = CS % IGNORABILITY )

    end associate !-- nB

  end subroutine InitializeAllocate_RS


  subroutine SetStream ( RS, S, nS )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    class ( StreamForm ), intent ( inout ) :: &
      S
    integer ( KDI ), intent ( in ) :: &
      nS  !-- nStages

    integer ( KDI ) :: &
      iC, &  !-- iChart
      iS, &  !-- iStage
      iD, &  !-- iDimension
      nD     !-- nDimensions
    character ( 1 ) :: &
      StageNumber, &
      DimensionNumber

    associate ( A  =>  RS % Atlas )
    nD  =  A % Chart ( 1 ) % Element % nDimensions
    do iC  =  2, A % nCharts
      nD  =  max ( nD, A % Chart ( iC ) % Element % nDimensions )
    end do

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
                 nFieldsOption = RS % nFields, &
                 IgnorabilityOption = RS % IGNORABILITY + 1 )
        call S % AddFieldSet ( SDC )
        end associate !-- SDC
      end do !-- iD
    end do !-- iS

    end associate !-- A

    select case ( trim ( RS % ReconstructedSet ) )
    case ( 'FLUXES' )
      associate &
        ( RBS  =>  RS % Reconstruction_BS, &
          RFS  =>  RS % Reconstruction_FS, &
          RES  =>  RS % Reconstruction_ES )
      call RBS % SetStream ( S, nS )
      call RFS % SetStream ( S, nS )
      call RES % SetStream ( S, nS )
      end associate !-- RBS, etc.
    case ( 'PRIMITIVE' )
      associate ( RPS  =>  RS % Reconstruction_PS )
      call RPS % SetStream ( S, nS )
      end associate !-- RPS
    case default
      call Show ( 'ReconstructedSet not recognized', CONSOLE % ERROR )
      call Show ( RS % ReconstructedSet, 'ReconstructedSet', CONSOLE % ERROR )
      call Show ( 'RiemannSolver_HLL__Form', 'module', CONSOLE % ERROR )
      call Show ( 'SetStream', 'subroutine', CONSOLE % ERROR )
      call PROGRAM_HEADER % Abort ( )
    end select !-- ReconstructedSet

  end subroutine SetStream


  subroutine Show_FS ( FS )

    class ( RiemannSolver_HLL_Form ), intent ( in ) :: &
      FS

    call FS % FieldSetForm % Show ( )
    call Show ( FS % ReconstructedSet, 'ReconstructedSet', FS % IGNORABILITY )

    select case ( trim ( FS % ReconstructedSet ) )
    case ( 'FLUXES' )
      call FS % FluxSet % Show ( )
      call FS % EigenspeedSet % Show ( )
      call FS % Reconstruction_BS % Show ( )
      call FS % Reconstruction_FS % Show ( )
      call FS % Reconstruction_ES % Show ( )
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


  function Timer ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
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


  function Timer_CFP ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  RS % iTimer_CFP )

    if ( iT == 0 ) then
      TimerName  =  trim ( RS % Name ) // '_CFP' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_CFP


  function Timer_A ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  RS % iTimer_K )

    if ( iT == 0 ) then
      TimerName  =  trim ( RS % Name ) // '_A' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_A


  function Timer_K ( RS, LevelOption ) result ( T )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ), optional :: &
      LevelOption
    type ( TimerForm ), pointer :: &
      T

    character ( LDL ) :: &
      TimerName

    associate ( iT  =>  RS % iTimer_K )

    if ( iT == 0 ) then
      TimerName  =  trim ( RS % Name ) // '_K' 
      if ( present ( LevelOption ) ) then
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, LevelOption )
      else
        call PROGRAM_HEADER % AddTimer ( TimerName, iT, Level = 1 )
      end if
    end if

    T  =>  PROGRAM_HEADER % TimerPointer ( iT )

    end associate !-- iT

  end function Timer_K


  ! subroutine Compute ( RS, iC, iD, T_Option, iS_Option )

  !   class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
  !     RS
  !   integer ( KDI ), intent ( in ) :: &
  !     iC, &  !-- iChart
  !     iD     !-- iDimensions
  !   type ( TimerForm ), intent ( in ), optional :: &
  !     T_Option
  !   integer ( KDI ), intent ( in ), optional :: &
  !     iS_Option  !-- iStage_Option

  !   call Show ( 'Computing ' // trim ( RS % Type ), RS % IGNORABILITY + 3 )
  !   call Show ( RS % Name, 'Name', RS % IGNORABILITY + 3 )

  !   select case ( trim ( RS % ReconstructedSet ) )
  !   case ( 'FLUXES' )
  !     call ComputeWithReconstructedFluxes ( RS, iC, iD, T_Option, iS_Option )
  !   case ( 'PRIMITIVE' )
  !     call ComputeWithReconstructedPrimitive ( RS, iC, iD, T_Option, iS_Option )
  !   case default
  !     call Show ( 'ReconstructedSet not recognized', CONSOLE % ERROR )
  !     call Show ( RS % ReconstructedSet, 'ReconstructedSet', CONSOLE % ERROR )
  !     call Show ( 'RiemannSolver_HLL__Form', 'module', CONSOLE % ERROR )
  !     call Show ( 'Compute', 'subroutine', CONSOLE % ERROR )
  !     call PROGRAM_HEADER % Abort ( )
  !   end select !-- ReconstructedSet

  !   if ( allocated ( RS % StageDimension ) .and. present ( iS_Option ) ) then
  !     associate ( SDC  =>  RS % StageDimension ( iS_Option, iD ) % Element )
  !     call RS % Copy ( SDC )
  !     end associate !-- SDC
  !   end if

  ! end subroutine Compute


  subroutine Prepare ( RS, iC, iD, T_Option, iS_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
      RS
    integer ( KDI ), intent ( in ) :: &
      iC, &  !-- iChart
      iD     !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option  !-- iStage_Option

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
      T_RPS  =>  RPS    % Timer     ( LevelOption = T_Option % Level + 1 )
      T_CFP  =>   RS    % Timer_CFP ( LevelOption = T_Option % Level + 1 )
      T_ES   =>   ES_IL % Timer     ( LevelOption = T_Option % Level + 1 )
      T_A    =>   RS    % Timer_A   ( LevelOption = T_Option % Level + 1 )
    else
      T_RPS  =>  null ( )
      T_CFP  =>  null ( )
      T_ES   =>  null ( )
      T_A    =>  null ( )
    end if !-- T_Option

    if ( associated ( T_RPS ) ) call T_RPS % Start ( )
    call RPS % Compute ( iC, iD, iS_Option )
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
    
    call ComputeAlphaKernel &
           ( EP_IL, EP_IR, EM_IL, EM_IR, &
             RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, RSV, &
             UseDeviceOption = RS % DeviceMemory )

    end associate !-- RSV, etc.
    end associate !-- RSS, etc.
    if ( associated ( T_A ) ) call T_A % Stop ( )

    end associate !-- CS, etc.

  end subroutine Prepare


  subroutine Compute ( RS, DC, iC, iD, T_Option, iS_Option )

    class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
         RS
    class ( DivergenceContribution_CS_Form ), intent ( inout ) :: &
      DC
    integer ( KDI ), intent ( in ) :: &
      iC, &   !-- iChart
      iD      !-- iDimensions
    type ( TimerForm ), intent ( in ), optional :: &
      T_Option
    integer ( KDI ), intent ( in ), optional :: &
      iS_Option  !-- iStage_Option

    integer ( KDI ) :: &
      iF 
    integer ( KDI ), dimension ( RS % CurrentSet % nBalanced ) :: &
      iaFluxes
    type ( TimerForm ), pointer :: &
      T_F, &
      T_K
    
    associate &
      (  CS     =>  RS % CurrentSet, &
         CS_IL  =>  RS % CurrentSet_IL, &
         CS_IR  =>  RS % CurrentSet_IR, &
         FS_IL  =>  RS % FluxSet_IL, &
         FS_IR  =>  RS % FluxSet_IR, &
         ES_IL  =>  RS % EigenspeedSet_IL, &
         ES_IR  =>  RS % EigenspeedSet_IR, &
        RPS     =>  RS % Reconstruction_PS )

    if ( present ( T_Option ) ) then
      T_F   =>   DC % Timer   ( LevelOption = T_Option % Level + 1 )
      T_K   =>   RS % Timer_K ( LevelOption = T_Option % Level + 1 )
    else
      T_F   =>  null ( )
      T_K    =>  null ( )
    end if !-- T_Option

    if ( associated ( T_F ) ) call T_F % Start ( )
    call DC % ComputeFluxes ( FS_IL, CS_IL, iC, iD )
    call DC % ComputeFluxes ( FS_IR, CS_IR, iC, iD )
    if ( associated ( T_F ) ) call T_F % Stop ( )

    if ( associated ( T_K ) ) call T_K % Start ( )
    associate &
      ( RSS     =>  RS % Storage ( iC ), &
        CSS_IL  =>  RS % CurrentSet_IL % Storage ( iC ), &
        CSS_IR  =>  RS % CurrentSet_IR % Storage ( iC ), &
        RFS_IL  =>  RS % FluxSet_IL % Storage ( iC ), &
        RFS_IR  =>  RS % FluxSet_IR % Storage ( iC ) )
    associate &
      (   RSV  =>  RSS % Value, &
         F_IL  =>  RFS_IL % Value, &
         F_IR  =>  RFS_IR % Value, &
         U_IL  =>  CSS_IL % Value, &
         U_IR  =>  CSS_IR % Value )
    
    iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    call RFS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
    call RFS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    
    call ComputeFluxesKernel &
           ( RSV, F_IL, F_IR, U_IL, U_IR, CS % iaBalanced, iaFluxes, &
             RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, &
             UseDeviceOption = RS % DeviceMemory )

    call RFS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call RFS_IL % ReassociateHost ( AssociateVariablesOption = .true. )
    call CSS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
    call CSS_IL % ReassociateHost ( AssociateVariablesOption = .true. )

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
    if ( allocated ( RS % Reconstruction_PS ) ) &
      deallocate ( RS % Reconstruction_PS )
    if ( allocated ( RS % EigenspeedSet_IR ) ) &
      deallocate ( RS % EigenspeedSet_IR )
    if ( allocated ( RS % EigenspeedSet_IL ) ) &
      deallocate ( RS % EigenspeedSet_IL )
    if ( allocated ( RS % EigenspeedSet ) ) &
      deallocate ( RS % EigenspeedSet )
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
    if ( allocated ( RS % BalancedSet ) ) &
      deallocate ( RS % BalancedSet )
    if ( allocated ( RS % PrimitiveSet ) ) &
      deallocate ( RS % PrimitiveSet )
    if ( allocated ( RS % iaAbundanceFluxes ) ) &
      deallocate ( RS % iaAbundanceFluxes )
    if ( allocated ( RS % iaAbundances ) ) &
      deallocate ( RS % iaAbundances )

    nullify ( RS % CurrentSet )

  end subroutine Finalize


  ! subroutine ComputeWithReconstructedFluxes ( RS, iC, iD, T_Option, iS_Option )

  !   class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
  !     RS
  !   integer ( KDI ), intent ( in ) :: &
  !     iC, &  !-- iChart
  !     iD     !-- iDimensions
  !   type ( TimerForm ), intent ( in ), optional :: &
  !     T_Option
  !   integer ( KDI ), intent ( in ), optional :: &
  !     iS_Option  !-- iStage_Option

  !   integer ( KDI ) :: &
  !     iF
  !   integer ( KDI ), dimension ( 1 : RS % CurrentSet % nBalanced ) :: &
  !     iaFluxes 
  !   type ( TimerForm ), pointer :: &
  !     T_FS, &
  !     T_ES, &
  !     T_RBS, &
  !     T_RFS, &
  !     T_RES, &
  !     T_K
    
  !   associate &
  !     (  CS  =>  RS % CurrentSet, &
  !        FS  =>  RS % FluxSet, &
  !        ES  =>  RS % EigenspeedSet, &
  !       RBS  =>  RS % Reconstruction_BS, &
  !       RFS  =>  RS % Reconstruction_FS, &
  !       RES  =>  RS % Reconstruction_ES )

  !   if ( present ( T_Option ) ) then
  !     T_FS   =>   FS % Timer ( LevelOption = T_Option % Level + 1 )
  !     T_ES   =>   ES % Timer ( LevelOption = T_Option % Level + 1 )
  !     T_RBS  =>  RBS % Timer ( LevelOption = T_Option % Level + 1 )
  !     T_RFS  =>  RFS % Timer ( LevelOption = T_Option % Level + 1 )
  !     T_RES  =>  RES % Timer ( LevelOption = T_Option % Level + 1 )
  !     T_K    =>   RS % Timer_K ( LevelOption = T_Option % Level + 1 )
  !   else
  !     T_FS   =>  null ( )
  !     T_ES   =>  null ( )
  !     T_RBS  =>  null ( )
  !     T_RFS  =>  null ( )
  !     T_RES  =>  null ( )
  !     T_K    =>  null ( )
  !   end if !-- T_Option

  !   if ( associated ( T_FS ) ) call T_FS % Start ( )
  !   call FS % Compute ( iC, iD )
  !   if ( associated ( T_FS ) ) call T_FS % Stop ( )

  !   if ( associated ( T_ES ) ) call T_ES % Start ( )
  !   call ES % Compute ( iC, iD )
  !   if ( associated ( T_ES ) ) call T_ES % Stop ( )

  !   if ( associated ( T_RBS ) ) call T_RBS % Start ( )
  !   call RBS % Compute ( iC, iD, iS_Option )
  !   if ( associated ( T_RBS ) ) call T_RBS % Stop ( )

  !   if ( associated ( T_RFS ) ) call T_RFS % Start ( )
  !   call RFS % Compute ( iC, iD, iS_Option )
  !   if ( associated ( T_RFS ) ) call T_RFS % Stop ( )

  !   if ( associated ( T_RES ) ) call T_RES % Start ( )
  !   call RES % Compute ( iC, iD, iS_Option )
  !   if ( associated ( T_RES ) ) call T_RES % Stop ( )

  !   if ( associated ( T_K ) ) call T_K % Start ( )
  !   associate &
  !     ( RSS     =>  RS % Storage ( iC ), &
  !       RBS_IL  =>  RBS % Output_IL % Storage ( iC ), &
  !       RBS_IR  =>  RBS % Output_IR % Storage ( iC ), &
  !       RFS_IL  =>  RFS % Output_IL % Storage ( iC ), &
  !       RFS_IR  =>  RFS % Output_IR % Storage ( iC ), &
  !       RES_IL  =>  RES % Output_IL % Storage ( iC ), &
  !       RES_IR  =>  RES % Output_IR % Storage ( iC ) )
  !   associate &
  !     (   RSV  =>  RSS % Value, &
  !        F_IL  =>  RFS_IL % Value, &
  !        F_IR  =>  RFS_IR % Value, &
  !        U_IL  =>  RBS_IL % Value, &
  !        U_IR  =>  RBS_IR % Value, &
  !       EP_IL  =>  RES_IL % Value ( :, ES % EIGENSPEED_FAST_PLUS_U ), &
  !       EP_IR  =>  RES_IR % Value ( :, ES % EIGENSPEED_FAST_PLUS_U ), &
  !       EM_IL  =>  RES_IL % Value ( :, ES % EIGENSPEED_FAST_MINUS_U ), &
  !       EM_IR  =>  RES_IR % Value ( :, ES % EIGENSPEED_FAST_MINUS_U ) )
    
  !   iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
  !   call RBS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RBS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RFS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RFS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
        
  !   call ComputeKernel &
  !          ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
  !            RS % FluxSet % iaSelected, iaFluxes, RS % ALPHA_PLUS_U, &
  !            RS % ALPHA_MINUS_U, RSV, UseDeviceOption = RS % DeviceMemory )
    
  !   call RFS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RFS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RBS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RBS_IL % ReassociateHost ( AssociateVariablesOption = .false. )

  !   end associate !-- F_I, etc.
  !   end associate !-- RSV, etc.
  !   if ( associated ( T_K ) ) call T_K % Stop ( )

  !   end associate !-- CS, etc.

  !   if ( allocated ( RS % StageDimension ) .and. present ( iS_Option ) ) then
  !     associate ( SDC  =>  RS % StageDimension ( iS_Option, iD ) % Element )
  !     call RS % Copy ( SDC )
  !     end associate !-- SDC
  !   end if

  ! end subroutine ComputeWithReconstructedFluxes


  ! subroutine ComputeWithReconstructedPrimitive &
  !              ( RS, iC, iD, T_Option, iS_Option )

  !   class ( RiemannSolver_HLL_Form ), intent ( inout ) :: &
  !     RS
  !   integer ( KDI ), intent ( in ) :: &
  !     iC, &  !-- iChart
  !     iD     !-- iDimensions
  !   type ( TimerForm ), intent ( in ), optional :: &
  !     T_Option
  !   integer ( KDI ), intent ( in ), optional :: &
  !     iS_Option  !-- iStage_Option

  !   integer ( KDI ) :: &
  !     iF 
  !   integer ( KDI ), dimension ( RS % CurrentSet % nBalanced ) :: &
  !     iaFluxes
  !   type ( TimerForm ), pointer :: &
  !     T_RPS, &
  !     T_CFP, &
  !     T_FS, &
  !     T_ES, &
  !     T_K
    
  !   associate &
  !     (  CS     =>  RS % CurrentSet, &
  !        CS_IL  =>  RS % CurrentSet_IL, &
  !        CS_IR  =>  RS % CurrentSet_IR, &
  !        FS_IL  =>  RS % FluxSet_IL, &
  !        FS_IR  =>  RS % FluxSet_IR, &
  !        ES_IL  =>  RS % EigenspeedSet_IL, &
  !        ES_IR  =>  RS % EigenspeedSet_IR, &
  !       RPS     =>  RS % Reconstruction_PS )

  !   if ( present ( T_Option ) ) then
  !     T_RPS  =>  RPS    % Timer     ( LevelOption = T_Option % Level + 1 )
  !     T_CFP  =>   RS    % Timer_CFP ( LevelOption = T_Option % Level + 1 )
  !     T_FS   =>   FS_IL % Timer     ( LevelOption = T_Option % Level + 1 )
  !     T_ES   =>   ES_IL % Timer     ( LevelOption = T_Option % Level + 1 )
  !     T_K    =>   RS    % Timer_K   ( LevelOption = T_Option % Level + 1 )
  !   else
  !     T_RPS  =>  null ( )
  !     T_CFP  =>  null ( )
  !     T_FS   =>  null ( )
  !     T_ES   =>  null ( )
  !     T_K    =>  null ( )
  !   end if !-- T_Option

  !   if ( associated ( T_RPS ) ) call T_RPS % Start ( )
  !   call RPS % Compute ( iC, iD, iS_Option )
  !   if ( associated ( T_RPS ) ) call T_RPS % Stop ( )

  !   if ( associated ( T_CFP ) ) call T_CFP % Start ( )
  !   call CS % ComputeFromPrimitive ( CS_IL )
  !   call CS % ComputeFromPrimitive ( CS_IR )
  !   if ( associated ( T_CFP ) ) call T_CFP % Stop ( )

  !   if ( associated ( T_FS ) ) call T_FS % Start ( )
  !   call FS_IL % Compute ( iC, iD )
  !   call FS_IR % Compute ( iC, iD )
  !   if ( associated ( T_FS ) ) call T_FS % Stop ( )

  !   if ( associated ( T_ES ) ) call T_ES % Start ( )
  !   call ES_IL % Compute ( iC, iD )
  !   call ES_IR % Compute ( iC, iD )
  !   if ( associated ( T_ES ) ) call T_ES % Stop ( )

  !   if ( associated ( T_K ) ) call T_K % Start ( )
  !   associate &
  !     ( RSS     =>  RS % Storage ( iC ), &
  !       CSS_IL  =>  RS % CurrentSet_IL % Storage ( iC ), &
  !       CSS_IR  =>  RS % CurrentSet_IR % Storage ( iC ), &
  !       RFS_IL  =>  RS % FluxSet_IL % Storage ( iC ), &
  !       RFS_IR  =>  RS % FluxSet_IR % Storage ( iC ), &
  !       RES_IL  =>  RS % EigenspeedSet_IL % Storage ( iC ), &
  !       RES_IR  =>  RS % EigenspeedSet_IR % Storage ( iC ) )
  !   associate &
  !     (   RSV  =>  RSS % Value, &
  !        F_IL  =>  RFS_IL % Value, &
  !        F_IR  =>  RFS_IR % Value, &
  !        U_IL  =>  CSS_IL % Value, &
  !        U_IR  =>  CSS_IR % Value, &
  !       CS_IL  =>  CSS_IL % Value, &
  !       CS_IR  =>  CSS_IR % Value, &
  !       EP_IL  =>  RES_IL % Value ( :, ES_IL % EIGENSPEED_FAST_PLUS_U ), &
  !       EP_IR  =>  RES_IR % Value ( :, ES_IR % EIGENSPEED_FAST_PLUS_U ), &
  !       EM_IL  =>  RES_IL % Value ( :, ES_IL % EIGENSPEED_FAST_MINUS_U ), &
  !       EM_IR  =>  RES_IR % Value ( :, ES_IR % EIGENSPEED_FAST_MINUS_U ) )
    
  !   iaFluxes = [ ( iF, iF = 1, CS % nBalanced ) ]
    
  !   call CSS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call CSS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RFS_IL % ReassociateHost ( AssociateVariablesOption = .false. )
  !   call RFS_IR % ReassociateHost ( AssociateVariablesOption = .false. )
    
  !   call ComputeKernel &
  !          ( F_IL, F_IR, U_IL, U_IR, EP_IL, EP_IR, EM_IL, EM_IR, &
  !            CS % iaBalanced, iaFluxes, RS % ALPHA_PLUS_U, &
  !            RS % ALPHA_MINUS_U, RSV, UseDeviceOption = RS % DeviceMemory )

  !   ! if ( allocated ( RS % iaAbundances ) ) &
  !   !   call ComputeAbundancesKernel &
  !   !            ( RSV, CS_IL, CS_IR, &
  !   !              RS % iaAbundanceFluxes, RS % iaAbundances, &
  !   !              RS % ALPHA_PLUS_U, RS % ALPHA_MINUS_U, &
  !   !              RS % iFiducialDensityFlux, &
  !   !              UseDeviceOption = RS % DeviceMemory )

  !   call RFS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
  !   call RFS_IL % ReassociateHost ( AssociateVariablesOption = .true. )
  !   call CSS_IR % ReassociateHost ( AssociateVariablesOption = .true. )
  !   call CSS_IL % ReassociateHost ( AssociateVariablesOption = .true. )

  !   end associate !-- F_I, etc.
  !   end associate !-- RSV, etc.
  !   if ( associated ( T_K ) ) call T_K % Stop ( )

  !   end associate !-- CS, etc.

  !   if ( allocated ( RS % StageDimension ) .and. present ( iS_Option ) ) then
  !     associate ( SDC  =>  RS % StageDimension ( iS_Option, iD ) % Element )
  !     call RS % Copy ( SDC )
  !     end associate !-- SDC
  !   end if
    
  ! end subroutine ComputeWithReconstructedPrimitive


end module RiemannSolver_HLL__Form
